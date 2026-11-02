from __future__ import annotations

import csv, hashlib, json, os, platform, shutil, subprocess, sys
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .analysis import rank_poses_by_energy
from .validation import prepare_native_ligand, run_native_redocking
from .contacts import detect_contacts
from .assay import correlate_with_ranking, extract_sdf_assay
from .core import (box_from_reference, clean_receptor_pdb, pdbqt_atom_types, prepare_ligand,
                   prepare_receptor, resolve_executable, run_command, split_and_minimize_sdf,
                   write_dpf, write_gpf)
from .figures import activity_scatter, pymol_binding_figures
from .intake import choose_reference_ligand, hetero_groups, protein_intake_summary, write_group, write_reference_candidates, receptor_completeness, write_receptor_completeness
from .errors import SIDockError
from .production import PRIMARY_INTERACTION_ENGINE, SIDOCK_INTERACTION_STATUS
from .model import DockingParameters, Project
from .stages import (run_pose_extraction_stage, run_cluster_stage, run_plip_stage,
                     run_figures_stage, run_finalization_stage)
from .dashboard import write_results, write_dashboard
from .compliance import evaluate as run_compliance
from .state import update_state, ensure_compounds, update_compound
from . import __version__
from .orchestrator import StageRunner, StageResult
from .executors import ProductionContext, build_production_executors, production_context_from_project
from .reference import (select_reference_instance, write_reference_selection,
                        load_local_ccd, load_local_mmcif, map_instance_to_ccd,
                        sha256_file, reference_config_fingerprint)


def log(msg: str) -> None:
    print(f"[SIDock Ubuntu] {msg}", flush=True)


def enforce_redocking_validation(result: dict[str, object], allow_failed: bool = False) -> None:
    if result.get('status') != 'PASS' and not allow_failed:
        raise RuntimeError('Native-ligand redocking validation failed; normal ligand docking is blocked. Use --allow-failed-redocking-validation only for explicit testing.')


def prepare_project(project_dir: Path, receptor: Path, sdf: Path, padding: float, ga_runs: int, jobs: int, reference_ligand: str | None = None, allow_incomplete_peripheral: bool = False, config: dict | None = None) -> Project:
    project=Project(str(project_dir.resolve()), str(receptor.resolve()), ligand_source=str(sdf.resolve()), docking=DockingParameters(ga_runs=ga_runs, parallel_jobs=jobs))
    project.ensure_layout()
    summary=protein_intake_summary(receptor)
    (project.path/'intake.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    candidates = write_reference_candidates(receptor, project.path/'validation/reference_ligand_candidates.csv')
    ref_cfg=(config or {}).get('reference',{})
    ccd_path=ref_cfg.get('ccd_path'); mmcif_path=ref_cfg.get('mmcif_path')
    ccd_sources={}
    for row in candidates:
        src=load_local_ccd(ccd_path,row.get('ligand_code')) if ccd_path else {'available':False,'parse_status':'MISSING','warnings':['CCD_SOURCE_NOT_CONFIGURED']}
        if src.get('available'):
            row.update(map_instance_to_ccd(row,src.get('chemistry')))
            row['ccd_source_hash']=src.get('source_hash')
        else:
            row.setdefault('ccd_available',False); row.setdefault('mapping_fraction',None); row.setdefault('selection_warning','CCD_UNAVAILABLE')
        ccd_sources[row.get('ligand_code')]=src
    mmcif_source=load_local_mmcif(mmcif_path) if mmcif_path else {'available':False,'parse_status':'MISSING','warnings':['MMCIF_SOURCE_NOT_CONFIGURED']}
    # The legacy intake CSV remains for compatibility; the canonical selection
    # record is chemistry-aware and independent of any docking outcome.
    try:
        selection = select_reference_instance(candidates, reference_ligand)
    except RuntimeError as exc:
        raise SIDockError('REFERENCE_LIGAND_UNRESOLVED', str(exc)) from exc
    write_reference_selection(project.path, candidates, override=reference_ligand,
                              input_hashes={'receptor': _sha256(receptor),
                                            'mmcif': mmcif_source.get('source_hash'),
                                            'ccd': sorted({v.get('source_hash') for v in ccd_sources.values() if v.get('source_hash')}),
                                            'reference_config': reference_config_fingerprint(ref_cfg)},
                              provenance={'ccd':ccd_sources,'mmcif':mmcif_source})
    selected = project_dir / 'validation' / 'reference_ligand_resolution.json'
    ref = None
    if selected.exists():
        key = json.loads(selected.read_text(encoding='utf-8')).get('selected_instance', '')
        for candidate in hetero_groups(receptor):
            if candidate.key == key:
                ref = candidate; break
    if reference_ligand:
        ref = next((g for g in hetero_groups(receptor) if g.key == reference_ligand), None)
        if ref is None: raise SIDockError('REFERENCE_LIGAND_UNRESOLVED', f'Explicit reference ligand not found: {reference_ligand}')
    if ref is None:
        selected_key=selection['selected'].get('key')
        ref=next((g for g in hetero_groups(receptor) if g.key==selected_key),None)
    ref = ref or choose_reference_ligand(receptor)
    if not ref:
        raise SIDockError('REFERENCE_LIGAND_UNRESOLVED', 'No suitable co-crystallized ligand detected.')
    completeness = receptor_completeness(receptor, ref)
    write_receptor_completeness(completeness, project.path/'qc/receptor_completeness.csv')
    bad_site=[r for r in completeness if r['classification']=='INCOMPLETE_BINDING_SITE']
    bad_peripheral=[r for r in completeness if r['classification']=='INCOMPLETE_PERIPHERAL']
    if bad_site: raise SIDockError('BINDING_SITE_RESIDUE_INCOMPLETE', 'Incomplete canonical residue within 5 Å of reference ligand.', {'residues':bad_site})
    if bad_peripheral and not allow_incomplete_peripheral: raise SIDockError('BINDING_SITE_RESIDUE_INCOMPLETE', 'Incomplete peripheral residue within 5–8 Å requires --allow-incomplete-peripheral.', {'residues':bad_peripheral})
    (project.path/'validation/reference_ligand_resolution.json').write_text(json.dumps({'selected_instance': ref.key, 'selection_method': 'explicit_cli' if reference_ligand else 'deterministic_candidate_ranking', 'candidates': candidates}, indent=2), encoding='utf-8')
    (project.path/'validation/receptor_completeness_policy.json').write_text(json.dumps({'allowed_bad_residues':[r for r in completeness if r['classification']=='INCOMPLETE_REMOTE' or (allow_incomplete_peripheral and r['classification']=='INCOMPLETE_PERIPHERAL')], 'allow_incomplete_peripheral':allow_incomplete_peripheral},indent=2),encoding='utf-8')
    ref_path=write_group(ref, project.path/'receptor/reference_ligand.pdb')
    project.grid=box_from_reference(ref_path,padding=padding)
    project.save(); return project


def _sha256(path: Path) -> str:
    digest=hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(1024*1024), b''): digest.update(block)
    return digest.hexdigest()


def _version(executable: str | None, flag: str='--version') -> str:
    if not executable: return 'not found'
    try: return subprocess.run([executable, flag], capture_output=True, text=True, timeout=10).stdout.strip().splitlines()[0]
    except Exception: return 'unknown'


def _write_provenance(project: Project, receptor: Path, sdf: Path, seed: int) -> dict[str, object]:
    """Write provenance early so failed validation runs remain auditable."""
    provenance={'sidock_version':__version__,'python_version':platform.python_version(),'rdkit_version':__import__('rdkit').__version__,'meeko_version':getattr(__import__('meeko'),'__version__','unknown'),'autodock_version':_version(resolve_executable('autodock4','AutoDock4')),'autogrid_version':_version(resolve_executable('autogrid4','AutoGrid4')),'pymol_version':_version(shutil.which('pymol'),'-h'),'operating_system':platform.platform(),'command_line':' '.join(sys.argv),'random_seed':seed,'utc_datetime':datetime.now(timezone.utc).isoformat(),'input_sha256':{'receptor':_sha256(receptor),'ligands':_sha256(sdf)}}
    (project.path/'provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    return provenance


def run_all(project_dir: Path, receptor: Path, sdf: Path, *, padding: float=8.0, ga_runs: int=10, parallel_jobs: int=1, top_n: int=20, allow_failed_redocking_validation: bool = False, seed: int = 20250301, reference_ligand: str | None = None, allow_incomplete_peripheral: bool = False) -> dict[str, object]:
    receptor=receptor.resolve(); sdf=sdf.resolve()
    log('Creating project and detecting co-crystallized ligand')
    project=prepare_project(project_dir,receptor,sdf,padding,ga_runs,parallel_jobs,reference_ligand,allow_incomplete_peripheral)
    update_state(project.path, 'INPUT', 'PASS', inputs={'receptor': _sha256(receptor), 'ligands': _sha256(sdf)})
    project.docking.random_seed=seed
    _write_provenance(project, receptor, sdf, seed)
    selected = project.path / 'validation' / 'reference_ligand_resolution.json'
    ref = None
    if selected.exists():
        key = json.loads(selected.read_text(encoding='utf-8')).get('selected_instance', '')
        ref = next((g for g in hetero_groups(receptor) if g.key == key), None)
    ref = ref or choose_reference_ligand(receptor); assert ref and project.grid
    update_state(project.path, 'REFERENCE', 'PASS')
    log('Reading experimental IC50 metadata')
    assay=extract_sdf_assay(sdf,project.path/'analysis/assay.csv')
    (project.path/'qc/ligand_preparation.csv').write_text('compound_id,input_record_index,heavy_atoms,formal_charge,connected_components,3D_input,embedding_status,minimization_status,Meeko_status,PDBQT_charge_sum,active_torsions,warnings,failure_reason,final_status\n' + '\n'.join(f"{r.get('ligand','')},{r.get('input_record_index','')},,,,,,,,,,,{r.get('failure_reason','')},{r.get('final_status','')}" for r in assay) + '\n', encoding='utf-8')
    log('Preparing receptor')
    clean=clean_receptor_pdb(receptor,project.path/'receptor/receptor_clean.pdb',ref.resname,ref.chain,ref.resseq)
    receptor_log=project.path/'logs/prepare_receptor.log'
    try:
        rec=prepare_receptor(clean,project.path/'receptor/receptor.pdbqt',receptor_log)
    except RuntimeError:
        log('WARNING receptor preparation found unmatched residues; retrying with explicit Meeko --allow_bad_res. See receptor preparation audit.')
        rec=prepare_receptor(clean,project.path/'receptor/receptor.pdbqt',receptor_log,allow_bad_res=True)
        (project.path/'validation/receptor_preparation_fallback.json').write_text(json.dumps({'fallback':'Meeko --allow_bad_res','reason':'Initial template matching failed for incomplete/non-template residues; removed residues are listed in prepare_receptor.log.','scientific_review_required':True},indent=2),encoding='utf-8')
    project.receptor_pdbqt=str(rec)
    update_state(project.path, 'RECEPTOR', 'PASS', outputs={'receptor_pdbqt': _sha256(rec)})
    native_validation_dir = project.path/'analysis/validation'
    native_prepared = prepare_native_ligand(project.path/'receptor/reference_ligand.pdb', native_validation_dir)
    log('Correcting, embedding and minimizing ligands')
    normalized=split_and_minimize_sdf(sdf,project.path/'ligands/source')
    ensure_compounds(project.path, [p.stem for p in normalized])
    prepared=[]
    for src in normalized:
        dst=project.path/'ligands/prepared'/f'{src.stem}.pdbqt'
        try:
            prepare_ligand(src,dst,project.path/'logs'/f'prepare_{src.stem}.log'); prepared.append(dst)
            update_compound(project.path, src.stem, 'preparation', 'PASS', outputs={str(dst): _sha256(dst)})
        except Exception as exc:
            log(f'WARNING {src.stem}: {exc}')
            update_compound(project.path, src.stem, 'preparation', 'FAIL', error_code='E401', error_message=str(exc))
    if not prepared: raise RuntimeError('No ligands were successfully prepared.')
    project.prepared_ligands=[str(p) for p in prepared]; project.save()
    log('Running AutoGrid4')
    grid=project.path/'grid'; work=project.path/'docking/work'; dlgdir=project.path/'docking/dlg'; dpfdir=project.path/'docking/dpf'
    rec_grid=grid/rec.name; shutil.copy2(rec,rec_grid)
    ligand_types=sorted({t for p in prepared for t in pdbqt_atom_types(p)} | set(pdbqt_atom_types(native_prepared['pdbqt'])))
    gpf=grid/f'{rec_grid.stem}.gpf'; write_gpf(rec_grid,ligand_types,project.grid,gpf)
    autogrid=resolve_executable('autogrid4','AutoGrid4')
    if run_command([autogrid,'-p',gpf.name,'-l',f'{rec_grid.stem}.glg'],grid,project.path/'logs/autogrid.log') != 0:
        raise RuntimeError('AutoGrid4 failed; see logs/autogrid.log')
    for f in grid.glob(f'{rec_grid.stem}*.map*'): shutil.copy2(f,work/f.name)
    fld=grid/f'{rec_grid.stem}.maps.fld'
    if fld.exists(): shutil.copy2(fld,work/fld.name)
    log('Running native-ligand redocking validation')
    validation = run_native_redocking(native_pdb=project.path/'receptor/reference_ligand.pdb', receptor_pdbqt=rec, grid_dir=grid, work_dir=work/'native_validation', output_dir=project.path/'analysis/validation', grid=project.grid, params=project.docking, reference_info={'resname': ref.resname, 'chain': ref.chain, 'resseq': ref.resseq})
    update_state(project.path, 'VALIDATION', 'PASS' if validation.get('status') == 'PASS' else 'FAIL', warnings=[] if validation.get('status') == 'PASS' else ['Native redocking threshold exceeded'])
    provenance_path = project.path/'provenance.json'
    provenance = json.loads(provenance_path.read_text(encoding='utf-8')) if provenance_path.exists() else {}
    chemistry = validation.get('files', {}).get('reference_chemistry')
    if chemistry:
        provenance['reference_chemistry'] = chemistry
        provenance_path.write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    enforce_redocking_validation(validation, allow_failed_redocking_validation)
    autodock=resolve_executable('autodock4','AutoDock4')
    def dock(lig: Path):
        local=work/lig.name; shutil.copy2(lig,local)
        dpf=dpfdir/f'{lig.stem}.dpf'; write_dpf(rec_grid.stem,local,pdbqt_atom_types(local),project.docking,dpf)
        shutil.copy2(dpf,work/dpf.name)
        dlg=dlgdir/f'{lig.stem}.dlg'
        code=run_command([autodock,'-p',dpf.name,'-l',str(dlg.resolve())],work,project.path/'logs'/f'autodock_{lig.stem}.log')
        return lig.stem,code
    log('Running AutoDock4')
    with ThreadPoolExecutor(max_workers=max(1,parallel_jobs)) as pool:
        futures=[pool.submit(dock,p) for p in prepared]
        for f in as_completed(futures):
            name,code=f.result(); log(f'{name}: '+('done' if code==0 else 'failed'))
            if code != 0:
                raise RuntimeError(f'AutoDock4 failed for {name}; see logs/autodock_{name}.log')
    update_state(project.path, 'DOCKING', 'PASS', outputs={p.stem: _sha256(p) for p in dlgdir.glob('*.dlg')})
    compound_results=[]; ranking_rows=[]; cluster_rows=[]; plip_files=[]
    assay_by_id={str(r.get('ligand','')): r for r in assay}
    for dlg in sorted(dlgdir.glob('*.dlg')):
        cid=dlg.stem
        record=assay_by_id.get(cid, {'ligand':cid, 'final_status':'ACCEPTED'})
        try:
            poses=run_pose_extraction_stage(dlg, project.path/'analysis', top_n)
            update_compound(project.path,cid,'pose_extraction','PASS')
            clusters=run_cluster_stage(dlg, project.path/'analysis', top_n)
            update_compound(project.path,cid,'clustering','PASS')
            cluster_rows.extend(clusters.clusters)
            ranking_rows.append({'rank':0,'ligand':cid,'energy_kcal_mol':f'{poses.lowest_energy:.3f}' if poses.lowest_energy is not None else '', 'pose_file':str(poses.selected_pose_path or ''), 'complex_file':str(next((p['complex'] for p in poses.poses if p.get('model') == poses.selected_model),''))})
            plip_path=project.path/'analysis'/f'interactions_plip_{cid}.csv'
            complex_path=next((p['complex'] for p in poses.poses if p.get('model') == poses.selected_model), None)
            if complex_path is None: raise RuntimeError('E701 PLIP_FAILED')
            plip=run_plip_stage(complex_path, project.path/'analysis'/f'plip_{cid}', plip_path, compound_id=cid, pose_id=str(poses.selected_model))
            plip_files.append(plip.normalized_interactions_path)
            update_compound(project.path,cid,'plip','PASS')
            try:
                run_figures_stage(project.path/'analysis', project.path/'analysis/figures', limit=5)
                update_compound(project.path,cid,'figures','PASS')
            except Exception as exc:
                update_compound(project.path,cid,'figures','FAIL',error_code='E801',error_message=str(exc))
            record.update({'lowest_docking_energy':poses.lowest_energy,'cluster_count':clusters.cluster_count,'largest_cluster_population':clusters.largest_cluster_population,'selected_pose':str(poses.selected_pose_path or ''),'preparation_status':'PASS','PLIP_H_bonds':'','PLIP_hydrophobic':''})
        except Exception as exc:
            code='E701' if 'PLIP' in str(exc) else ('E602' if 'POSE' in str(exc) else ('E601' if 'DLG' in str(exc) else 'E801'))
            stage='plip' if code=='E701' else ('pose_extraction' if code=='E602' else ('clustering' if code=='E601' else 'figures'))
            update_compound(project.path,cid,stage,'FAIL',error_code=code,error_message=str(exc))
            record.update({'preparation_status':'PASS','status':'FAIL','failure_reason':str(exc)})
        compound_results.append(record)
    ranking_rows.sort(key=lambda r: float(r['energy_kcal_mol']) if r['energy_kcal_mol'] else 1e99)
    for rank,row in enumerate(ranking_rows,1): row['rank']=rank
    ranking_path=project.path/'analysis/ranking.csv'
    with ranking_path.open('w',newline='',encoding='utf-8') as fh:
        fields=['rank','ligand','energy_kcal_mol','pose_file','complex_file']; w=csv.DictWriter(fh,fieldnames=fields); w.writeheader(); w.writerows(ranking_rows)
    clusters_path=project.path/'analysis/clusters.csv'
    with clusters_path.open('w',newline='',encoding='utf-8') as fh:
        fields=['ligand','cluster_rank','cluster_size','cluster_percentage','lowest_binding_energy','mean_binding_energy','run_number','pose_model']; w=csv.DictWriter(fh,fieldnames=fields); w.writeheader(); w.writerows(cluster_rows)
    plip_combined=project.path/'analysis/interactions_plip.csv'
    if plip_files:
        with plip_combined.open('w',encoding='utf-8') as out:
            for i,path in enumerate(plip_files):
                text=path.read_text(encoding='utf-8')
                out.write(text if i==0 else '\n'.join(text.splitlines()[1:])+'\n')
    stats=correlate_with_ranking(assay,ranking_path,project.path/'analysis/assay_statistics.json')
    evidence={'schema_version':'0.4','interaction_engine':PRIMARY_INTERACTION_ENGINE,'sidock_internal_parser_status':SIDOCK_INTERACTION_STATUS,'validation':{'native_redocking_RMSD':validation.get('rmsd_A'),'sampling_result':validation.get('status')},'compounds':compound_results,'docking':{'ranking':str(ranking_path)},'clusters':str(clusters_path),'interactions':{'source':'PLIP','path':str(plip_combined)},'assay':assay,'statistics':stats,'warnings':[],'failures':[r for r in compound_results if r.get('status')=='FAIL'],'provenance':str(project.path/'provenance.json'),'grid':project.grid.__dict__,'assay_statistics':stats,'redocking_validation':validation,'files':{'ranking':str(ranking_path),'clusters':str(clusters_path)}}
    finalization=run_finalization_stage(project.path,evidence,write_results_fn=write_results,write_dashboard_fn=write_dashboard,compliance_fn=run_compliance)
    update_state(project.path, 'FINALIZATION', 'PASS', outputs={'evidence': _sha256(project.path/'analysis/evidence.json')})
    log('Pipeline complete')
    return {'project':project.path,'analysis':{'ranking':ranking_path,'clusters':clusters_path,'interactions':plip_combined},'statistics':stats,'finalization':finalization}

def run_stage_runner(project_dir: Path, receptor: Path, sdf: Path, *, config: dict, padding: float=8.0,
                     parallel_jobs: int=1, top_n: int=20, seed: int=20250301,
                     reference_ligand: str | None=None, allow_incomplete_peripheral: bool=False) -> dict[str, object]:
    """Bootstrap-only entrypoint; StageRunner owns all global execution."""
    project_dir=project_dir.resolve(); receptor=receptor.resolve(); sdf=sdf.resolve()
    project_json=project_dir/'project.json'
    if not project_json.exists():
        project=Project(str(project_dir),str(receptor),ligand_source=str(sdf),docking=DockingParameters(ga_runs=config['validation']['ga_runs'],parallel_jobs=1,random_seed=seed))
        project.ensure_layout(); project.save()
    ctx=production_context_from_project(project_dir,config) if project_json.exists() else ProductionContext(project_dir,config,receptor,sdf)
    ctx.rehydrate_global_outputs()
    registry=build_production_executors(ctx)
    runner=StageRunner(project_dir,config,registry)
    global_result=runner.run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    return {'project':project_dir,'global':global_result,'status':'BLOCKED' if any(s=='BLOCKED' for _,s in global_result) else 'RUNNING'}
