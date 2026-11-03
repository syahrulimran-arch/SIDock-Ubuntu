"""Thin, validated-backend stage seams.

These functions delegate to existing scientific implementations. They do not
change command construction or scientific algorithms. PLIP is intentionally
not represented here until a real source callable is recovered.
"""
from __future__ import annotations
from pathlib import Path
from .core import box_from_reference, clean_receptor_pdb, prepare_receptor, split_and_minimize_sdf, prepare_ligand, write_gpf, write_dpf, run_command
from .analysis import extract_docking_poses, analyze_autodock_clusters
from .validation import prepare_native_ligand, run_native_redocking
from .figures import pymol_binding_figures
from .plip_runner import run_plip_analysis

def run_reference_stage(pdb: Path, project: Path, reference_group):
    from .intake import write_group
    out=project/'receptor/reference_ligand.pdb'; write_group(reference_group,out)
    return {'reference_pdb':out,'grid':box_from_reference(out)}

def run_receptor_stage(pdb: Path, project: Path, reference_group, *, allow_bad_res=False):
    clean=clean_receptor_pdb(pdb,project/'receptor/receptor_clean.pdb',reference_group.resname,reference_group.chain,reference_group.resseq)
    rec=prepare_receptor(clean,project/'receptor/receptor.pdbqt',project/'logs/prepare_receptor.log',allow_bad_res=allow_bad_res)
    return {'clean':clean,'receptor_pdbqt':rec}

def run_compound_preparation(sdf: Path, project: Path):
    sources=split_and_minimize_sdf(sdf,project/'ligands/source'); outputs=[]
    for source in sources:
        target=project/'ligands/prepared'/f'{source.stem}.pdbqt'
        outputs.append(prepare_ligand(source,target,project/'logs'/f'prepare_{source.stem}.log'))
    return outputs

def run_validation_stage(native_pdb, receptor_pdbqt, grid_dir, work_dir, output_dir, grid, params, reference_info, ccd_path=None):
    native=prepare_native_ligand(native_pdb,output_dir,ccd_path=ccd_path)
    return run_native_redocking(native_pdb=native_pdb,receptor_pdbqt=receptor_pdbqt,grid_dir=grid_dir,work_dir=work_dir,output_dir=output_dir,grid=grid,params=params,reference_info=reference_info,ccd_path=ccd_path)

def run_grid_stage(receptor_pdbqt: Path, grid_dir: Path, grid, ligand_types: list[str], autogrid: str):
    """Generate one role-specific AutoGrid grid and persist exact provenance."""
    import hashlib
    import json
    import shutil

    grid_dir.mkdir(parents=True,exist_ok=True)
    local_receptor=grid_dir/'receptor.pdbqt'
    if receptor_pdbqt.resolve()!=local_receptor.resolve():
        shutil.copy2(receptor_pdbqt,local_receptor)

    gpf=grid_dir/'receptor.gpf'
    write_gpf(local_receptor,ligand_types,grid,gpf)
    glg=grid_dir/'receptor.glg'
    code=run_command([autogrid,'-p',gpf.name,'-l',glg.name],grid_dir,grid_dir/'autogrid.log')
    if code!=0:
        raise RuntimeError('E501 AUTOGRID_FAILED')

    maps=sorted(p for p in grid_dir.glob('receptor.*.map*') if p.is_file())
    role=grid_dir.name if grid_dir.name in {'production','validation'} else 'legacy'
    receptor_hash=hashlib.sha256(local_receptor.read_bytes()).hexdigest()
    provenance={
        'role': role,
        'receptor_pdbqt': str(local_receptor.resolve()),
        'receptor_sha256': receptor_hash,
        'ligand_types': sorted(set(ligand_types)),
        'center': [grid.center_x,grid.center_y,grid.center_z],
        'dimensions': [grid.npts_x,grid.npts_y,grid.npts_z],
        'spacing': grid.spacing,
        'gpf': str(gpf.resolve()),
        'glg': str(glg.resolve()),
        'maps': [str(p.resolve()) for p in maps],
    }
    provenance_path=grid_dir/'grid_provenance.json'
    tmp=provenance_path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(provenance,indent=2,sort_keys=True),encoding='utf-8')
    tmp.replace(provenance_path)

    return {
        'gpf':gpf,
        'glg':glg,
        'maps':maps,
        'center':provenance['center'],
        'dimensions':provenance['dimensions'],
        'spacing':provenance['spacing'],
        'grid_role':role,
        'grid_provenance':provenance_path,
    }

def run_compound_docking(ligand_pdbqt: Path, receptor_stem: str, work_dir: Path, dpf_dir: Path, dlg_dir: Path, ligand_types: list[str], params, autodock: str):
    """Run one AutoDock job in a compound-private execution directory.

    AutoDock4 itself is a single CPU process. Isolating each compound's working
    directory allows many AutoDock processes to run concurrently without racing
    on ligand, DPF, map, or log files. Production GA runs are read from the
    current run provenance and are independent of the number of CPU workers.
    """
    import copy
    import json
    import shutil

    dpf_dir.mkdir(parents=True,exist_ok=True); dlg_dir.mkdir(parents=True,exist_ok=True); work_dir.mkdir(parents=True,exist_ok=True)
    project_root=dpf_dir.parent.parent
    run_dir=work_dir/ligand_pdbqt.stem; run_dir.mkdir(parents=True,exist_ok=True)
    local_ligand=run_dir/ligand_pdbqt.name
    shutil.copy2(ligand_pdbqt,local_ligand)

    grid_source=project_root/'grid'/'production'
    if not grid_source.is_dir(): grid_source=project_root/'grid'
    for artifact in grid_source.glob(f'{receptor_stem}.*.map*'):
        if artifact.is_file(): shutil.copy2(artifact,run_dir/artifact.name)
    fld=grid_source/f'{receptor_stem}.maps.fld'
    if fld.is_file(): shutil.copy2(fld,run_dir/fld.name)

    effective_params=copy.copy(params)
    settings_path=project_root/'provenance'/'parallel_execution.json'
    if settings_path.is_file():
        try:
            settings=json.loads(settings_path.read_text(encoding='utf-8'))
            ga_runs=int(settings.get('ga_runs_per_compound',effective_params.ga_runs))
            if ga_runs >= 1: effective_params.ga_runs=ga_runs
        except (OSError,ValueError,TypeError,json.JSONDecodeError):
            pass

    dpf=dpf_dir/f'{ligand_pdbqt.stem}.dpf'; write_dpf(receptor_stem,local_ligand,ligand_types,effective_params,dpf)
    shutil.copy2(dpf,run_dir/dpf.name)
    dlg=dlg_dir/f'{ligand_pdbqt.stem}.dlg'
    code=run_command([autodock,'-p',dpf.name,'-l',str(dlg.resolve())],run_dir,run_dir/f'{ligand_pdbqt.stem}.log')
    if code!=0 or not dlg.exists(): raise RuntimeError('E502 AUTODOCK_FAILED')
    return {'dpf':dpf,'dlg':dlg,'work_dir':run_dir,'ga_runs':effective_params.ga_runs}

def _single_dlg(dlg_or_dir: Path) -> Path | None:
    if dlg_or_dir.is_file():
        return dlg_or_dir
    return next(iter(sorted(dlg_or_dir.glob('*.dlg'))),None)

def run_pose_extraction_stage(dlg_dir: Path, analysis_dir: Path, top_n=20):
    """Extract poses using the existing analysis implementation."""
    dlg=_single_dlg(dlg_dir)
    if dlg is None: raise RuntimeError('E602 POSE_MAPPING_FAILED')
    receptor=dlg_dir.parent.parent/'receptor'/'receptor.pdbqt'
    return extract_docking_poses(dlg,analysis_dir/'poses',analysis_dir/'complexes',receptor)

def run_cluster_stage(dlg_dir: Path, analysis_dir: Path, top_n=20):
    """Read one DLG-derived cluster table and persist it per compound."""
    dlg=_single_dlg(dlg_dir)
    if dlg is None: raise RuntimeError('E601 DLG_PARSE_FAILED')
    result=analyze_autodock_clusters(dlg,analysis_dir/f'clusters_{dlg.stem}.csv')
    return result

def run_pose_and_cluster_stage(dlg_dir: Path, analysis_dir: Path, top_n=20):
    pose=run_pose_extraction_stage(dlg_dir,analysis_dir,top_n)
    cluster=run_cluster_stage(dlg_dir,analysis_dir,top_n)
    return {**pose, **cluster}

def run_plip_stage(complex_pdb: Path, work_dir: Path, normalized_csv: Path, *, compound_id='', pose_id=''):
    """Use the recovered Open Babel + PLIP production path; never contacts.py."""
    return run_plip_analysis(complex_pdb, work_dir, normalized_csv, compound_id=compound_id, pose_id=pose_id)

def run_finalization_stage(project: Path, evidence: dict | None = None, *, write_results_fn, write_dashboard_fn, compliance_fn=None):
    """Finalize machine-readable project outputs; no manuscript generation."""
    if evidence is None:
        from .finalization import load_finalization_input, snapshot_evidence
        evidence=snapshot_evidence(load_finalization_input(project))
    analysis=project/'analysis'; analysis.mkdir(parents=True,exist_ok=True)
    evidence_tmp=analysis/'evidence.json.tmp'; evidence_tmp.write_text(__import__('json').dumps(evidence,indent=2,sort_keys=True),encoding='utf-8'); evidence_tmp.replace(analysis/'evidence.json')
    write_results_fn(project)
    # Pure post-processing: reads evidence/DLG/PLIP; never executes docking.
    from .scientific_report import write_scientific_report
    write_scientific_report(project)
    write_dashboard_fn(project)
    return compliance_fn(project) if compliance_fn else {'status':'NOT_RUN'}

def run_figures_stage(analysis_dir: Path, figures_dir: Path, limit=5):
    return pymol_binding_figures(analysis_dir,figures_dir,limit=limit)
