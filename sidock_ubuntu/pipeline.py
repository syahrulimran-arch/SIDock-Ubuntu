from __future__ import annotations

import csv, json, shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .analysis import analyze_docking_results
from .assay import correlate_with_ranking, extract_sdf_assay
from .core import (box_from_reference, clean_receptor_pdb, pdbqt_atom_types, prepare_ligand,
                   prepare_receptor, resolve_executable, run_command, split_and_minimize_sdf,
                   write_dpf, write_gpf)
from .figures import activity_scatter, pymol_binding_figures
from .intake import choose_reference_ligand, protein_intake_summary, write_group
from .model import DockingParameters, Project


def log(msg: str) -> None:
    print(f"[SIDock Ubuntu] {msg}", flush=True)


def prepare_project(project_dir: Path, receptor: Path, sdf: Path, padding: float, ga_runs: int, jobs: int) -> Project:
    project=Project(str(project_dir.resolve()), str(receptor.resolve()), ligand_source=str(sdf.resolve()), docking=DockingParameters(ga_runs=ga_runs, parallel_jobs=jobs))
    project.ensure_layout()
    summary=protein_intake_summary(receptor)
    (project.path/'intake.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    ref=choose_reference_ligand(receptor)
    if not ref:
        raise RuntimeError('No suitable co-crystallized ligand detected. Blind pocket selection is disabled in v0.1.')
    ref_path=write_group(ref, project.path/'receptor/reference_ligand.pdb')
    project.grid=box_from_reference(ref_path,padding=padding)
    project.save(); return project


def run_all(project_dir: Path, receptor: Path, sdf: Path, *, padding: float=8.0, ga_runs: int=10, parallel_jobs: int=1, top_n: int=20) -> dict[str, object]:
    receptor=receptor.resolve(); sdf=sdf.resolve()
    log('Creating project and detecting co-crystallized ligand')
    project=prepare_project(project_dir,receptor,sdf,padding,ga_runs,parallel_jobs)
    ref=choose_reference_ligand(receptor); assert ref and project.grid
    log('Reading experimental IC50 metadata')
    assay=extract_sdf_assay(sdf,project.path/'analysis/assay.csv')
    log('Preparing receptor')
    clean=clean_receptor_pdb(receptor,project.path/'receptor/receptor_clean.pdb',ref.resname,ref.chain,ref.resseq)
    rec=prepare_receptor(clean,project.path/'receptor/receptor.pdbqt',project.path/'logs/prepare_receptor.log')
    project.receptor_pdbqt=str(rec)
    log('Correcting, embedding and minimizing ligands')
    normalized=split_and_minimize_sdf(sdf,project.path/'ligands/source')
    prepared=[]
    for src in normalized:
        dst=project.path/'ligands/prepared'/f'{src.stem}.pdbqt'
        try:
            prepare_ligand(src,dst,project.path/'logs'/f'prepare_{src.stem}.log'); prepared.append(dst)
        except Exception as exc:
            log(f'WARNING {src.stem}: {exc}')
    if not prepared: raise RuntimeError('No ligands were successfully prepared.')
    project.prepared_ligands=[str(p) for p in prepared]; project.save()
    log('Running AutoGrid4')
    grid=project.path/'grid'; work=project.path/'docking/work'; dlgdir=project.path/'docking/dlg'; dpfdir=project.path/'docking/dpf'
    rec_grid=grid/rec.name; shutil.copy2(rec,rec_grid)
    ligand_types=sorted({t for p in prepared for t in pdbqt_atom_types(p)})
    gpf=grid/f'{rec_grid.stem}.gpf'; write_gpf(rec_grid,ligand_types,project.grid,gpf)
    autogrid=resolve_executable('autogrid4','AutoGrid4')
    if run_command([autogrid,'-p',gpf.name,'-l',f'{rec_grid.stem}.glg'],grid,project.path/'logs/autogrid.log') != 0:
        raise RuntimeError('AutoGrid4 failed; see logs/autogrid.log')
    for f in grid.glob(f'{rec_grid.stem}*.map*'): shutil.copy2(f,work/f.name)
    fld=grid/f'{rec_grid.stem}.maps.fld'
    if fld.exists(): shutil.copy2(fld,work/fld.name)
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
    analysis=analyze_docking_results(dlgdir,project.path/'analysis',top_n)
    stats=correlate_with_ranking(assay,Path(analysis['ranking']),project.path/'analysis/assay_statistics.json')
    activity_scatter(stats,project.path/'analysis/figures/docking_vs_activity.png')
    pymol_binding_figures(project.path/'analysis',project.path/'analysis/figures',limit=5)
    evidence={'rules':['Only discuss generated values.','Do not invent residues/interactions/scores.','Docking is predictive, not experimental proof.'],
              'grid':project.grid.__dict__,'assay_statistics':stats,'files':{k:str(v) for k,v in analysis.items() if isinstance(v,Path)}}
    (project.path/'analysis/evidence.json').write_text(json.dumps(evidence,indent=2),encoding='utf-8')
    (project.path/'analysis/CODEX_DISCUSSION_PROMPT.md').write_text('# SIDock Ubuntu writing task\nRead evidence.json and generated CSV files. Write Results_and_Discussion.md. Every numerical or interaction claim must be traceable to generated outputs. Never invent data. Distinguish docking predictions from experimental evidence and state limitations.\n',encoding='utf-8')
    log('Pipeline complete')
    return {'project':project.path,'analysis':analysis,'statistics':stats}
