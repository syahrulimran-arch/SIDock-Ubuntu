"""Disk-only reconstruction of project results for finalization."""
from __future__ import annotations
import json
import csv
from dataclasses import dataclass, field
from pathlib import Path
from .state import project_status
from .analysis import parse_dlg_models, parse_dlg_clusters
from .stage_outputs import load_validation_output, load_reference_output, load_receptor_output, load_grid_output

@dataclass
class CompoundResultSnapshot:
    compound_id: str
    stages: dict
    warnings: list[str] = field(default_factory=list)
    errors: list[dict] = field(default_factory=list)

@dataclass
class ProjectResultSnapshot:
    project: Path
    global_outputs: dict
    compounds: list[CompoundResultSnapshot]
    warnings: list[str] = field(default_factory=list)
    status: str = 'PENDING'

def _global(project: Path) -> dict:
    out={}
    for key, loader in (('reference',load_reference_output),('receptor',load_receptor_output),('grid',load_grid_output),('validation',load_validation_output)):
        try: out[key]=loader(project)
        except (FileNotFoundError, ValueError, TypeError): out[key]=None
    return out

def _validation_policy(project: Path) -> str:
    path=project/'provenance'/'validation_policy.json'
    try:
        mode=str(json.loads(path.read_text(encoding='utf-8')).get('mode','required')).lower()
    except Exception:
        mode='required'
    return mode if mode in {'required','if_available','skip'} else 'required'

def _target_summary(project: Path, outputs: dict) -> dict:
    """Recover receptor and reference identity from persisted project metadata."""
    pdb_id=None; receptor_source=None
    try:
        project_data=json.loads((project/'project.json').read_text(encoding='utf-8'))
        receptor_source=project_data.get('receptor_source')
    except Exception:
        project_data={}
    try:
        gui=json.loads((project/'gui_inputs.json').read_text(encoding='utf-8'))
        pdb_id=gui.get('pdb_id') or None
    except Exception:
        pass
    ref=outputs.get('reference')
    selected=getattr(ref,'selected_ligand',None) if ref is not None else None
    if not selected:
        for candidate in (project/'validation'/'reference_selection.json', project/'validation'/'reference_ligand_resolution.json'):
            try:
                data=json.loads(candidate.read_text(encoding='utf-8'))
                selected=data.get('selected_instance') or data.get('selected_ligand') or data.get('reference_ligand')
                if selected: break
            except Exception:
                continue
    receptor_label=pdb_id or (Path(receptor_source).name if receptor_source else None)
    return {
        'pdb': receptor_label,
        'pdb_id': pdb_id,
        'receptor_source': receptor_source,
        'reference_ligand': selected,
        'reference_residue_name': getattr(ref,'residue_name',None) if ref is not None else None,
        'reference_chain': getattr(ref,'chain',None) if ref is not None else None,
        'reference_residue_number': getattr(ref,'residue_number',None) if ref is not None else None,
    }

def load_finalization_input(project: Path) -> ProjectResultSnapshot:
    state=project_status(project); compounds=[]
    for cid, stages in state.get('compounds',{}).items():
        errors=[]; warnings=[]
        for name, rec in stages.items():
            if not isinstance(rec,dict): continue
            if rec.get('error_code') or rec.get('error_message'):
                errors.append({'stage':name,'error_code':rec.get('error_code'),'message':rec.get('error_message')})
            warnings.extend(rec.get('warnings',[]) or [])
        compounds.append(CompoundResultSnapshot(cid,stages,warnings,errors))
    outputs=_global(project); validation=outputs.get('validation'); mode=_validation_policy(project)
    validation_scientific=getattr(validation,'scientific_status',None) if validation is not None else None
    validation_warnings=state.get('stages',{}).get('VALIDATION',{}).get('warnings',[]) or []
    unvalidated=(mode=='skip' or 'UNVALIDATED_DOCKING' in validation_warnings or (mode=='if_available' and validation_scientific not in {None,'PASS'}))
    warnings=list(validation_warnings)
    if validation_scientific=='FAIL' and mode=='required': status='BLOCKED'
    elif any(any(s.get('status')=='FAIL' for s in c.stages.values() if isinstance(s,dict)) for c in compounds): status='PARTIAL'
    elif compounds and all(all(s.get('status')=='PASS' for s in c.stages.values() if isinstance(s,dict)) for c in compounds): status='COMPLETE_WITH_WARNINGS' if unvalidated else 'COMPLETE'
    else: status='PENDING'
    return ProjectResultSnapshot(project,outputs,compounds,warnings,status)

def _cluster_rows(project: Path, cid: str) -> list[dict]:
    """Load per-compound cluster evidence, falling back to immutable DLG parsing."""
    candidates=(project/'analysis'/f'clusters_{cid}.csv', project/'analysis'/'clusters.csv')
    for path in candidates:
        if not path.exists(): continue
        try:
            with path.open(newline='',encoding='utf-8') as fh:
                rows=[r for r in csv.DictReader(fh) if r.get('ligand')==cid]
            if rows: return rows
        except (OSError,csv.Error):
            pass
    dlg=project/'docking'/'dlg'/f'{cid}.dlg'
    if dlg.exists():
        try: return parse_dlg_clusters(dlg)
        except (OSError,ValueError,TypeError): return []
    return []

def _disk_compound_summary(project: Path, cid: str, stages: dict) -> dict:
    """Reconstruct the canonical summary solely from persisted stage artifacts."""
    out = {'lowest_docking_energy': None, 'selected_pose': None,
           'cluster_count': None, 'largest_cluster_population': None,
           'PLIP_H_bonds': None, 'PLIP_hydrophobic': None,
           'PLIP_pi_stack': None, 'PLIP_salt_bridge': None}
    ranking = project/'analysis'/'ranking.csv'
    if ranking.exists():
        with ranking.open(newline='', encoding='utf-8') as fh:
            row = next((r for r in csv.DictReader(fh) if r.get('ligand') == cid), None)
        if row:
            out['lowest_docking_energy'] = float(row['energy_kcal_mol']) if row.get('energy_kcal_mol') else None
            out['selected_pose'] = row.get('pose_file') or None
    if out['lowest_docking_energy'] is None:
        dlg = project/'docking'/'dlg'/f'{cid}.dlg'
        if dlg.exists():
            models = [m for m in parse_dlg_models(dlg) if m.get('energy') is not None]
            if models:
                best=min(models, key=lambda m:(float(m['energy']), int(m['model'])))
                out['lowest_docking_energy']=float(best['energy'])
                pose=project/'analysis'/'poses'/f"{cid}_model_{best['model']}.pdbqt"
                out['selected_pose']=str(pose) if pose.exists() else None
    rows=_cluster_rows(project,cid)
    if rows:
        out['cluster_count']=len(rows)
        sizes=[]
        for r in rows:
            try: sizes.append(int(r.get('cluster_size',0)))
            except (TypeError,ValueError): pass
        out['largest_cluster_population']=max(sizes) if sizes else None
    plip = project/'analysis'/f'interactions_plip_{cid}.csv'
    if plip.exists() and stages.get('plip',{}).get('status') == 'PASS':
        counts={}
        with plip.open(newline='', encoding='utf-8') as fh:
            for row in csv.DictReader(fh):
                kind=(row.get('interaction_type') or '').strip()
                if kind: counts[kind]=counts.get(kind,0)+1
        out.update({'PLIP_H_bonds': counts.get('hydrogen_bond',0),
                    'PLIP_hydrophobic': counts.get('hydrophobic',0),
                    'PLIP_pi_stack': counts.get('pi-stack',0),
                    'PLIP_salt_bridge': counts.get('salt_bridge',0)})
    return out

def snapshot_evidence(snapshot: ProjectResultSnapshot) -> dict:
    compounds=[]
    for c in snapshot.compounds:
        stages=c.stages
        disk=_disk_compound_summary(snapshot.project, c.compound_id, stages)
        if snapshot.status=='BLOCKED': compound_status='BLOCKED'
        else: compound_status=('COMPLETE' if all(s.get('status')=='PASS' for s in stages.values() if isinstance(s,dict)) else 'FAILED')
        compounds.append({'compound_id':c.compound_id,'status':compound_status,**({k:None for k in disk} if compound_status=='BLOCKED' else disk),'warnings':c.warnings,'errors':c.errors})
    validation=snapshot.global_outputs.get('validation'); mode=_validation_policy(snapshot.project)
    scientific=getattr(validation,'scientific_status',None) if validation else None
    execution=getattr(validation,'execution_status',None) if validation else None
    docking_validation_status='VALIDATED' if scientific=='PASS' else 'UNVALIDATED'
    val={'sampling_result':scientific or ('UNVALIDATED' if mode=='skip' else None),
         'native_redocking_RMSD_A':getattr(validation,'minimum_rmsd',None) if validation else None,
         'rmsd_threshold':getattr(validation,'rmsd_threshold',None) if validation else None,
         'execution_status':execution or ('SKIPPED' if mode=='skip' else None),
         'validation_policy':mode,
         'docking_validation_status':docking_validation_status}
    blocked_reason='NATIVE_REDOCKING_RMSD_THRESHOLD' if snapshot.status=='BLOCKED' else None
    return {'schema_version':'0.6','project_status':snapshot.status,
            'target':_target_summary(snapshot.project,snapshot.global_outputs),
            'interaction_engine':'PLIP','sidock_internal_parser_status':'EXPERIMENTAL_SECONDARY',
            'validation':val,'blocked_by':'VALIDATION' if snapshot.status=='BLOCKED' else None,
            'blocked_reason':blocked_reason,'compounds':compounds,
            'failures':[x for x in compounds if x['status']=='FAILED'],'warnings':snapshot.warnings}
