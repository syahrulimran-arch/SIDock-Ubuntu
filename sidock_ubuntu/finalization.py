"""Disk-only reconstruction of project results for finalization."""
from __future__ import annotations
import json
from dataclasses import dataclass, field
from pathlib import Path
from .state import project_status
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
    validation=out_validation=_global(project).get('validation')
    if validation is not None and getattr(validation,'scientific_status',None)=='FAIL': status='BLOCKED'
    elif any(any(s.get('status')=='FAIL' for s in c.stages.values() if isinstance(s,dict)) for c in compounds): status='PARTIAL'
    elif compounds and all(all(s.get('status')=='PASS' for s in c.stages.values() if isinstance(s,dict)) for c in compounds): status='COMPLETE'
    else: status='PENDING'
    return ProjectResultSnapshot(project,_global(project),compounds,[],status)

def snapshot_evidence(snapshot: ProjectResultSnapshot) -> dict:
    compounds=[]
    for c in snapshot.compounds:
        stages=c.stages; plip=stages.get('plip',{}); docking=stages.get('docking',{}); pose=stages.get('pose_extraction',{}); cluster=stages.get('clustering',{})
        plip_ok=plip.get('status')=='PASS'
        if snapshot.status=='BLOCKED': compound_status='BLOCKED'
        else: compound_status=('COMPLETE' if all(s.get('status')=='PASS' for s in stages.values() if isinstance(s,dict)) else 'FAILED')
        compounds.append({'compound_id':c.compound_id,'status':compound_status,'lowest_docking_energy':None if compound_status=='BLOCKED' else docking.get('lowest_energy'),'selected_pose':None if compound_status=='BLOCKED' else pose.get('selected_pose'),'cluster_count':None if compound_status=='BLOCKED' else cluster.get('cluster_count'),'largest_cluster_population':None if compound_status=='BLOCKED' else cluster.get('largest_cluster_population'),'PLIP_H_bonds':None if (not plip_ok or compound_status=='BLOCKED') else plip.get('hbond_count'),'PLIP_hydrophobic':None if (not plip_ok or compound_status=='BLOCKED') else plip.get('hydrophobic_count'),'warnings':c.warnings,'errors':c.errors})
    validation=snapshot.global_outputs.get('validation')
    val={'sampling_result':getattr(validation,'scientific_status',None),'native_redocking_RMSD_A':getattr(validation,'minimum_rmsd',None),'rmsd_threshold':getattr(validation,'rmsd_threshold',None)} if validation else {}
    blocked_reason='NATIVE_REDOCKING_RMSD_THRESHOLD' if snapshot.status=='BLOCKED' else None
    return {'schema_version':'0.4','project_status':snapshot.status,'interaction_engine':'PLIP','sidock_internal_parser_status':'EXPERIMENTAL_SECONDARY','validation':val,'blocked_by':'VALIDATION' if snapshot.status=='BLOCKED' else None,'blocked_reason':blocked_reason,'compounds':compounds,'failures':[x for x in compounds if x['status']=='FAILED'],'warnings':snapshot.warnings}
