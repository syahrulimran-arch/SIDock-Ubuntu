"""Repair/rebuild post-processing from immutable docking artifacts.

This module never runs AutoGrid or AutoDock. It reuses existing DLG files to
regenerate pose extraction, cluster tables, PLIP interactions and final report
artifacts for completed projects.
"""
from __future__ import annotations
import argparse
import json
import shutil
from pathlib import Path

from .config import load_config
from .executors import ProductionContext, build_production_executors
from .finalization import load_finalization_input, snapshot_evidence
from .orchestrator import StageRunner
from .stage_outputs import load_receptor_output
from .stages import run_pose_extraction_stage, run_cluster_stage, run_plip_stage
from .state import project_status, sha256, update_compound, update_state


def _hash_output(path: Path | None) -> dict[str,str]:
    return {str(path):sha256(path)} if path is not None and Path(path).is_file() else {}


def rebuild_results(project: Path) -> dict:
    project=Path(project).resolve()
    if not (project/'project_state.json').exists():
        raise FileNotFoundError(f'Not a SIDock project: {project}')
    state=project_status(project)
    compounds=sorted(state.get('compounds',{}))
    if not compounds:
        raise RuntimeError('No compound execution records found in project state.')
    analysis=project/'analysis'; analysis.mkdir(parents=True,exist_ok=True)
    summary={'project':str(project),'docking_rerun':False,'compounds':{},'finalization':None}

    for cid in compounds:
        dlg=project/'docking'/'dlg'/f'{cid}.dlg'
        row={'dlg':str(dlg),'pose_extraction':'NOT_RUN','clustering':'NOT_RUN','plip':'NOT_RUN'}
        if not dlg.is_file():
            row.update({'pose_extraction':'FAILED','clustering':'FAILED','plip':'FAILED','error':'DLG_MISSING'})
            summary['compounds'][cid]=row
            continue
        try:
            pose=run_pose_extraction_stage(dlg,analysis)
            if pose.selected_pose_path is None:
                raise RuntimeError('E602 POSE_MAPPING_FAILED')
            update_compound(project,cid,'pose_extraction','PASS',outputs=_hash_output(Path(pose.selected_pose_path)))
            row['pose_extraction']='PASS'
        except Exception as exc:
            update_compound(project,cid,'pose_extraction','FAIL',error_code='E602',error_message=str(exc))
            row.update({'pose_extraction':'FAILED','error':str(exc)})
            summary['compounds'][cid]=row
            continue

        try:
            cluster=run_cluster_stage(dlg,analysis)
            cluster_path=Path(cluster.cluster_table_path) if cluster.cluster_table_path else None
            update_compound(project,cid,'clustering','PASS',outputs=_hash_output(cluster_path))
            row.update({'clustering':'PASS','cluster_count':cluster.cluster_count,'largest_cluster_population':cluster.largest_cluster_population})
        except Exception as exc:
            update_compound(project,cid,'clustering','FAIL',error_code='E601',error_message=str(exc))
            row.update({'clustering':'FAILED','error':str(exc)})
            summary['compounds'][cid]=row
            continue

        try:
            selected=next((p for p in pose.poses if p.get('model')==pose.selected_model),None)
            complex_path=Path(selected['complex']) if selected and selected.get('complex') else None
            if complex_path is None or not complex_path.is_file():
                raise RuntimeError('E701 PLIP_FAILED: selected pose complex missing')
            work=analysis/f'plip_{cid}'
            if work.exists(): shutil.rmtree(work)
            normalized=analysis/f'interactions_plip_{cid}.csv'
            if normalized.exists(): normalized.unlink()
            plip=run_plip_stage(complex_path,work,normalized,compound_id=cid,pose_id=str(pose.selected_model))
            update_compound(project,cid,'plip','PASS',outputs=_hash_output(Path(plip.normalized_interactions_path)))
            row['plip']='PASS'
        except Exception as exc:
            update_compound(project,cid,'plip','FAIL',error_code='E701',error_message=str(exc))
            row.update({'plip':'FAILED','error':str(exc)})
        summary['compounds'][cid]=row

    update_state(project,'FINALIZATION','STALE',warnings=['POSTPROCESSING_REBUILT_FROM_EXISTING_DLG'])
    config=load_config()
    ctx=ProductionContext(project,config)
    runner=StageRunner(project,config,build_production_executors(ctx))
    final=runner.run_global(('FINALIZATION',))
    summary['finalization']=final
    evidence=snapshot_evidence(load_finalization_input(project))
    summary['project_status']=evidence.get('project_status')
    summary['target']=evidence.get('target')
    summary['validation']=evidence.get('validation')
    return summary


def main() -> None:
    parser=argparse.ArgumentParser(description='Rebuild SIDock post-processing from existing DLG files without rerunning docking.')
    parser.add_argument('project',type=Path)
    args=parser.parse_args()
    result=rebuild_results(args.project)
    print(json.dumps(result,indent=2,default=str))


if __name__=='__main__':
    main()
