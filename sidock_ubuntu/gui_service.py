"""Thin service boundary for the v0.5 local GUI.

All scientific execution is delegated to the existing pipeline/StageRunner;
this module only owns safe storage, state reads, and process coordination.
"""
from __future__ import annotations

import copy
import json
import shutil
import threading
import csv
import os, multiprocessing
from pathlib import Path

from .config import load_config, VALIDATION_MODES
from .model import DockingParameters, Project
from .pipeline import run_stage_runner
from .orchestrator import available_cpu_count
from .state import STAGES, project_status, derived_project_status, invalidate_from, _atomic_write
from .job_store import ACTIVE, create_job, load_job, update_job, project_lock, process_alive, process_identity, list_jobs
from .project_lock import ProjectExecutionLocked, project_execution_lock

MAX_UPLOAD_BYTES = int(os.environ.get('SIDOCK_GUI_MAX_UPLOAD_BYTES', str(100 * 1024 * 1024)))

def validate_uploaded_inputs(protein: Path, compounds: Path, mmcif: Path | None = None, ccd: Path | None = None) -> None:
    """Use the installed scientific parsers for upload admission only."""
    if protein.stat().st_size == 0 or not any(x.startswith(('ATOM', 'HETATM', 'HEADER', 'MODEL')) for x in protein.read_text(errors='replace').splitlines()):
        raise ValueError('E001 INVALID_PDB: no structural records found')
    try:
        from rdkit import Chem
        supplier = Chem.SDMolSupplier(str(compounds), sanitize=True, removeHs=False)
        if not any(mol is not None for mol in supplier): raise ValueError('E002 INVALID_SDF: no valid molecules found')
    except ValueError: raise
    except Exception as exc: raise ValueError(f'E002 INVALID_SDF: {exc}')
    if mmcif:
        try:
            import gemmi; gemmi.cif.read_file(str(mmcif))
        except Exception as exc: raise ValueError(f'E001 INVALID_MMCIF: {exc}')
    if ccd:
        try:
            if ccd.suffix.lower() == '.json': json.loads(ccd.read_text())
            else:
                import gemmi; gemmi.cif.read_file(str(ccd))
        except Exception as exc: raise ValueError(f'E102 INVALID_CCD: {exc}')

def get_doctor_status() -> dict:
    from .cli import doctor_status
    return doctor_status()

def get_config(path: Path | None = None) -> dict:
    return load_config(path)

def validate_config_values(values: dict) -> dict:
    from .config import deep_merge, DEFAULT_CONFIG, validate_config
    try:
        resolved=deep_merge(DEFAULT_CONFIG, values or {}); errors=validate_config(resolved)
    except (TypeError, ValueError, KeyError) as exc:
        return {'valid':False,'errors':[str(exc)]}
    return {'valid':not errors,'errors':errors,'config':resolved}

def get_project_config(project: Path) -> dict:
    path=project/'provenance'/'resolved_config.yaml'
    if not path.exists(): return {'config':None,'path':None,'source':'not persisted'}
    try:
        import yaml
        return {'config':yaml.safe_load(path.read_text(encoding='utf-8')) or {},'path':str(path),'source':'project resolved configuration'}
    except Exception as exc:
        return {'config':None,'path':str(path),'source':'invalid persisted configuration','error':str(exc)}

def input_compound_count(project: Path) -> tuple[int, list[str]]:
    """Count valid molecules using the same RDKit admission path as uploads."""
    try:
        from rdkit import Chem
        data=json.loads((Path(project)/'project.json').read_text(encoding='utf-8'))
        source=Path(data.get('ligand_source',''))
        if not source.exists(): return 0, []
        molecules=[m for m in Chem.SDMolSupplier(str(source), sanitize=True, removeHs=False) if m is not None]
        ids=[m.GetProp('_Name') if m.HasProp('_Name') and m.GetProp('_Name') else str(i+1) for i,m in enumerate(molecules)]
        return len(molecules), ids
    except Exception:
        return 0, []

_jobs: dict[str, multiprocessing.Process] = {}
_jobs_lock = threading.Lock()


def safe_project_name(name: str) -> str:
    if "/" in str(name) or "\\" in str(name) or ".." in str(name):
        raise ValueError("E001 INVALID_PROJECT_NAME")
    cleaned = "".join(ch for ch in str(name).strip() if ch.isalnum() or ch in "-_ ").strip()
    cleaned = cleaned.replace(" ", "-")
    if not cleaned or cleaned in {".", ".."} or cleaned.startswith("."):
        raise ValueError("E001 INVALID_PROJECT_NAME")
    return cleaned[:80]


def workspace_path(workspace: Path, name: str) -> Path:
    root = workspace.resolve()
    target = (root / safe_project_name(name)).resolve()
    if root not in target.parents:
        raise ValueError("E001 PROJECT_PATH_OUTSIDE_WORKSPACE")
    return target


def _copy_input(source: Path, destination: Path, suffixes: tuple[str, ...]) -> Path:
    source = Path(source).resolve()
    if not source.is_file() or source.suffix.lower() not in suffixes:
        raise ValueError("E001 INVALID_INPUT_FILE")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return destination


def create_project(workspace: Path, name: str, protein: Path, compounds: Path,
                   mmcif: Path | None = None, ccd: Path | None = None,
                   validation_structure: Path | None = None, pdb_id: str | None = None) -> Path:
    project = workspace_path(workspace, name)
    if project.exists():
        raise FileExistsError("E003 PROJECT_ALREADY_EXISTS")
    project.mkdir(parents=True)
    try:
        protein_dst = _copy_input(protein, project / "input" / "protein.pdb", (".pdb", ".ent"))
        sdf_dst = _copy_input(compounds, project / "input" / "compounds.sdf", (".sdf",))
        mmcif_dst = _copy_input(mmcif, project / "input" / "structure.cif", (".cif", ".mmcif")) if mmcif else None
        ccd_dst = _copy_input(ccd, project / "input" / "ccd.json", (".json", ".cif")) if ccd else None
        validation_dst = _copy_input(validation_structure, project / "input" / "validation_structure.pdb", (".pdb", ".ent", ".cif", ".mmcif")) if validation_structure else None
        defaults=load_config(); dock_defaults=defaults.get('docking',{})
        model = Project(str(project), str(protein_dst), ligand_source=str(sdf_dst), docking=DockingParameters(ga_runs=int(dock_defaults.get('ga_runs',10))))
        model.ensure_layout(); model.save()
        from .intake import protein_intake_summary
        (project / "intake.json").write_text(json.dumps(protein_intake_summary(protein_dst, mmcif=mmcif_dst), indent=2), encoding="utf-8")
        from .structure_context import structure_context
        context=structure_context(protein_dst, validation_dst)
        (project / "provenance" / "structure_context.json").write_text(json.dumps(context, indent=2), encoding="utf-8")
        state = {"stages": {stage: {"status": "PENDING"} for stage in STAGES}, "compounds": {}, "updated": None}
        (project / "project_state.json").write_text(json.dumps(state, indent=2), encoding="utf-8")
        default_mode=defaults.get('validation',{}).get('mode','required')
        (project / "gui_inputs.json").write_text(json.dumps({"mmcif": str(mmcif_dst) if mmcif_dst else None, "ccd": str(ccd_dst) if ccd_dst else None, "validation_structure": str(validation_dst) if validation_dst else None, "pdb_id": pdb_id or None, "validation_mode": default_mode, "ga_runs": int(dock_defaults.get('ga_runs',10)), "parallel_jobs": dock_defaults.get('parallel_jobs','auto')}, indent=2), encoding="utf-8")
        return project
    except Exception:
        shutil.rmtree(project, ignore_errors=True)
        raise

def attach_project_inputs(project: Path, *, validation_structure: Path | None = None,
                          mmcif: Path | None = None, ccd: Path | None = None,
                          confirm_completed: bool = False) -> dict:
    """Atomically attach authoritative structure inputs and invalidate dependents."""
    project=Path(project).resolve()
    state=project_status(project)
    if state.get('stages',{}).get('FINALIZATION',{}).get('status')=='PASS' and not confirm_completed:
        raise ValueError('E103 COMPLETED_PROJECT_INPUT_CHANGE_REQUIRES_CONFIRMATION')
    data=json.loads((project/'project.json').read_text())
    protein=Path(data['receptor_source']); compounds=Path(data['ligand_source'])
    if not any((validation_structure,mmcif,ccd)): raise ValueError('E001 NO_INPUT_CHANGE')
    temp=project/'input'/'.pending-inputs'; temp.mkdir(parents=True,exist_ok=True)
    try:
        def stage(source, name, suffixes):
            if source is None: return None
            source=Path(source)
            if source.suffix.lower() not in suffixes or not source.is_file(): raise ValueError('E001 INVALID_INPUT_FILE')
            target=temp/name; shutil.copy2(source,target); return target
        v=stage(validation_structure,'validation_structure',('.pdb','.ent','.cif','.mmcif'))
        m=stage(mmcif,'structure.cif',('.cif','.mmcif'))
        c=stage(ccd,'ccd.json',('.json','.cif'))
        validate_uploaded_inputs(protein,compounds,m,c)
        for staged, final in ((v,project/'input'/'validation_structure.pdb'),(m,project/'input'/'structure.cif'),(c,project/'input'/'ccd.json')):
            if staged:
                final.parent.mkdir(parents=True,exist_ok=True); os.replace(staged,final)
        inputs=json.loads((project/'gui_inputs.json').read_text()) if (project/'gui_inputs.json').exists() else {}
        if v: inputs['validation_structure']=str(project/'input'/'validation_structure.pdb')
        if m: inputs['mmcif']=str(project/'input'/'structure.cif')
        if c: inputs['ccd']=str(project/'input'/'ccd.json')
        _atomic_write(project/'gui_inputs.json',inputs)
        from .intake import protein_intake_summary
        _atomic_write(project/'intake.json',protein_intake_summary(protein,mmcif=Path(inputs['mmcif']) if inputs.get('mmcif') else None))
        from .structure_context import structure_context
        validation=Path(inputs['validation_structure']) if inputs.get('validation_structure') else None
        _atomic_write(project/'provenance'/'structure_context.json',structure_context(protein,validation))
        state=invalidate_from(state,'REFERENCE','AUTHORITATIVE_STRUCTURE_INPUT_CHANGED')
        for record in state.get('compounds',{}).values():
            for stage in record.values():
                if isinstance(stage,dict): stage['status']='STALE'; stage.setdefault('warnings',[]).append('AUTHORITATIVE_STRUCTURE_INPUT_CHANGED')
        _atomic_write(project/'project_state.json',state)
        return {'project':str(project),'updated':['validation_structure' if v else None,'mmcif' if m else None,'ccd' if c else None]}
    finally:
        shutil.rmtree(temp,ignore_errors=True)


def reference_view(project: Path) -> dict:
    validation = project / "validation"
    result = {}
    for filename, key in (("reference_ligand_candidates.csv", "candidates_csv"), ("reference_selection.json", "selection")):
        path = validation / filename
        if path.exists():
            result[key] = path.read_text(encoding="utf-8") if key.endswith("csv") else json.loads(path.read_text(encoding="utf-8"))
    try:
        result['structure_intake']=json.loads((project/'intake.json').read_text())
    except Exception:
        result['structure_intake']={}
    try:
        result['structure_context']=json.loads((project/'provenance'/'structure_context.json').read_text())
    except Exception:
        result['structure_context']={}
    return result

def set_reference_override(project: Path, instance: str) -> dict:
    view=reference_view(project); candidates=view.get('selection',{}).get('candidate_ranking',[])
    if candidates and not any(x.get('key')==instance or x.get('instance_id')==instance for x in candidates):
        raise ValueError('E101 INVALID_REFERENCE_OVERRIDE')
    path=project/'.gui_reference_override'; path.write_text(str(instance),encoding='utf-8')
    return {'selection_method':'EXPLICIT_OVERRIDE','reference_ligand':instance}


def _gui_inputs(project: Path) -> dict:
    try: return json.loads((project/'gui_inputs.json').read_text(encoding='utf-8'))
    except Exception: return {}


def _stored_validation_mode(project: Path) -> str:
    mode=str(_gui_inputs(project).get('validation_mode','')).lower()
    if mode not in VALIDATION_MODES:
        mode=str(load_config().get('validation',{}).get('mode','required')).lower()
    return mode if mode in VALIDATION_MODES else 'required'


def _stored_ga_runs(project: Path) -> int:
    value=_gui_inputs(project).get('ga_runs',load_config().get('docking',{}).get('ga_runs',10))
    try: return max(1,int(value))
    except (TypeError,ValueError): return 10


def _stored_parallel_jobs(project: Path):
    value=_gui_inputs(project).get('parallel_jobs',load_config().get('docking',{}).get('parallel_jobs','auto'))
    if str(value).lower()=='auto': return 'auto'
    try: return max(1,int(value))
    except (TypeError,ValueError): return 'auto'


def dashboard_view(project: Path, *, doctor: dict | None = None) -> dict:
    state = project_status(project)
    stored_final=state.get('stages',{}).get('FINALIZATION',{}).get('status','PENDING')
    final=derived_project_status(state)
    job=load_job(project); active=bool(job and job.get('status') in ACTIVE)
    ref=reference_view(project); selection=ref.get('selection',{}) if isinstance(ref.get('selection'),dict) else {}
    confidence=selection.get('confidence'); reference_status='RESOLVED' if selection.get('selected_instance') else ('LOW' if confidence=='LOW' else ('NOT_RUN' if not ref else 'UNRESOLVED'))
    validation=state.get('stages',{}).get('VALIDATION',{}); compounds=get_compounds(project); statuses=[x['status'] for x in compounds]
    input_count, _ = input_compound_count(project)
    locked=False
    try:
        with project_execution_lock(project, {'interface':'GUI','action':'PROBE'}): pass
    except ProjectExecutionLocked: locked=True
    dependency_ready=True if doctor is None else doctor.get('status') in {'READY','READY_WITH_WARNINGS'}
    intake={}; context={}
    try: intake=json.loads((project/'intake.json').read_text())
    except Exception: pass
    try: context=json.loads((project/'provenance'/'structure_context.json').read_text())
    except Exception: pass
    structure_class=intake.get('structure_classification'); policy=_stored_validation_mode(project)
    validation_attached=bool(context.get('validation_structure',{}).get('path')) and context.get('validation_structure',{}).get('path') != context.get('production_receptor',{}).get('path')
    if active or locked: action, url, message='RUNNING', None, 'Project currently executing.'
    elif reference_status=='LOW': action, url, message='SELECT_REFERENCE', f'/projects/{project.name}/reference', 'Reference ligand selection requires an explicit choice.'
    elif structure_class in {'PROCESSED_STRIPPED','APO_NO_REFERENCE','ION_ONLY'} and not validation_attached: action, url, message='UPLOAD_VALIDATION_STRUCTURE', f'/projects/{project.name}/reference-inputs', 'No reference ligand/site is available to define the docking box. Attach a matching holo/reference structure. Native redocking itself may still be optional.'
    elif context.get('validation_mode')=='SURROGATE_HOLO_VALIDATION': action, url, message='NO_NATIVE_VALIDATION', f'/projects/{project.name}/reference-inputs', 'The attached structure is related but not equivalent and cannot be silently used as the native reference/grid source.'
    elif final=='BLOCKED': action, url, message='VIEW_VALIDATION', f'/projects/{project.name}/results', 'Strict native redocking validation blocked this run.'
    elif job and job.get('status') in {'INTERRUPTED','ORPHANED'} or stored_final in {'STALE'}: action, url, message='RESUME', f'/projects/{project.name}/resume', 'Persisted scientific state can be resumed.'
    elif final in {'COMPLETE','COMPLETE_WITH_WARNINGS','PARTIAL'}: action, url, message='VIEW_RESULTS', f'/projects/{project.name}/results', 'View persisted results.'
    elif job and job.get('status')=='FAILED': action, url, message='RESUME', f'/projects/{project.name}/resume', 'The previous execution failed and may be resumed.'
    elif not dependency_ready: action, url, message='NONE', '/doctor', 'Execution environment is not ready.'
    else: action, url, message='RUN', f'/projects/{project.name}/run', 'Ready to start SIDock execution.'
    return {'project':str(project.resolve()),'project_status':final,'operational_job_status':job.get('status') if job else None,'reference_status':reference_status,'reference_confidence':confidence,'validation_status':validation.get('status','NOT_RUN'),'validation_rmsd':validation.get('rmsd_A') or validation.get('minimum_rmsd'),'validation_threshold':validation.get('threshold') or validation.get('rmsd_threshold'),'structure_classification':structure_class,'validation_mode':context.get('validation_mode'),'validation_policy':policy,'ga_runs':_stored_ga_runs(project),'parallel_jobs':_stored_parallel_jobs(project),'detected_cpus':available_cpu_count(),'compound_total':len(compounds),'input_compound_count':input_count,'compound_success':statuses.count('PASS'),'compound_failed':statuses.count('FAIL'),'compound_blocked':statuses.count('BLOCKED'),'primary_action':action,'primary_action_url':url,'explanatory_message':message,'warnings':_evidence(project).get('warnings',[]),'job':job,'active_job':active or locked,'reference':ref}

def project_view(project: Path) -> dict:
    return {"project": str(project.resolve()), "state": project_status(project), "reference": reference_view(project), 'dashboard':dashboard_view(project)}

def _evidence(project: Path) -> dict:
    path=project/'analysis'/'evidence.json'
    return json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}

def _compound_status(record: dict) -> str:
    statuses=[v.get('status','PENDING') for v in record.values() if isinstance(v,dict) and 'status' in v]
    if 'FAIL' in statuses: return 'FAIL'
    if 'BLOCKED' in statuses: return 'BLOCKED'
    if 'STALE' in statuses: return 'STALE'
    if statuses and all(x=='PASS' for x in statuses): return 'PASS'
    if 'RUNNING' in statuses: return 'RUNNING'
    return 'PENDING'

def get_compounds(project: Path) -> list[dict]:
    state=project_status(project); evidence={x.get('compound_id'):x for x in _evidence(project).get('compounds',[])}
    rows=[]
    for cid, record in sorted(state.get('compounds',{}).items()):
        ev=evidence.get(cid,{})
        rows.append({'compound_id':cid,'status':_compound_status(record),'stages':record,'docking_energy':ev.get('lowest_docking_energy'),'error_code':next((v.get('error_code') for v in record.values() if isinstance(v,dict) and v.get('error_code')),None),'warnings':ev.get('warnings',[])})
    return rows

def get_compound_page(project: Path, *, page=1, page_size=25, status='', search='', sort='compound_id') -> dict:
    page=max(1,int(page)); page_size=int(page_size) if int(page_size) in {25,50,100} else 25
    rows=get_compounds(project)
    if status and status != 'ALL': rows=[x for x in rows if x['status']==status]
    if search: rows=[x for x in rows if search.lower() in str(x['compound_id']).lower()]
    if sort=='energy': rows.sort(key=lambda x: (x.get('docking_energy') is None, x.get('docking_energy') if x.get('docking_energy') is not None else 0, str(x['compound_id'])))
    elif sort=='status': rows.sort(key=lambda x: (str(x['status']),str(x['compound_id'])))
    else: rows.sort(key=lambda x: str(x['compound_id']))
    total=len(rows); start=(page-1)*page_size
    return {'items':rows[start:start+page_size],'page':page,'page_size':page_size,'total':total,'pages':max(1,(total+page_size-1)//page_size)}

def get_compound_detail(project: Path, compound_id: str) -> dict:
    for row in get_compounds(project):
        if row['compound_id']==compound_id:
            evidence=next((x for x in _evidence(project).get('compounds',[]) if x.get('compound_id')==compound_id),{})
            return {**row,'evidence':evidence}
    raise KeyError(compound_id)

def get_project_results(project: Path) -> dict:
    state=project_status(project); evidence=_evidence(project)
    compounds=get_compounds(project); statuses=[x['status'] for x in compounds]
    project_status_value=derived_project_status(state)
    if evidence.get('project_status') in {'COMPLETE_WITH_WARNINGS','BLOCKED','PARTIAL'}:
        project_status_value=evidence['project_status']
    elif project_status_value == 'INCOMPLETE' and evidence.get('project_status'):
        project_status_value=evidence['project_status']
    return {'project':str(project.resolve()),'project_status':project_status_value,'validation':evidence.get('validation',{}),'reference':reference_view(project),'compound_counts':{x:statuses.count(x) for x in ('PASS','FAIL','BLOCKED','STALE','PENDING','RUNNING')},'compounds':compounds,'warnings':evidence.get('warnings',[]),'finalized':state.get('stages',{}).get('FINALIZATION',{}).get('status')=='PASS'}

def resolve_project_artifact(project: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute() or '..' in Path(relative).parts:
        raise ValueError('E001 INVALID_ARTIFACT_PATH')
    root=project.resolve(); target=(root/relative).resolve()
    if root != target and root not in target.parents: raise ValueError('E001 ARTIFACT_OUTSIDE_PROJECT')
    if not target.is_file(): raise FileNotFoundError('ARTIFACT_NOT_FOUND')
    return target


def _worker_entry(project_str, job_id, config, reference_ligand):
    project=Path(project_str); log=project/'provenance/gui_worker.log'
    try:
        with log.open('a',encoding='utf-8') as fh: fh.write(f'job {job_id} STARTING pid={os.getpid()}\n')
        update_job(project,status='STARTING',pid=os.getpid(),process_start_identity=process_identity(os.getpid()),started_at=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),backend_state_at_start=project_status(project))
        update_job(project,status='RUNNING')
        data=json.loads((project/'project.json').read_text(encoding='utf-8')); override=reference_ligand or ((project/'.gui_reference_override').read_text(encoding='utf-8').strip() if (project/'.gui_reference_override').exists() else None)
        run_stage_runner(project,Path(data['receptor_source']),Path(data['ligand_source']),config=config or load_config(),reference_ligand=override)
        update_job(project,status='SUCCEEDED',exit_code=0,finished_at=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),backend_state_at_finish=project_status(project))
        with log.open('a',encoding='utf-8') as fh: fh.write(f'job {job_id} SUCCEEDED\n')
    except ProjectExecutionLocked as exc:
        update_job(project,status='FAILED',exit_code=75,last_error={'type':type(exc).__name__,'message':'PROJECT_LOCKED'},finished_at=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),backend_state_at_finish=project_status(project))
        with log.open('a',encoding='utf-8') as fh: fh.write('PROJECT_LOCKED\n')
    except Exception as exc:
        update_job(project,status='FAILED',exit_code=1,last_error={'type':type(exc).__name__,'message':str(exc)},finished_at=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),backend_state_at_finish=project_status(project))
        with log.open('a',encoding='utf-8') as fh: fh.write(f'{type(exc).__name__}: {exc}\n')

def reconcile_jobs(workspace: Path) -> list[dict]:
    results=[]
    root=Path(workspace).resolve()
    for project in sorted(root.iterdir()) if root.exists() else []:
        if not project.is_dir(): continue
        job=load_job(project)
        if not job or job.get('status') not in ACTIVE: continue
        if process_alive(job.get('pid'),job.get('process_start_identity')): results.append(job); continue
        state=project_status(project); terminal={'COMPLETE','COMPLETE_WITH_WARNINGS','PARTIAL','BLOCKED'}
        statuses={x.get('status') for x in state.get('stages',{}).values()}
        update_job(project,status='SUCCEEDED' if statuses & terminal else 'INTERRUPTED',finished_at=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),backend_state_at_finish=state,last_error=None if statuses & terminal else {'type':'WorkerExited','message':'Worker process was not found during reconciliation'})
        results.append(load_job(project))
    return results

def get_job(project):
    job=load_job(project)
    return {'job':job,'project_status':get_project_results(project).get('project_status','PENDING'),'can_resume':bool(job and job.get('status') in {'INTERRUPTED','FAILED','ORPHANED'})} if job else {'job':None,'project_status':get_project_results(project).get('project_status','PENDING'),'can_resume':True}

def get_job_history(project): return list_jobs(project)

def get_worker_log(project, limit=12000):
    path=Path(project)/'provenance/gui_worker.log'
    if not path.exists(): return ''
    return path.read_text(encoding='utf-8',errors='replace')[-int(limit):]

def start_run(project: Path, *, config: dict | None = None, reference_ligand: str | None = None, validation_mode: str | None = None, ga_runs: int | None = None, parallel_jobs=None, action='RUN') -> dict:
    project=project.resolve()
    if get_doctor_status().get('status') == 'NOT_READY':
        raise RuntimeError('Execution environment is not ready; see Doctor for missing dependencies.')
    data=json.loads((project/'project.json').read_text(encoding='utf-8'))
    if not Path(data.get('receptor_source','')).is_file(): raise RuntimeError('Required protein input is missing.')
    count, _ = input_compound_count(project)
    if count < 1: raise RuntimeError('No valid compounds were found in the project SDF input.')
    run_config=copy.deepcopy(config or load_config())
    mode=str(validation_mode or _stored_validation_mode(project) or run_config.get('validation',{}).get('mode','required')).lower()
    if mode not in VALIDATION_MODES: raise RuntimeError('E001 INVALID_VALIDATION_MODE')
    run_config.setdefault('validation',{})['mode']=mode

    try: effective_ga=max(1,int(ga_runs if ga_runs is not None else _stored_ga_runs(project)))
    except (TypeError,ValueError): raise RuntimeError('E001 INVALID_GA_RUNS')
    requested_parallel=parallel_jobs if parallel_jobs is not None else _stored_parallel_jobs(project)
    if str(requested_parallel).lower()=='auto': effective_parallel='auto'
    else:
        try: effective_parallel=max(1,int(requested_parallel))
        except (TypeError,ValueError): raise RuntimeError('E001 INVALID_PARALLEL_JOBS')
    run_config.setdefault('docking',{})['ga_runs']=effective_ga
    run_config['docking']['parallel_jobs']=effective_parallel

    inputs=_gui_inputs(project); inputs.update({'validation_mode':mode,'ga_runs':effective_ga,'parallel_jobs':effective_parallel}); _atomic_write(project/'gui_inputs.json',inputs)
    try: intake=json.loads((project/'intake.json').read_text())
    except Exception: intake={}
    try: context=json.loads((project/'provenance'/'structure_context.json').read_text())
    except Exception: context={}
    if intake.get('structure_classification') in {'PROCESSED_STRIPPED','APO_NO_REFERENCE','ION_ONLY'} and not context.get('validation_structure',{}).get('path'):
        raise RuntimeError('E101 NO_REFERENCE_OR_GRID_DEFINITION: no usable ligand/site is available to define the docking box. Attach a matching holo/reference structure. Native redocking validation may still be skipped.')
    if context.get('validation_mode')=='SURROGATE_HOLO_VALIDATION':
        raise RuntimeError('E102 SURROGATE_STRUCTURE_NOT_NATIVE: a related but different structure cannot be silently used as the native reference/grid source')
    with _jobs_lock:
        try:
            with project_lock(project,'reservation'):
                job=create_job(project,action=action)
                proc=multiprocessing.Process(target=_worker_entry,args=(str(project),job['job_id'],run_config,reference_ligand),daemon=False); proc.start()
                update_job(project,pid=proc.pid,process_start_identity=process_identity(proc.pid))
                _jobs[str(project)]=proc
        except BlockingIOError as exc: raise RuntimeError('E999 PROJECT_ALREADY_RUNNING') from exc
    return {**job,'status':'STARTING','pid':proc.pid,'validation_mode':mode,'ga_runs':effective_ga,'parallel_jobs':effective_parallel,'detected_cpus':available_cpu_count()}


def resume_project(project: Path, *, failed_only=False, compound=None) -> dict:
    """Delegate resume to the same CLI/backend machinery in a worker."""
    return start_run(project, config=load_config(), validation_mode=_stored_validation_mode(project), ga_runs=_stored_ga_runs(project), parallel_jobs=_stored_parallel_jobs(project), action='RESUME')
