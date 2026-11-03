from __future__ import annotations
import hashlib, json, copy, threading
from datetime import datetime, timezone
from pathlib import Path

STAGES = ("INPUT", "REFERENCE", "RECEPTOR", "VALIDATION", "LIGAND_PREPARATION", "GRID", "DOCKING", "POSE_EXTRACTION", "CLUSTERING", "PLIP", "FIGURES", "FINALIZATION")
COMPOUND_STAGES = ("preparation", "docking", "pose_extraction", "clustering", "plip", "figures")
DEPENDENCIES = {"DOCKING": ["INPUT", "REFERENCE", "RECEPTOR", "VALIDATION", "LIGAND_PREPARATION", "GRID"], "POSE_EXTRACTION": ["DOCKING"], "CLUSTERING": ["POSE_EXTRACTION"], "PLIP": ["CLUSTERING"], "FIGURES": ["PLIP"], "FINALIZATION": ["FIGURES"]}
CONFIG_OWNERS = {"INPUT": ("input",), "REFERENCE": ("reference",), "RECEPTOR": ("receptor",), "VALIDATION": ("validation",), "LIGAND_PREPARATION": ("ligand_preparation",), "GRID": ("grid",), "DOCKING": ("docking",), "PLIP": ("interaction",), "FIGURES": ("figures",), "FINALIZATION": ("output",)}
_STATE_LOCK = threading.RLock()

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()

def config_hash(config: dict, stage: str) -> str:
    keys=CONFIG_OWNERS.get(stage, ())
    payload={k:config.get(k) for k in keys}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()

def compound_template() -> dict:
    return {s:{"status":"PENDING","started_at":None,"completed_at":None,"input_hashes":{},"output_hashes":{},"config_hash":None,"error_code":None,"error_message":None,"warnings":[],"attempts":[]} for s in COMPOUND_STAGES}

def ensure_compounds(project: Path, compound_ids: list[str]) -> dict:
    with _STATE_LOCK:
        path=project/'project_state.json'; state=project_status(project)
        state.setdefault('compounds', {})
        for cid in compound_ids: state['compounds'].setdefault(cid, compound_template())
        _atomic_write(path,state); return state

def update_compound(project: Path, compound: str, stage: str, status: str, *, config=None, inputs=None, outputs=None, error_code=None, error_message=None, warnings=None) -> dict:
    with _STATE_LOCK:
        state=project_status(project); state.setdefault('compounds',{}).setdefault(compound,compound_template())
        record=state['compounds'][compound][stage]; now=datetime.now(timezone.utc).isoformat()
        record.update({'status':status,'started_at':record.get('started_at') or now,'completed_at':now if status in {'PASS','FAIL','BLOCKED'} else None,'input_hashes':inputs or record.get('input_hashes',{}),'output_hashes':outputs or record.get('output_hashes',{}),'config_hash':config_hash(config, stage.upper()) if config else record.get('config_hash'),'error_code':error_code,'error_message':error_message,'warnings':warnings or record.get('warnings',[])})
        record.setdefault('attempts',[]).append({'attempt':len(record['attempts'])+1,'status':status,'error_code':error_code,'error_message':error_message,'timestamp':now})
        _atomic_write(project/'project_state.json',state); return state

def invalidate_from(state: dict, stage: str, reason: str) -> dict:
    names=list(STAGES); start=names.index(stage)
    for downstream in names[start:]:
        if downstream in state.get('stages',{}):
            state['stages'][downstream]['status']='STALE'; state['stages'][downstream].setdefault('warnings',[]).append(reason)
    return state

def verify_output_hashes(project: Path, stage_record: dict) -> bool:
    for name, expected in stage_record.get('output_hashes',{}).items():
        p=Path(name)
        if not p.is_absolute(): p=project/name
        if not p.exists() or sha256(p)!=expected: return False
    return True

def verify_input_hashes(project: Path, stage_record: dict) -> bool:
    """Verify persisted source fingerprints when paths are available."""
    for name, expected in stage_record.get('input_hashes',{}).items():
        if name in {'reference_config'}: continue
        p=Path(name)
        if not p.is_absolute(): p=project/name
        if expected is None or (not p.exists()): return False
        if sha256(p)!=expected: return False
    return True

def invalidate_for_config_change(state: dict, old_config: dict, new_config: dict) -> dict:
    for stage in STAGES:
        old=config_hash(old_config,stage); new=config_hash(new_config,stage)
        if old != new: invalidate_from(state,stage,f'CONFIG_CHANGED:{stage}')
    return state

def _atomic_write(path: Path, value: dict) -> None:
    with _STATE_LOCK:
        if path.exists(): path.replace(path.with_name('project_state.previous.json'))
        tmp=path.with_name(path.name+'.tmp'); tmp.write_text(json.dumps(value,indent=2,default=str),encoding='utf-8'); tmp.replace(path)

def update_state(project: Path, stage: str, status: str, *, inputs=None, outputs=None, errors=None, warnings=None, result_metadata=None) -> dict:
    with _STATE_LOCK:
        path = project / "project_state.json"
        state = json.loads(path.read_text()) if path.exists() else {"stages": {}, "updated": None}
        old_metadata=state.get("stages",{}).get(stage,{}).get("result_metadata")
        state.setdefault("stages", {})
        state["stages"][stage] = {"status": status, "timestamp": datetime.now(timezone.utc).isoformat(),
                                   "input_hashes": inputs or {}, "output_hashes": outputs or {},
                                   "errors": errors or [], "warnings": warnings or []}
        if result_metadata is not None:
            state["stages"][stage]["result_metadata"] = result_metadata
        elif old_metadata is not None:
            state["stages"][stage]["result_metadata"] = old_metadata
        state["updated"] = datetime.now(timezone.utc).isoformat()
        _atomic_write(path,state)
        return state

def project_status(project: Path) -> dict:
    with _STATE_LOCK:
        state_path = project / "project_state.json"
        if state_path.exists(): return json.loads(state_path.read_text())
        result = {"stages": {}, "updated": None}
        for stage in STAGES: result["stages"][stage] = {"status": "PENDING"}
        return result

def aggregate_compound_stage_status(state: dict, stage: str) -> str:
    """Derive a project-level compound-stage status without creating attempts."""
    compounds = state.get('compounds', {}) or {}
    validation = state.get('stages', {}).get('VALIDATION', {}).get('status')
    if validation in {'FAIL', 'BLOCKED', 'STALE'}:
        return 'BLOCKED'
    if not compounds:
        return 'PENDING'
    statuses = [c.get(stage, {}).get('status', 'PENDING') for c in compounds.values()]
    if all(x == 'PASS' for x in statuses): return 'PASS'
    if any(x == 'RUNNING' for x in statuses): return 'RUNNING'
    if any(x == 'PASS' for x in statuses) and any(x in {'FAIL', 'BLOCKED'} for x in statuses): return 'PARTIAL'
    if all(x == 'FAIL' for x in statuses): return 'FAIL'
    if any(x == 'BLOCKED' for x in statuses): return 'BLOCKED'
    return 'PENDING'

def derived_stage_statuses(state: dict) -> dict:
    result = dict(state.get('stages', {}))
    for stage in ('LIGAND_PREPARATION', 'DOCKING', 'POSE_EXTRACTION', 'CLUSTERING', 'PLIP', 'FIGURES'):
        result[stage] = dict(result.get(stage, {}))
        result[stage]['status'] = aggregate_compound_stage_status(state, 'preparation' if stage == 'LIGAND_PREPARATION' else stage.lower())
    return result

def derived_project_status(state: dict) -> str:
    stages = derived_stage_statuses(state)
    if not state.get('compounds'):
        stored = stages.get('FINALIZATION', {}).get('status')
        if stored in {'COMPLETE','COMPLETE_WITH_WARNINGS','PARTIAL','BLOCKED','STALE'}: return stored
    if stages.get('VALIDATION', {}).get('status') in {'FAIL', 'BLOCKED'}: return 'BLOCKED'
    final = stages.get('FINALIZATION', {}).get('status')
    required=('REFERENCE','RECEPTOR','GRID','VALIDATION','LIGAND_PREPARATION','DOCKING','POSE_EXTRACTION','CLUSTERING','PLIP','FIGURES')
    if final == 'PASS' and all(stages.get(x, {}).get('status') == 'PASS' for x in required):
        validation_warnings=stages.get('VALIDATION',{}).get('warnings',[]) or []
        return 'COMPLETE_WITH_WARNINGS' if 'UNVALIDATED_DOCKING' in validation_warnings else 'COMPLETE'
    if any(v.get('status') == 'FAIL' for v in stages.values() if isinstance(v, dict)): return 'FAILED'
    return 'INCOMPLETE'