from __future__ import annotations
import json
import time
from datetime import datetime, timezone
from dataclasses import dataclass
from pathlib import Path
from .state import COMPOUND_STAGES, config_hash, project_status, update_state, update_compound, verify_output_hashes, verify_input_hashes

@dataclass
class StageResult:
    status: str
    outputs: dict | None = None
    error_code: str | None = None
    error_message: str | None = None

class StageRunner:
    """Dependency-aware execution shell with injectable scientific stage callbacks."""
    def __init__(self, project: Path, config: dict, executors: dict | None = None):
        self.project=project; self.config=config; self.executors=executors or {}

    def _validation_mode(self) -> str:
        mode=str((self.config or {}).get('validation',{}).get('mode','if_available')).lower()
        return mode if mode in {'required','if_available','skip'} else 'if_available'

    def _read_validation_policy(self) -> str | None:
        try:
            mode=str(json.loads((self.project/'provenance'/'validation_policy.json').read_text(encoding='utf-8')).get('mode','')).lower()
            return mode if mode in {'required','if_available','skip'} else None
        except Exception:
            return None

    def _write_validation_policy(self, mode: str) -> None:
        path=self.project/'provenance'/'validation_policy.json'; path.parent.mkdir(parents=True,exist_ok=True)
        payload={'mode':mode,'rmsd_threshold':(self.config or {}).get('validation',{}).get('rmsd_threshold',2.0)}
        tmp=path.with_suffix('.json.tmp'); tmp.write_text(json.dumps(payload,indent=2,sort_keys=True),encoding='utf-8'); tmp.replace(path)

    def _persist_skipped_validation(self) -> Path:
        from .stage_outputs import ValidationStageOutput, save_validation_output
        from .state import sha256
        warnings=['UNVALIDATED_DOCKING','NATIVE_VALIDATION_SKIPPED_BY_POLICY']
        value=ValidationStageOutput(
            native_dlg=None,
            minimum_rmsd=None,
            median_rmsd=None,
            lowest_energy_rmsd=None,
            rmsd_threshold=(self.config or {}).get('validation',{}).get('rmsd_threshold',2.0),
            scientific_status='UNVALIDATED',
            execution_status='SKIPPED',
            native_like_pose_count=None,
            artifacts=[],
            warnings=warnings,
        )
        meta=save_validation_output(self.project,value)
        update_state(self.project,'VALIDATION','PASS',outputs={str(meta):sha256(meta)},warnings=warnings,result_metadata=str(meta.relative_to(self.project)))
        return meta

    def _trace(self, scope, stage, status, compound=None, attempt=None, error_code=None, started=None):
        path=self.project/'provenance/execution_trace.jsonl'; path.parent.mkdir(parents=True,exist_ok=True)
        row={'timestamp':datetime.now(timezone.utc).isoformat(),'scope':scope,'compound':compound,'stage':stage,'attempt':attempt,'executor':getattr(self.executors.get(stage.lower()),'__name__',stage),'status':status,'duration':(time.time()-started) if started else None,'error_code':error_code}
        with path.open('a',encoding='utf-8') as fh: fh.write(json.dumps(row)+'\n')

    def run_global(self, stages=('REFERENCE','RECEPTOR','GRID','VALIDATION','FINALIZATION')):
        """Execute reusable global callbacks while preserving dependency gates."""
        state=project_status(self.project); results=[]; passed=set(); downstream_dirty=False
        deps={'RECEPTOR':'REFERENCE','GRID':'RECEPTOR','VALIDATION':'GRID','FINALIZATION':'VALIDATION'}
        for stage in stages:
            dep=deps.get(stage)
            if dep and dep not in passed and state.get('stages',{}).get(dep,{}).get('status')!='PASS':
                update_state(self.project,stage,'BLOCKED',errors=[f'DEPENDENCY_{dep}']); results.append((stage,'BLOCKED')); break

            previous_validation_policy=None
            if stage=='VALIDATION':
                mode=self._validation_mode(); previous_validation_policy=self._read_validation_policy()
                if mode=='skip':
                    self._write_validation_policy(mode)
                    started=time.time(); self._persist_skipped_validation()
                    trace_path=self.project/'provenance/execution_trace.jsonl'; prior=0
                    if trace_path.exists():
                        prior=sum(1 for line in trace_path.read_text(encoding='utf-8').splitlines() if json.loads(line).get('scope')=='GLOBAL' and json.loads(line).get('stage')=='VALIDATION')
                    self._trace('GLOBAL','VALIDATION','SKIPPED',attempt=prior+1,started=started)
                    results.append(('VALIDATION','PASS')); passed.add('VALIDATION'); downstream_dirty=True
                    state=project_status(self.project)
                    continue

            current=state.get('stages',{}).get(stage,{}).get('status','PENDING')
            reusable=current=='PASS' and verify_output_hashes(self.project,state['stages'][stage])
            if reusable and stage=='REFERENCE': reusable=verify_input_hashes(self.project,state['stages'][stage])
            metadata_loaders={'REFERENCE':'load_reference_output','RECEPTOR':'load_receptor_output','GRID':'load_grid_output','VALIDATION':'load_validation_output','FINALIZATION':'load_finalization_output'}
            if reusable and stage in metadata_loaders and not state['stages'][stage].get('result_metadata'):
                reusable=False
            if reusable and stage in metadata_loaders:
                try:
                    module=__import__('sidock_ubuntu.finalization_outputs',fromlist=['load_finalization_output']) if stage=='FINALIZATION' else __import__('sidock_ubuntu.stage_outputs',fromlist=[metadata_loaders[stage]])
                    loaded=getattr(module,metadata_loaders[stage])(self.project)
                    reusable=getattr(loaded,'load_status','VALID')=='VALID'
                except (FileNotFoundError, ValueError, TypeError):
                    reusable=False
            if reusable and stage=='VALIDATION':
                reusable=previous_validation_policy==self._validation_mode()
            if reusable and not downstream_dirty:
                passed.add(stage); results.append((stage,'REUSED')); continue
            if stage=='VALIDATION': self._write_validation_policy(self._validation_mode())
            fn=self.executors.get(stage.lower()) or self.executors.get(stage)
            if fn is None: update_state(self.project,stage,'BLOCKED',errors=['MISSING_EXECUTOR']); results.append((stage,'BLOCKED')); break
            started=time.time(); update_state(self.project,stage,'RUNNING')
            try: raw=fn(); r=raw if isinstance(raw,StageResult) else StageResult('PASS',raw if isinstance(raw,dict) else {})
            except Exception as exc: r=StageResult('FAIL',error_code='E999',error_message=str(exc))
            persisted=project_status(self.project).get('stages',{}).get(stage,{})
            warnings=[]
            if stage=='VALIDATION' and self._validation_mode()=='if_available' and r.status=='PASS':
                out=r.outputs or {}; scientific=out.get('scientific_status',out.get('status'))
                if scientific not in {None,'PASS'}:
                    warnings=['UNVALIDATED_DOCKING','NATIVE_VALIDATION_FAILED_CONTINUING_BY_POLICY']
            update_state(self.project,stage,r.status,outputs=persisted.get('output_hashes',{}),errors=[x for x in (getattr(r,'error_code',None),getattr(r,'error_message',None)) if x],warnings=warnings,result_metadata=persisted.get('result_metadata'))
            trace_path=self.project/'provenance/execution_trace.jsonl'; prior=0
            if trace_path.exists():
                prior=sum(1 for line in trace_path.read_text(encoding='utf-8').splitlines() if json.loads(line).get('scope')=='GLOBAL' and json.loads(line).get('stage')==stage)
            self._trace('GLOBAL',stage,r.status,attempt=prior+1,error_code=getattr(r,'error_code',None),started=started); results.append((stage,r.status))
            if r.status!='PASS': break
            passed.add(stage)
            downstream_dirty=True
            state=project_status(self.project)
        return results

    def plan(self, compounds: list[str], *, failed_only=False, force=False) -> dict:
        state=project_status(self.project); result={}
        for compound in compounds:
            rec=state.get('compounds',{}).get(compound,{})
            result[compound]=[s for s in COMPOUND_STAGES if (not failed_only or rec.get(s,{}).get('status') in {'FAIL','BLOCKED'}) and (force or rec.get(s,{}).get('status','PENDING')!='PASS')]
        return result

    def _run_compound_unlocked(self, compound: str, *, failed_only=False, force=False) -> list[tuple[str,str]]:
        state=project_status(self.project); rec=state.get('compounds',{}).get(compound)
        if rec is None: return [('compound','MISSING')]
        global_validation=state.get('stages',{}).get('VALIDATION'); mode=self._validation_mode()
        # Validation execution errors still block. A scientific RMSD failure is
        # only a hard gate in required mode; other modes preserve it as an
        # explicit UNVALIDATED warning while allowing exploratory docking.
        if global_validation is not None and global_validation.get('status') in {'FAIL','BLOCKED','STALE'}:
            return [('preparation','BLOCKED')]
        if global_validation is not None and global_validation.get('status')=='PASS' and mode!='skip':
            try:
                from .stage_outputs import load_validation_output
                validation=load_validation_output(self.project)
                if getattr(validation,'load_status','STALE')!='VALID':
                    return [('preparation','BLOCKED')]
                if mode=='required' and validation.scientific_status!='PASS':
                    return [('preparation','BLOCKED')]
            except (FileNotFoundError, ValueError, TypeError):
                return [('preparation','BLOCKED')]
        result=[]
        for stage in COMPOUND_STAGES:
            status=rec.get(stage,{}).get('status','PENDING')
            if failed_only and status not in {'FAIL','BLOCKED'}: continue
            if status=='PASS' and not force and verify_output_hashes(self.project,rec.get(stage,{})): continue
            if stage!='preparation':
                prev=COMPOUND_STAGES[COMPOUND_STAGES.index(stage)-1]
                if rec.get(prev,{}).get('status')!='PASS': result.append((stage,'BLOCKED')); break
            fn=self.executors.get(stage)
            if fn is None: result.append((stage,'BLOCKED')); break
            update_compound(self.project,compound,stage,'RUNNING',config=self.config)
            try:
                r=fn(compound); r=r if isinstance(r,StageResult) else StageResult('PASS',r or {})
            except Exception as exc: r=StageResult('FAIL',error_code='E999',error_message=str(exc))
            update_compound(self.project,compound,stage,r.status,config=self.config,outputs=r.outputs,error_code=r.error_code,error_message=r.error_message)
            result.append((stage,r.status)); rec=project_status(self.project)['compounds'][compound]
            self._trace('COMPOUND',stage,r.status,compound=compound,error_code=r.error_code,started=None)
            if r.status!='PASS': break
        return result

    def run_compound(self, compound: str, *, failed_only=False, force=False) -> list[tuple[str,str]]:
        from .project_lock import project_execution_lock, ProjectExecutionLocked
        with project_execution_lock(self.project, {'interface':'CLI','action':'COMPOUND','compound_id':compound}):
            return self._run_compound_unlocked(compound,failed_only=failed_only,force=force)

def config_change_audit(project: Path, old: dict, new: dict) -> dict:
    owners=('INPUT','REFERENCE','RECEPTOR','VALIDATION','LIGAND_PREPARATION','GRID','DOCKING','PLIP','FIGURES','FINALIZATION')
    names=('INPUT','REFERENCE','RECEPTOR','VALIDATION','LIGAND_PREPARATION','GRID','DOCKING','POSE_EXTRACTION','CLUSTERING','PLIP','FIGURES','FINALIZATION')
    changes={}; invalid=[]
    for stage in owners:
        a,b=config_hash(old,stage),config_hash(new,stage)
        if a!=b: changes[stage]={'old_hash':a,'new_hash':b}; invalid.append(stage)
    expanded=[]
    for stage in invalid:
        for name in names[names.index(stage):]:
            if name not in expanded: expanded.append(name)
    result={'changed':bool(changes),'changes':changes,'invalidated_stages':expanded}
    (project/'provenance').mkdir(parents=True,exist_ok=True)
    (project/'provenance/config_change_audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result
