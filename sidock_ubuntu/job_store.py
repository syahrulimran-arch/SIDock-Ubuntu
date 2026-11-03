"""Atomic operational job records and local project locks."""
from __future__ import annotations
import json, os, socket, time, uuid
from datetime import datetime, timezone
from pathlib import Path
from .project_lock import project_execution_lock, process_identity, ProjectExecutionLocked

SCHEMA_VERSION='0.5-gui-job-1'
ACTIVE={'QUEUED','STARTING','RUNNING','CANCEL_REQUESTED'}

def now(): return datetime.now(timezone.utc).isoformat()
def _write(path:Path,data:dict):
    path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_suffix(path.suffix+'.tmp'); tmp.write_text(json.dumps(data,indent=2,sort_keys=True,default=str),encoding='utf-8'); tmp.replace(path)
def job_path(project): return Path(project)/'provenance/gui_job.json'
def history_path(project): return Path(project)/'provenance/gui_job_history.jsonl'
def load_job(project):
    path=job_path(project)
    if not path.exists(): return None
    try: return json.loads(path.read_text(encoding='utf-8'))
    except (OSError,json.JSONDecodeError): return None
def update_job(project, **changes):
    current=load_job(project) or {'schema_version':SCHEMA_VERSION,'job_id':str(uuid.uuid4()),'created_at':now()}; current.update(changes); _write(job_path(project),current)
    history_path(project).parent.mkdir(parents=True,exist_ok=True)
    with history_path(project).open('a',encoding='utf-8') as fh: fh.write(json.dumps(current,sort_keys=True,default=str)+'\n')
    return current
def create_job(project, action='RUN', **fields):
    existing=load_job(project)
    if existing and existing.get('status') in ACTIVE: raise RuntimeError('E999 PROJECT_ALREADY_RUNNING')
    record={'schema_version':SCHEMA_VERSION,'job_id':str(uuid.uuid4()),'project_id':Path(project).name,'action':action,'status':'QUEUED','created_at':now(),'started_at':None,'finished_at':None,'pid':None,'host':socket.gethostname(),'process_start_identity':None,'requested_by':'gui','resume_mode':fields.pop('resume_mode',None),'compound_id':fields.pop('compound_id',None),'failed_only':bool(fields.pop('failed_only',False)),'last_error':None,'exit_code':None,'attempt':int(existing.get('attempt',0) if existing else 0)+1,'backend_state_at_start':None,'backend_state_at_finish':None,**fields}
    _write(job_path(project),record); history_path(project).parent.mkdir(parents=True,exist_ok=True)
    with history_path(project).open('a',encoding='utf-8') as fh: fh.write(json.dumps(record,sort_keys=True)+'\n')
    return record
def list_jobs(project):
    path=history_path(project); result=[]
    if not path.exists(): return result
    for line in path.read_text(encoding='utf-8').splitlines():
        if not line.strip(): continue
        try: result.append(json.loads(line))
        except json.JSONDecodeError: continue
    return result

def _linux_process_state(pid: int) -> str | None:
    """Return the one-letter Linux process state from /proc, when available."""
    try:
        # /proc/<pid>/stat field 3 is the state. The command name in field 2 may
        # contain spaces, so split only after the final closing parenthesis.
        text=Path(f'/proc/{int(pid)}/stat').read_text(encoding='utf-8',errors='replace')
        tail=text.rsplit(')',1)[1].strip().split()
        return tail[0] if tail else None
    except (FileNotFoundError,PermissionError,OSError,ValueError,IndexError):
        return None

def process_alive(pid, identity=None):
    if not pid:
        return False
    try: pid=int(pid)
    except (TypeError,ValueError): return False
    current_identity=process_identity(pid)
    if not current_identity or (identity is not None and current_identity != str(identity)):
        return False
    # A zombie has exited and cannot perform work even though /proc/<pid> still
    # exists until its parent reaps it. Treat it as dead so GUI reconciliation
    # does not falsely report RUNNING forever.
    if _linux_process_state(pid) == 'Z':
        return False
    return True

def project_lock(project, job_id, *, blocking=False):
    return project_execution_lock(project, {'interface':'GUI','job_id':job_id})
