"""Backend-neutral advisory lock for scientific project writers."""
from __future__ import annotations
import fcntl, json, os, socket
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

LOCK_NAME='provenance/project_execution.lock'
class ProjectExecutionLocked(BlockingIOError):
    """Operational conflict: another process owns the project lock."""
    pass
def lock_path(project): return Path(project)/LOCK_NAME
def process_identity(pid):
    try: return Path(f'/proc/{pid}/stat').read_text().split()[21]
    except (FileNotFoundError,IndexError): return None
@contextmanager
def project_execution_lock(project, owner_metadata=None):
    path=lock_path(project); path.parent.mkdir(parents=True,exist_ok=True); fh=path.open('a+')
    try:
        try: fcntl.flock(fh.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as exc: raise ProjectExecutionLocked('PROJECT_EXECUTION_LOCKED') from exc
        metadata={'pid':os.getpid(),'host':socket.gethostname(),'acquired_at':datetime.now(timezone.utc).isoformat(),**(owner_metadata or {})}
        fh.seek(0); fh.truncate(); fh.write(json.dumps(metadata)); fh.flush(); yield fh
    finally:
        try: fcntl.flock(fh.fileno(),fcntl.LOCK_UN)
        finally: fh.close()
