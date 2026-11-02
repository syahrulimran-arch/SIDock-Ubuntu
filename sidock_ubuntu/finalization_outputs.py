from __future__ import annotations
from dataclasses import dataclass, asdict, field
from pathlib import Path
import hashlib, json

@dataclass
class FinalizationStageOutput:
    results_csv: str
    evidence_json: str
    report_html: str
    project_status: str
    warnings: list[str] = field(default_factory=list)
    artifacts: list[dict] = field(default_factory=list)

def _hash(p: Path):
    h=hashlib.sha256();
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def save_finalization_output(project: Path, value: FinalizationStageOutput) -> Path:
    out=project/'validation/finalization_stage.json'; out.parent.mkdir(parents=True,exist_ok=True)
    payload={'schema_version':'0.4-finalization-1','result':asdict(value)}
    tmp=out.with_suffix('.json.tmp'); tmp.write_text(json.dumps(payload,indent=2,sort_keys=True),encoding='utf-8'); tmp.replace(out)
    return out

def _load(project: Path):
    path=project/'validation/finalization_stage.json'
    if not path.exists(): raise FileNotFoundError(path)
    data=json.loads(path.read_text(encoding='utf-8')); result=data.get('result',data)
    value=FinalizationStageOutput(**{k:result[k] for k in FinalizationStageOutput.__dataclass_fields__ if k in result})
    issues=[]
    for item in value.artifacts:
        p=project/item['path']
        if not p.exists(): issues.append(('MISSING',str(p)))
        elif _hash(p)!=item['sha256']: issues.append(('STALE',str(p)))
    value.load_status='VALID' if not issues else issues[0][0]; value.load_issues=issues
    return value

def load_finalization_output(project: Path): return _load(project)
