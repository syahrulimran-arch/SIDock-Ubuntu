"""Persisted, restart-safe contracts for global production stages."""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "0.4-global-stage-1"

class StageOutputError(ValueError):
    """Persisted stage metadata is missing required structure."""

@dataclass
class StageArtifact:
    path: str
    sha256: str
    required: bool = True

@dataclass
class ReferenceStageOutput:
    selected_ligand: str
    residue_name: str
    chain: str
    residue_number: str
    reference_structure: str
    reference_chemistry: str | None = None
    ccd_id: str | None = None
    selection_metadata: str | None = None
    warnings: list[str] = field(default_factory=list)
    artifacts: list[StageArtifact] = field(default_factory=list)
    source_metadata: dict[str, Any] = field(default_factory=dict)
    input_hashes: dict[str, str] = field(default_factory=dict)

@dataclass
class ReceptorStageOutput:
    prepared_receptor_pdb: str
    prepared_receptor_pdbqt: str
    completeness_report: str | None = None
    fallback_metadata: str | None = None
    retained_chains: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    artifacts: list[StageArtifact] = field(default_factory=list)

@dataclass
class GridStageOutput:
    gpf: str
    glg: str
    maps: list[str] = field(default_factory=list)
    center: list[float] = field(default_factory=list)
    dimensions: list[int] = field(default_factory=list)
    spacing: float | None = None
    warnings: list[str] = field(default_factory=list)
    artifacts: list[StageArtifact] = field(default_factory=list)

@dataclass
class ValidationStageOutput:
    native_dlg: str | None
    minimum_rmsd: float | None
    median_rmsd: float | None = None
    lowest_energy_rmsd: float | None = None
    rmsd_threshold: float | None = None
    scientific_status: str = "NOT_APPLICABLE"
    execution_status: str = "SUCCESS"
    native_like_pose_count: int | None = None
    artifacts: list[StageArtifact] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

def _hash(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(1024*1024), b''): h.update(block)
    return h.hexdigest()

def _rel(project: Path, value: Path | str | None) -> str | None:
    if value is None: return None
    p=Path(value)
    try: return str(p.resolve().relative_to(project.resolve()))
    except ValueError: return str(p)

def artifact(project: Path, path: Path | str, required: bool=True) -> StageArtifact:
    p=Path(path); return StageArtifact(_rel(project,p) or str(p), _hash(p), required) if p.exists() else StageArtifact(_rel(project,p) or str(p), '', required)

def _write(project: Path, name: str, value: Any) -> Path:
    out=project/name; out.parent.mkdir(parents=True,exist_ok=True)
    payload={'schema_version':SCHEMA_VERSION,'stage':name.split('/')[-1].replace('_stage.json','').upper(),'result':asdict(value)}
    tmp=out.with_suffix(out.suffix+'.tmp'); tmp.write_text(json.dumps(payload,indent=2,sort_keys=True),encoding='utf-8'); tmp.replace(out)
    return out

def save_reference_output(project: Path, value: ReferenceStageOutput) -> Path: return _write(project,'validation/reference_stage.json',value)
def save_receptor_output(project: Path, value: ReceptorStageOutput) -> Path: return _write(project,'validation/receptor_stage.json',value)
def save_grid_output(project: Path, value: GridStageOutput) -> Path: return _write(project,'validation/grid_stage.json',value)
def save_validation_output(project: Path, value: ValidationStageOutput) -> Path: return _write(project,'validation/validation_stage.json',value)

def _load(project: Path, name: str, cls):
    path=project/name
    if not path.exists(): raise FileNotFoundError(f'MISSING_STAGE_METADATA: {path}')
    data=json.loads(path.read_text(encoding='utf-8'))
    result=data.get('result',data)
    result=dict(result)
    required={ReferenceStageOutput: ('selected_ligand','residue_name','chain','residue_number','reference_structure'),
              ReceptorStageOutput: ('prepared_receptor_pdb','prepared_receptor_pdbqt'),
              GridStageOutput: ('gpf','glg'),
              ValidationStageOutput: ('minimum_rmsd','rmsd_threshold','scientific_status','execution_status')}.get(cls,())
    missing=[key for key in required if key not in result]
    if missing: raise StageOutputError(f'INCOMPLETE_STAGE_METADATA: {name}: {", ".join(missing)}')
    result['artifacts']=[StageArtifact(**x) if isinstance(x,dict) else x for x in result.get('artifacts',[])]
    value=cls(**{k:v for k,v in result.items() if k in cls.__dataclass_fields__})
    statuses=[]
    for item in value.artifacts:
        p=project/item.path if not Path(item.path).is_absolute() else Path(item.path)
        if item.required and not p.exists(): statuses.append(('MISSING',str(p)))
        elif p.exists() and item.sha256 and _hash(p)!=item.sha256: statuses.append(('STALE',str(p)))
    value.load_status='VALID' if not statuses else statuses[0][0]
    value.load_issues=statuses
    return value

def load_reference_output(project: Path) -> ReferenceStageOutput: return _load(project,'validation/reference_stage.json',ReferenceStageOutput)
def load_receptor_output(project: Path) -> ReceptorStageOutput: return _load(project,'validation/receptor_stage.json',ReceptorStageOutput)
def load_grid_output(project: Path) -> GridStageOutput: return _load(project,'validation/grid_stage.json',GridStageOutput)
def load_validation_output(project: Path) -> ValidationStageOutput: return _load(project,'validation/validation_stage.json',ValidationStageOutput)
