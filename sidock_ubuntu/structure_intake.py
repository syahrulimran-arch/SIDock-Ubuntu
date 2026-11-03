"""Deterministic semantic classification of uploaded structure files."""
from __future__ import annotations
from dataclasses import asdict, dataclass
from pathlib import Path
import json

WATERS = {'HOH','WAT','H2O','DOD'}
IONS = {'CA','ZN','MG','MN','FE','NA','K','CL','CU','CO','NI','CD','HG','BR','I','LI','CS','SR','BA','AL','AG','AU','PB'}
AMINO = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','ASX','GLX','SEC','PYL'}
NUCLEOTIDES = {'A','C','G','U','DA','DC','DG','DT'}

@dataclass(frozen=True)
class StructureClassification:
    classification: str
    candidates: tuple[dict, ...]
    ions: tuple[dict, ...]
    waters: tuple[dict, ...]
    excluded: tuple[dict, ...]
    processed_indicators: tuple[str, ...]
    reason: str

def _records(path: Path):
    grouped={}
    for line in path.read_text(encoding='utf-8', errors='replace').splitlines():
        if line.startswith(('ATOM','HETATM')):
            name=line[17:20].strip().upper(); chain=line[21:22].strip(); resseq=line[22:27].strip()
            grouped.setdefault((name,chain,resseq),[]).append(line)
    return grouped

def classify_structure(path: Path, *, mmcif: Path | None = None) -> StructureClassification:
    groups=_records(Path(path)); candidates=[]; ions=[]; waters=[]; excluded=[]
    for (name,chain,resseq), lines in groups.items():
        heavy=sum(1 for x in lines if (x[76:78].strip().upper() or x[12:16].strip()[:1].upper()) not in {'H','D'})
        item={'component':name,'chain':chain,'residue_number':resseq,'heavy_atoms':heavy,'record_types':sorted({x[:6].strip() for x in lines})}
        if name in WATERS: waters.append({**item,'reason':'water'})
        elif name in IONS and heavy <= 2: ions.append({**item,'reason':'inorganic ion'})
        elif name in AMINO or name in NUCLEOTIDES: excluded.append({**item,'reason':'standard polymer residue'})
        elif heavy >= 3: candidates.append(item)
        else: excluded.append({**item,'reason':'too small for reference ligand'})
    text=Path(path).read_text(encoding='utf-8', errors='replace').upper()
    indicators=[]
    if any(x in text for x in ('DISCOVERY STUDIO','ACCELRYS')): indicators.append('ACCELRYS_DISCOVERY_STUDIO_REMARK')
    if not any(line.startswith(('HEADER','TITLE','COMPND','HETNAM')) for line in text.splitlines()): indicators.append('MISSING_DEPOSITED_METADATA')
    explicit_h=sum(1 for line in text.splitlines() if line.startswith(('ATOM','HETATM')) and (line[76:78].strip() or line[12:16].strip()[:1]) in {'H','D'})
    if explicit_h: indicators.append('EXPLICIT_HYDROGENS')
    if candidates:
        classification='HOLO_NATIVE'
        reason='organic nonpolymer reference candidate present'
    elif ions:
        classification='ION_ONLY'
        reason='only inorganic ion components were detected'
    elif indicators or not groups:
        classification='PROCESSED_STRIPPED' if indicators else 'APO_NO_REFERENCE'
        reason='no crystallographic organic reference ligand is present'
    else:
        classification='APO_NO_REFERENCE'; reason='no experimental reference ligand is present'
    if candidates and indicators: classification='PROCESSED_STRIPPED' if not any(x['component'] not in IONS for x in candidates) else classification
    return StructureClassification(classification,tuple(candidates),tuple(ions),tuple(waters),tuple(excluded),tuple(indicators),reason)

def save_classification(value: StructureClassification, destination: Path) -> Path:
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(asdict(value),indent=2),encoding='utf-8'); return destination
