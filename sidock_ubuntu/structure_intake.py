"""Deterministic semantic classification of uploaded structure files."""
from __future__ import annotations
from dataclasses import asdict, dataclass
from pathlib import Path
import json

def _mmcif_components(path: Path) -> dict:
    """Read the small, stable metadata subset needed for intake.

    This deliberately does not replace coordinate parsing.  It only supplies
    deposited component/entity identity when a local mmCIF is available.
    """
    result = {'components': [], 'source': str(Path(path).resolve())}
    try:
        import gemmi
        doc = gemmi.cif.read_file(str(path))
        block = doc.sole_block()
        for tag in ('_chem_comp.id', '_pdbx_nonpoly_scheme.mon_id', '_atom_site.label_comp_id'):
            values = block.find_values(tag)
            if values:
                result['components'].extend(str(v).strip().upper() for v in values if str(v).strip() not in {'.','?'})
        result['components'] = sorted(set(result['components']))
        return result
    except Exception:
        # Keep local legacy parsing usable when Gemmi cannot parse an optional
        # metadata file; the caller records that PDB fallback was used.
        return result

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
    mm_components = _mmcif_components(Path(mmcif))['components'] if mmcif and Path(mmcif).exists() else []
    # mmCIF component identity is authoritative for classification.  The PDB
    # remains the coordinate source; synthetic metadata-only rows are enough
    # to prevent an ambiguous PDB record from overriding deposited identity.
    if mm_components:
        pdb_by_name = {x['component']: x for x in candidates + ions + waters + excluded}
        for name in mm_components:
            if name in WATERS and name not in {x['component'] for x in waters}:
                waters.append({'component':name,'chain':'','residue_number':'','heavy_atoms':0,'record_types':[],'reason':'water'})
            elif name in IONS and name not in {x['component'] for x in ions}:
                ions.append({'component':name,'chain':'','residue_number':'','heavy_atoms':1,'record_types':[],'reason':'inorganic ion'})
            elif name not in WATERS and name not in IONS and name not in AMINO and name not in NUCLEOTIDES and name not in {x['component'] for x in candidates}:
                candidates.append({'component':name,'chain':'','residue_number':'','heavy_atoms':0,'record_types':[],'metadata_only':True})
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
