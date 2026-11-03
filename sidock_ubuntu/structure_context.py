"""Production/validation structure provenance and equivalence checks."""
from __future__ import annotations
from dataclasses import asdict, dataclass
from pathlib import Path
import hashlib, json, math

AA=set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())

def _hash(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def _residues(path):
    result={}
    for line in Path(path).read_text(errors='replace').splitlines():
        if not line.startswith('ATOM') or line[17:20].strip().upper() not in AA: continue
        chain=line[21:22].strip(); number=line[22:27].strip(); key=(chain,number)
        atom=line[12:16].strip().upper(); result.setdefault(key,{'name':line[17:20].strip().upper(),'atoms':{}})['atoms'][atom]=line
    return result

def _xyz(line): return tuple(float(line[a:b]) for a,b in ((30,38),(38,46),(46,54)))

@dataclass(frozen=True)
class StructureComparison:
    classification: str
    sequence_identity: float
    mapping_coverage: float
    backbone_rmsd: float | None
    mapped_residues: int
    production_residues: int
    validation_residues: int
    mutations: tuple[dict,...]
    missing_residues: tuple[str,...]
    binding_site_mutations: tuple[dict,...]
    reason: str

def compare_validation_receptor(production: Path, validation: Path, *, production_pdb_id=None, validation_pdb_id=None, binding_site=None) -> StructureComparison:
    if Path(production).resolve() == Path(validation).resolve():
        residues=_residues(production)
        return StructureComparison('SAME_STRUCTURE',1.0,1.0,0.0,len(residues),len(residues),len(residues),(),(),(),'production and validation structure are identical')
    p=_residues(production); v=_residues(validation); pairs=[]; mutations=[]; missing=[]
    for key, rec in p.items():
        other=v.get(key)
        if other is None: missing.append(f'{key[0]}:{key[1]}'); continue
        pairs.append((rec,other))
        if rec['name'] != other['name']: mutations.append({'chain':key[0],'residue_number':key[1],'production':rec['name'],'validation':other['name']})
    coverage=len(pairs)/max(1,len(p)); identity=sum(a['name']==b['name'] for a,b in pairs)/max(1,len(pairs))
    coords=[]
    for a,b in pairs:
        if 'CA' in a['atoms'] and 'CA' in b['atoms']:
            x,y=_xyz(a['atoms']['CA']),_xyz(b['atoms']['CA']); coords.append(sum((x[i]-y[i])**2 for i in range(3)))
    rmsd=math.sqrt(sum(coords)/len(coords)) if coords else None
    site=set(binding_site or ())
    site_mut=tuple(x for x in mutations if (x['chain'],x['residue_number']) in site)
    same_id=bool(production_pdb_id and validation_pdb_id and production_pdb_id.upper()==validation_pdb_id.upper())
    if coverage < .90 or identity < .95 or rmsd is None:
        cls='INCOMPATIBLE'; reason='insufficient residue mapping or structural coordinates'
    elif site_mut:
        cls='INCOMPATIBLE'; reason='binding-site mutations detected'
    elif identity >= .995 and rmsd <= 1.5 and (same_id or not production_pdb_id or not validation_pdb_id):
        cls='SAME_STRUCTURE' if same_id else 'PREPARED_DERIVATIVE_OF_SAME_STRUCTURE'; reason='high-identity, high-coverage mapped receptor'
    else:
        cls='RELATED_BUT_DIFFERENT_STRUCTURE'; reason='related receptor is not proven to be the same deposited structure'
    return StructureComparison(cls,identity,coverage,rmsd,len(pairs),len(p),len(v),tuple(mutations),tuple(missing),site_mut,reason)

def structure_context(production: Path, validation: Path | None = None, **ids) -> dict:
    validation=Path(validation or production)
    cmp=compare_validation_receptor(production,validation,**ids)
    return {'production_receptor':{'path':str(Path(production).resolve()),'sha256':_hash(production)},'validation_structure':{'path':str(validation.resolve()),'sha256':_hash(validation)},'equivalence':asdict(cmp),'validation_mode':('NATIVE_SAME_STRUCTURE_VALIDATED' if cmp.classification=='SAME_STRUCTURE' and validation==Path(production) else ('PAIRED_STRIPPED_SAME_STRUCTURE_VALIDATED' if cmp.classification=='PREPARED_DERIVATIVE_OF_SAME_STRUCTURE' else ('SURROGATE_HOLO_VALIDATION' if cmp.classification=='RELATED_BUT_DIFFERENT_STRUCTURE' else 'NO_NATIVE_VALIDATION')))}
