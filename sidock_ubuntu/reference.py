"""Deterministic chemistry-aware crystallographic reference selection."""
from __future__ import annotations
from collections import defaultdict
import csv, hashlib, json
from pathlib import Path
from typing import Any

RANKING_POLICY_VERSION='1'
RANKING_POLICY=('mapping_completeness','mapping_ambiguity','altloc_quality','occupancy_quality','rscc_rsr','binding_site_completeness','peripheral_completeness','contact_environment','identifier')

def _source(path, kind, *, parser):
    path=Path(path) if path else None
    base={'source_type':kind,'source_path':str(path) if path else None,'source_hash':None,'available':False,'parse_status':'MISSING','warnings':[]}
    if not path: return base
    if not path.exists(): base['parse_status']='MISSING'; base['warnings']=[f'{kind}_SOURCE_NOT_FOUND']; return base
    base['source_hash']=sha256_file(path)
    try:
        base.update(parser(path)); base['available']=True; base['parse_status']='PASS'
    except Exception as exc:
        base['parse_status']='INVALID'; base['warnings']=[f'{kind}_PARSE_FAILED: {exc}']
    return base

def load_local_ccd(path: Path|str|None, component_id: str|None=None) -> dict:
    """Load an offline CCD JSON (single component or component mapping)."""
    def parse(p):
        data=json.loads(p.read_text(encoding='utf-8'))
        if component_id and isinstance(data,dict) and component_id in data: data=data[component_id]
        if not isinstance(data,dict) or not (data.get('atoms') or data.get('atom_types')): raise ValueError('CCD component has no atoms')
        return {'component_id':component_id or data.get('component_id',data.get('id')),'chemistry':data}
    return _source(path,'CCD',parser=parse)

def load_local_mmcif(path: Path|str|None) -> dict:
    result=extract_mmcif_quality(Path(path)) if path and Path(path).exists() else {}
    return _source(path,'MMCIF',parser=lambda p:{'instances':extract_mmcif_quality(p),'structure_id':p.stem})

def sha256_file(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open('rb') as fh:
        for block in iter(lambda: fh.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def reference_config_fingerprint(config: dict) -> str:
    """Fingerprint only settings that can change reference selection."""
    ref=(config or {}).get('reference',{})
    return hashlib.sha256(_stable({'ranking_policy_version':RANKING_POLICY_VERSION,'ligand':ref.get('ligand'),'remote_enrichment':bool(ref.get('remote_enrichment',False)),
                                   'exclusions':ref.get('exclusions',[])}).encode()).hexdigest()

def chemical_identity_from_pdb_lines(lines: list[str]) -> tuple:
    """Compatibility helper: elemental identity, never a selection criterion."""
    return tuple(sorted((line[76:78].strip() or line[12:16].strip()[:1]).upper() for line in lines
                        if line.startswith(('ATOM','HETATM')) and (line[76:78].strip() or line[12:16].strip()[:1]).upper() not in {'H','D'}))

def extract_mmcif_quality(path: Path) -> dict[str,dict[str,Any]]:
    """Extract available local mmCIF atom/instance fields without fabricating metrics.

    Gemmi is optional at import time; an unavailable parser yields an empty
    mapping so PDB+CCD fallback remains usable.
    """
    try:
        import gemmi
        doc=gemmi.cif.read_file(str(path)); block=doc.sole_block(); out={}
        category=block.get_mmcif_category('_atom_site.')
        if not category: return out
        n=len(next(iter(category.values())))
        def val(i, tag):
            values=category.get(tag.removeprefix('_atom_site.'),[])
            return values[i] if i < len(values) else None
        for i in range(n):
            key=':'.join(str(val(i,t) or '') for t in ('_atom_site.label_comp_id','_atom_site.label_asym_id','_atom_site.auth_seq_id','_atom_site.pdbx_PDB_ins_code'))
            rec=out.setdefault(key,{'component_id':val(i,'_atom_site.label_comp_id'),'asym_id':val(i,'_atom_site.label_asym_id'),'residue_number':val(i,'_atom_site.auth_seq_id'),'insertion_code':val(i,'_atom_site.pdbx_PDB_ins_code'),'model':val(i,'_atom_site.pdbx_PDB_model_num'),'altloc':val(i,'_atom_site.label_alt_id'),'occupancies':[],'b_factors':[]})
            try: rec['occupancies'].append(float(val(i,'_atom_site.occupancy')))
            except (TypeError,ValueError): pass
            try: rec['b_factors'].append(float(val(i,'_atom_site.B_iso_or_equiv')))
            except (TypeError,ValueError): pass
        for rec in out.values():
            rec['occupancy_summary']=sum(rec['occupancies'])/len(rec['occupancies']) if rec['occupancies'] else None
            rec['b_factor_summary']=sum(rec['b_factors'])/len(rec['b_factors']) if rec['b_factors'] else None
            rec.pop('occupancies',None); rec.pop('b_factors',None); rec['rscc']=None; rec['rsr']=None
        # Deposited ligand-validation categories vary slightly between mmCIF
        # producers; accept the common spellings and leave absent values null.
        for prefix in ('_pdbx_validate_ligand.', '_pdbx_validate_ligand_atom.'):
            cat=block.get_mmcif_category(prefix)
            if not cat: continue
            n=len(next(iter(cat.values())))
            def cv(i,*names):
                for name in names:
                    vals=cat.get(name,[])
                    if i < len(vals) and vals[i] not in (None,'?','.'): return vals[i]
                return None
            for i in range(n):
                comp=cv(i,'comp_id','label_comp_id'); asym=cv(i,'asym_id','label_asym_id'); seq=cv(i,'auth_seq_id','seq_id','label_seq_id')
                for key,rec in out.items():
                    if comp and rec.get('component_id')!=comp: continue
                    if asym and rec.get('asym_id')!=asym: continue
                    if seq and str(rec.get('residue_number'))!=str(seq): continue
                    try: rec['rscc']=float(cv(i,'RSCC','rscc')) if cv(i,'RSCC','rscc') is not None else rec['rscc']
                    except (TypeError,ValueError): pass
                    try: rec['rsr']=float(cv(i,'RSR','rsr')) if cv(i,'RSR','rsr') is not None else rec['rsr']
                    except (TypeError,ValueError): pass
        return out
    except Exception:
        return {}

def _stable(x): return json.dumps(x,sort_keys=True,separators=(',',':'),default=str)

def canonical_chemical_graph(chemistry: dict[str,Any]|None) -> dict[str,Any]:
    chemistry=chemistry or {}; atoms=chemistry.get('atoms',chemistry.get('atom_types',[])); bonds=chemistry.get('bonds',[])
    aa=[]
    for a in atoms:
        if isinstance(a,str): aa.append({'id':a,'element':a[:1].upper()})
        else: aa.append({k:a.get(k) for k in ('id','element','charge','formal_charge','aromatic','stereo') if k in a})
    aa.sort(key=lambda a:str(a.get('id',''))); bb=[]
    for b in bonds:
        if isinstance(b,(list,tuple)): x,y,*rest=b; item={'a':min(str(x),str(y)),'b':max(str(x),str(y)),'order':rest[0] if rest else 1}
        else:
            x,y=b.get('a',b.get('atom1')),b.get('b',b.get('atom2')); item={'a':min(str(x),str(y)),'b':max(str(x),str(y)),'order':b.get('order',1)}
            for k in ('aromatic','stereo'):
                if k in b: item[k]=b[k]
        bb.append(item)
    bb.sort(key=lambda b:(b['a'],b['b'],str(b.get('order')))); return {'atoms':aa,'bonds':bb}

def chemical_graph_id(chemistry=None, *, graph=None):
    return 'CG-'+hashlib.sha256(_stable(canonical_chemical_graph(graph or chemistry)).encode()).hexdigest()[:16]

def map_instance_to_ccd(instance:dict, ccd:dict|None) -> dict:
    ccd=ccd or {}; observed=instance.get('atoms',instance.get('atom_names',[])); ref=ccd.get('atoms',ccd.get('atom_types',[]))
    names=lambda xs:{str(x.get('name',x.get('id',x))) if isinstance(x,dict) else str(x) for x in xs}
    obs,refids=names(observed),names(ref); mapped=instance.get('mapped_heavy_atoms',len(obs&refids) if refids else instance.get('heavy_atom_count',len(obs))); total=instance.get('total_ccd_heavy_atoms',len(refids)) or instance.get('heavy_atom_count',mapped)
    return {**instance,'ccd_available':bool(ccd),'ccd_id':instance.get('ccd_id',instance.get('ligand_code')),'mapped_heavy_atoms':int(mapped),'total_ccd_heavy_atoms':int(total),'mapping_fraction':float(mapped)/total if total else 0.0,'missing_heavy_atoms':int(instance.get('missing_heavy_atoms',max(0,total-mapped))),'extra_heavy_atoms':int(instance.get('extra_heavy_atoms',max(0,len(obs-refids)) if refids else 0)),'ambiguous_mapping':bool(instance.get('ambiguous_mapping',False)),'altloc_present':bool(instance.get('altloc_present',instance.get('alternate_location',False))),'chemical_group_id':instance.get('chemical_group_id') or (chemical_graph_id(ccd) if ccd else 'UNKNOWN')}

def group_equivalent_instances(candidates:list[dict]) -> dict[str,list[dict]]:
    groups=defaultdict(list)
    for row in candidates:
        gid=row.get('chemical_group_id'); chemistry=row.get('ccd_graph') or row.get('chemistry')
        if not gid: gid=chemical_graph_id(chemistry) if chemistry else ':'.join(map(str,(row.get('ligand_code',row.get('resname','')),row.get('heavy_atom_count',''),row.get('formal_charge','UNKNOWN'),row.get('chemical_fingerprint','UNKNOWN'))))
        groups[str(gid)].append({**row,'chemical_group_id':str(gid)})
    return dict(groups)

def _rank_key(r):
    fraction=r.get('mapping_fraction'); fraction=1.0 if fraction is None and r.get('CCD_mapping_complete') is True else (fraction or 0)
    occ=r.get('occupancy_summary',r.get('occupancy')) or 0; rscc=r.get('rscc',r.get('RSCC')); rsr=r.get('rsr',r.get('RSR')); local=r.get('local_binding_site_status',r.get('protein_contacts_5A',0)) or 0
    binding_bad=int(r.get('binding_site_incomplete_residue_count',0) or 0)+int(r.get('binding_site_incomplete_atom_count',0) or 0)
    peripheral_bad=int(r.get('peripheral_incomplete_residue_count',0) or 0)+int(r.get('peripheral_incomplete_atom_count',0) or 0)
    contact_fraction=r.get('contact_environment_complete_fraction')
    contact_fraction=1.0 if contact_fraction is None and r.get('contact_environment_complete') is True else (contact_fraction if contact_fraction is not None else 0.0)
    return (-float(fraction),bool(r.get('ambiguous_mapping',False)),bool(r.get('altloc_present',r.get('alternate_location',False))),-float(occ),-(float(rscc) if rscc is not None else -1e9),float(rsr) if rsr is not None else 1e9,binding_bad,peripheral_bad,-float(contact_fraction),-float(local),str(r.get('key',r.get('instance_id',''))))

def rank_instances(candidates:list[dict]): return sorted((dict(x) for x in candidates),key=_rank_key)

def confidence_for_candidates(candidates:list[dict]) -> str:
    groups=group_equivalent_instances(candidates)
    if len(groups)!=1: return 'LOW'
    if len(candidates)==1: return 'HIGH'
    quality=[x for x in candidates if float(x.get('mapping_fraction',1) or 0)>=1 and not x.get('ambiguous_mapping',False)]
    return 'HIGH' if len(quality)==1 else 'MODERATE'

def select_reference_instance(candidates:list[dict], override:str|None=None) -> dict:
    if not candidates: raise RuntimeError('E101 REFERENCE_LIGAND_UNRESOLVED')
    confidence=confidence_for_candidates(candidates)
    if override:
        selected=next((x for x in candidates if x.get('key')==override or x.get('instance_id')==override),None)
        if selected is None: raise RuntimeError('E101 INVALID_REFERENCE_OVERRIDE')
        # Missing CCD mapping is unknown, not a measured zero. A user may
        # explicitly select an experimental ligand whose CCD file was not
        # supplied; retain that lack of chemistry verification in metadata.
        mapping_fraction=selected.get('mapping_fraction')
        if mapping_fraction is not None:
            try:
                fraction=float(mapping_fraction)
            except (TypeError, ValueError):
                raise RuntimeError('E101 INVALID_REFERENCE_OVERRIDE: malformed CCD mapping fraction')
            if not 0 < fraction <= 1:
                raise RuntimeError('E101 INVALID_REFERENCE_OVERRIDE: invalid CCD mapping')
        selected={**selected,'selection_method':'EXPLICIT_OVERRIDE','selected':True}
    else:
        if confidence=='LOW': raise RuntimeError('E101 REFERENCE_SELECTION_LOW_CONFIDENCE')
        selected={**rank_instances(candidates)[0],'selection_method':'AUTOMATIC','selected':True}
    return {'selected':selected,'candidates':rank_instances(candidates),'confidence':confidence,'chemical_group_count':len(group_equivalent_instances(candidates))}

def write_reference_selection(project:Path,candidates:list[dict],*,override=None,input_hashes=None,provenance=None):
    selection=select_reference_instance(candidates,override); out=project/'validation'; out.mkdir(parents=True,exist_ok=True)
    fields=['component_id','ligand_code','chain_id','chain','residue_number','insertion_code','model','altloc','chemical_group_id','ccd_available','mapped_heavy_atoms','total_ccd_heavy_atoms','mapping_fraction','missing_heavy_atoms','extra_heavy_atoms','ambiguous_mapping','altloc_present','occupancy_summary','b_factor_summary','rscc','rsr','binding_site_complete','binding_site_incomplete_residue_count','binding_site_incomplete_atom_count','peripheral_incomplete_residue_count','peripheral_incomplete_atom_count','contact_environment_complete','contact_environment_complete_fraction','local_binding_site_status','local_peripheral_status','ranking_position','selection_confidence','selected','selection_reason','key']
    csv_path=out/'reference_ligand_candidates.csv'
    with csv_path.open('w',newline='',encoding='utf-8') as fh:
        w=csv.DictWriter(fh,fieldnames=fields,extrasaction='ignore'); w.writeheader()
        for i,row in enumerate(selection['candidates'],1):
            reason='explicit override' if override and row.get('key')==override else ('complete mapping/no ambiguity/highest local quality' if i==1 else 'candidate retained for audit')
            w.writerow({**row,'ranking_position':i,'selection_confidence':selection['confidence'],'selected':row.get('key')==selection['selected'].get('key'),'selection_reason':reason})
    payload={'schema_version':'0.4-reference-selection-1','reference_ranking_policy_version':RANKING_POLICY_VERSION,'selection_method':selection['selected']['selection_method'],'selected_instance':selection['selected'],'chemical_group_id':selection['selected'].get('chemical_group_id'),'confidence':selection['confidence'],'candidate_count':len(candidates),'chemical_group_count':selection['chemical_group_count'],'ranking_policy':list(RANKING_POLICY),'candidate_ranking':selection['candidates'],'selection_reason':('explicit override' if override else 'highest deterministic pre-docking quality'),'warnings':[],'override':override,'input_hashes':input_hashes or {},'provenance':provenance or {}}
    json_path=out/'reference_selection.json'; json_path.write_text(json.dumps(payload,indent=2,sort_keys=True,default=str),encoding='utf-8'); return {**payload,'candidate_csv':str(csv_path),'selection_json':str(json_path)}
