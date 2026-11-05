"""PLIP workflow and PDB-complex materialization.

PLIP is run as ``obabel <complex.pdb> -O <protonated.pdb> -h`` followed by
``plip -f <protonated.pdb> --nohydro -o <out> -x -q``.  Docking pose
artifacts are PDBQT, so they are first converted into a real PDB complex in
which the docked compound is an explicit HETATM residue.  This prevents a
successful PLIP invocation from silently treating the ligand as protein.
"""
from __future__ import annotations
import csv, shutil, subprocess, xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from .errors import SIDockError
from .dependencies import resolve_plip_executable

@dataclass
class PlipStageOutput:
    protonated_complex_path: Path
    raw_output_paths: list[Path]
    normalized_interactions_path: Path
    provenance_path: Path | None = None
    warnings: list[str] | None = None

SCHEMA=['compound_id','pose_id','interaction_type','protein_chain','protein_residue_name','protein_residue_number','protein_atom','ligand_atom','distance_A','angle_deg','PLIP_source']

def _obabel_convert(source: Path, destination: Path, *, runner=subprocess.run) -> Path:
    destination.parent.mkdir(parents=True,exist_ok=True)
    p=runner(['obabel',str(source),'-O',str(destination)],capture_output=True,text=True)
    destination.with_suffix('.conversion.log').write_text(
        str(getattr(p, 'stdout', '') or '') + str(getattr(p, 'stderr', '') or ''), encoding='utf-8')
    if p.returncode != 0 or not destination.exists():
        raise SIDockError('PLIP_FAILED','Open Babel format conversion failed',{'command':['obabel',str(source),'-O',str(destination)]})
    return destination

def materialize_plip_complex(source: Path, destination: Path, *, runner=subprocess.run) -> Path:
    """Return a PDB complex with the docked ligand explicitly marked HETATM.

    ``extract_docking_poses`` persists a receptor PDBQT followed by a docked
    ligand MODEL.  PLIP requires a chemically identifiable hetero residue, so
    the two regions are converted separately and recombined with a unique
    ligand residue identity (LIG, chain Z, residue 999).  Ligand CONECT records
    emitted by Open Babel are renumbered together with the ligand atoms.
    """
    source=Path(source); destination=Path(destination)
    destination.parent.mkdir(parents=True,exist_ok=True)
    if source.suffix.lower()=='.pdb':
        if source.resolve()!=destination.resolve(): shutil.copy2(source,destination)
        return destination
    if source.suffix.lower()!='.pdbqt':
        raise SIDockError('PLIP_FAILED',f'Unsupported PLIP complex format: {source.suffix}',{})

    lines=source.read_text(encoding='utf-8',errors='replace').splitlines()
    model_index=next((i for i,line in enumerate(lines) if line.startswith('MODEL')),None)
    if model_index is None:
        raise SIDockError('PLIP_FAILED','Docked PDBQT complex does not contain a ligand MODEL',{})
    receptor_lines=lines[:model_index]
    ligand_lines=lines[model_index:]
    if not any(line.startswith(('ATOM','HETATM')) for line in receptor_lines):
        raise SIDockError('PLIP_FAILED','PDBQT complex contains no receptor atoms',{})
    if not any(line.startswith(('ATOM','HETATM')) for line in ligand_lines):
        raise SIDockError('PLIP_FAILED','PDBQT complex contains no ligand atoms',{})

    scratch=destination.parent/'plip_input_parts'; scratch.mkdir(parents=True,exist_ok=True)
    receptor_pdbqt=scratch/'receptor.pdbqt'; ligand_pdbqt=scratch/'ligand.pdbqt'
    receptor_pdbqt.write_text('\n'.join(receptor_lines)+'\n',encoding='utf-8')
    ligand_pdbqt.write_text('\n'.join(ligand_lines)+'\n',encoding='utf-8')
    receptor_pdb=_obabel_convert(receptor_pdbqt,scratch/'receptor.pdb',runner=runner)
    ligand_pdb=_obabel_convert(ligand_pdbqt,scratch/'ligand.pdb',runner=runner)

    receptor_out=[]; max_serial=0
    for line in receptor_pdb.read_text(encoding='utf-8',errors='replace').splitlines():
        if line.startswith(('ATOM  ','HETATM')):
            if (line[17:20].strip(), line[21:22].strip(), line[22:26].strip()) == ('LIG', 'Z', '999'):
                raise SIDockError('PLIP_FAILED', 'Receptor already uses reserved docked-ligand identity LIG:Z:999', {})
            receptor_out.append(line.rstrip())
            try: max_serial=max(max_serial,int(line[6:11]))
            except (ValueError,IndexError): pass
        elif line.startswith('TER'):
            receptor_out.append(line.rstrip())

    ligand_raw=ligand_pdb.read_text(encoding='utf-8',errors='replace').splitlines()
    mapping={}; next_serial=max_serial+1
    for line in ligand_raw:
        if not line.startswith(('ATOM  ','HETATM')): continue
        try: old=int(line[6:11])
        except (ValueError,IndexError): continue
        mapping[old]=next_serial; next_serial+=1
    if not mapping:
        raise SIDockError('PLIP_FAILED','Open Babel produced no ligand atoms for PLIP',{})

    ligand_out=[]; conect=[]
    for line in ligand_raw:
        if line.startswith(('ATOM  ','HETATM')):
            try: old=int(line[6:11]); new=mapping[old]
            except (ValueError,IndexError,KeyError): continue
            chars=list(line.rstrip('\n').ljust(80))
            chars[0:6]=list('HETATM')
            chars[6:11]=list(f'{new:5d}')
            chars[16]=' '
            chars[17:20]=list('LIG')
            chars[21]='Z'
            chars[22:26]=list(f'{999:4d}')
            chars[26]=' '
            ligand_out.append(''.join(chars).rstrip())
        elif line.startswith('CONECT'):
            ids=[]
            for token in line.split()[1:]:
                try: ids.append(int(token))
                except ValueError: pass
            remapped=[mapping[x] for x in ids if x in mapping]
            if len(remapped)>=2:
                conect.append('CONECT'+''.join(f'{x:5d}' for x in remapped))

    destination.write_text('\n'.join(receptor_out+['TER']+ligand_out+conect+['END'])+'\n',encoding='utf-8')
    return destination

def prepare_protonated_complex(complex_pdb: Path, destination: Path, *, runner=subprocess.run) -> Path:
    destination.parent.mkdir(parents=True,exist_ok=True)
    p=runner(['obabel',str(complex_pdb),'-O',str(destination),'-h'],capture_output=True,text=True)
    destination.with_suffix('.protonation.log').write_text(
        str(getattr(p, 'stdout', '') or '') + str(getattr(p, 'stderr', '') or ''), encoding='utf-8')
    if p.returncode != 0 or not destination.exists(): raise SIDockError('PLIP_FAILED','Open Babel protonation failed',{'command':['obabel',str(complex_pdb),'-O',str(destination),'-h']})
    return destination

def run_plip(protonated: Path, output_dir: Path, *, runner=subprocess.run) -> list[Path]:
    output_dir.mkdir(parents=True,exist_ok=True)
    resolved=resolve_plip_executable()
    if not resolved['available']:
        raise SIDockError('PLIP_FAILED',resolved['message'],{})
    command=[resolved['path'],'-f',str(protonated),'--nohydro','-o',str(output_dir),'-x','-q']
    p=runner(command,capture_output=True,text=True)
    (output_dir / 'plip.log').write_text(
        str(getattr(p, 'stdout', '') or '') + str(getattr(p, 'stderr', '') or ''), encoding='utf-8')
    reports=sorted(set(output_dir.rglob('*_report.xml')) | set(output_dir.rglob('report.xml')))
    if p.returncode != 0 or not reports: raise SIDockError('PLIP_FAILED','PLIP failed or produced no XML report',{'command':command})
    return reports

def _text(node, name):
    x=node.find(name); return x.text.strip() if x is not None and x.text else ''


def _first_text(node, *names):
    return next((value for name in names if (value := _text(node, name))), '')


def _interaction_distance(item, kind):
    # Use the physical distance definition encoded by PLIP for this contact.
    # Never substitute H...A for donor...acceptor without an explicit label.
    if kind == 'hydrogen_bond':
        return _first_text(item, 'dist_d-a', 'dist_d_a', 'dist_d.a', 'dist')
    if kind == 'pi-stack':
        return _first_text(item, 'centdist', 'cent_dist', 'cent-dist', 'dist')
    return _first_text(item, 'dist')


def _interaction_atom_indices(item, kind):
    """Return PLIP atom *indices* only when their protein/ligand role is explicit."""
    if kind == 'hydrophobic':
        return _first_text(item, 'protcarbonidx'), _first_text(item, 'ligcarbonidx')
    if kind == 'hydrogen_bond':
        protein_is_donor = _text(item, 'protisdon').strip().lower()
        donor = _first_text(item, 'donoridx', 'don_idx')
        acceptor = _first_text(item, 'acceptoridx', 'acc_idx')
        if protein_is_donor in {'true', '1', 'yes'}:
            return donor, acceptor
        if protein_is_donor in {'false', '0', 'no'}:
            return acceptor, donor
    return '', ''


def normalize_plip_interactions(reports: list[Path], destination: Path, *, compound_id='', pose_id='', ligand_identity=None) -> Path:
    rows=[]
    type_map={'hydrogen_bonds':'hydrogen_bond','hydrophobic_interactions':'hydrophobic','pistacking':'pi-stack','pi_cation_interactions':'pi-cation','salt_bridges':'salt_bridge','halogen_bonds':'halogen_bond'}
    matched_site = False
    for report in reports:
        root=ET.parse(report).getroot()
        if ligand_identity is not None:
            sites = [site for site in root.findall('.//bindingsite')
                     if site.find('identifiers') is not None
                     if tuple(_text(site.find('identifiers'), key) for key in ('hetid', 'chain', 'position')) == tuple(ligand_identity)]
            matched_site = matched_site or bool(sites)
        else:
            sites = [root]
        for section, interaction_type in type_map.items():
            for parent in (parent for site in sites for parent in site.findall(f'.//{section}')):
                for item in list(parent):
                    chain=_text(item,'reschain'); resname=_text(item,'restype'); resnum=_text(item,'resnr')
                    protein_idx, ligand_idx = _interaction_atom_indices(item, interaction_type)
                    rows.append({'compound_id':compound_id,'pose_id':pose_id,'interaction_type':interaction_type,'protein_chain':chain,'protein_residue_name':resname,'protein_residue_number':resnum,'protein_atom':protein_idx,'ligand_atom':ligand_idx,'distance_A':_interaction_distance(item, interaction_type),'angle_deg':_first_text(item,'don_angle','angle'),'PLIP_source':'PLIP'})
    if ligand_identity is not None and not matched_site:
        raise SIDockError('PLIP_FAILED', 'PLIP did not identify the docked LIG:Z:999 residue; contacts are unavailable.', {})
    destination.parent.mkdir(parents=True,exist_ok=True)
    with destination.open('w',newline='',encoding='utf-8') as fh:
        w=csv.DictWriter(fh,fieldnames=SCHEMA); w.writeheader(); w.writerows(rows)
    return destination

def run_plip_analysis(complex_pdb: Path, work_dir: Path, normalized_csv: Path, *, compound_id='', pose_id='', provenance_path=None, runner=subprocess.run, ligand_identity=None) -> PlipStageOutput:
    work_dir.mkdir(parents=True,exist_ok=True)
    plip_input=materialize_plip_complex(complex_pdb,work_dir/'plip_input_complex.pdb',runner=runner)
    protonated=work_dir/f'{plip_input.stem}_protonated.pdb'; protonated=prepare_protonated_complex(plip_input,protonated,runner=runner)
    reports=run_plip(protonated,work_dir/'plip',runner=runner)
    normalize_plip_interactions(reports,normalized_csv,compound_id=compound_id,pose_id=pose_id,
                               ligand_identity=ligand_identity)
    return PlipStageOutput(protonated,reports,normalized_csv,provenance_path,[])
