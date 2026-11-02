"""Recovered PLIP workflow from the preserved production dry-run provenance.

The commands are intentionally kept identical to the recorded workflow:
``obabel <complex.pdb> -O <protonated.pdb> -h`` followed by
``plip -f <protonated.pdb> --nohydro -o <out> -x -q``.
"""
from __future__ import annotations
import csv, subprocess, xml.etree.ElementTree as ET
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

def prepare_protonated_complex(complex_pdb: Path, destination: Path, *, runner=subprocess.run) -> Path:
    destination.parent.mkdir(parents=True,exist_ok=True)
    p=runner(['obabel',str(complex_pdb),'-O',str(destination),'-h'],capture_output=True,text=True)
    if p.returncode != 0 or not destination.exists(): raise SIDockError('PLIP_FAILED','Open Babel protonation failed',{'command':['obabel',str(complex_pdb),'-O',str(destination),'-h']})
    return destination

def run_plip(protonated: Path, output_dir: Path, *, runner=subprocess.run) -> list[Path]:
    output_dir.mkdir(parents=True,exist_ok=True)
    resolved=resolve_plip_executable()
    if not resolved['available']:
        raise SIDockError('PLIP_FAILED',resolved['message'],{})
    command=[resolved['path'],'-f',str(protonated),'--nohydro','-o',str(output_dir),'-x','-q']
    p=runner(command,capture_output=True,text=True)
    reports=list(output_dir.rglob('*_report.xml'))
    if p.returncode != 0 or not reports: raise SIDockError('PLIP_FAILED','PLIP failed or produced no XML report',{'command':command})
    return reports

def _text(node, name):
    x=node.find(name); return x.text.strip() if x is not None and x.text else ''

def normalize_plip_interactions(reports: list[Path], destination: Path, *, compound_id='', pose_id='') -> Path:
    rows=[]
    type_map={'hydrogen_bonds':'hydrogen_bond','hydrophobic_interactions':'hydrophobic','pistacking':'pi-stack','pi_cation_interactions':'pi-cation','salt_bridges':'salt_bridge','halogen_bonds':'halogen_bond'}
    for report in reports:
        root=ET.parse(report).getroot()
        for section, interaction_type in type_map.items():
            parent=root.find(f'.//{section}')
            if parent is None: continue
            for item in list(parent):
                chain=_text(item,'reschain'); resname=_text(item,'restype'); resnum=_text(item,'resnr')
                rows.append({'compound_id':compound_id,'pose_id':pose_id,'interaction_type':interaction_type,'protein_chain':chain,'protein_residue_name':resname,'protein_residue_number':resnum,'protein_atom':_text(item,'protisdon') or _text(item,'protisacc'),'ligand_atom':_text(item,'ligisdon') or _text(item,'ligisacc'),'distance_A':_text(item,'dist'),'angle_deg':_text(item,'angle'),'PLIP_source':'PLIP'})
    destination.parent.mkdir(parents=True,exist_ok=True)
    with destination.open('w',newline='',encoding='utf-8') as fh:
        w=csv.DictWriter(fh,fieldnames=SCHEMA); w.writeheader(); w.writerows(rows)
    return destination

def run_plip_analysis(complex_pdb: Path, work_dir: Path, normalized_csv: Path, *, compound_id='', pose_id='', provenance_path=None, runner=subprocess.run) -> PlipStageOutput:
    protonated=work_dir/f'{complex_pdb.stem}_protonated.pdb'; protonated=prepare_protonated_complex(complex_pdb,protonated,runner=runner)
    reports=run_plip(protonated,work_dir/'plip',runner=runner)
    normalize_plip_interactions(reports,normalized_csv,compound_id=compound_id,pose_id=pose_id)
    return PlipStageOutput(protonated,reports,normalized_csv,provenance_path,[])
