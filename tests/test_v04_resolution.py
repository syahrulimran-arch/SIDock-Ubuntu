import json
from pathlib import Path
import pytest
from sidock_ubuntu.intake import (HeteroGroup, reference_candidates, receptor_completeness,
                                  rank_reference_candidates, write_reference_candidates)
from sidock_ubuntu.compliance import evaluate
from sidock_ubuntu.state import update_state, project_status

def atom(serial, name, res, chain, num, x, y, z, record='ATOM'):
    return f"{record:<6}{serial:5d} {name:^4s} {res:>3s} {chain}{num:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           {name[0]:>2s}"

def ligand(serial, res, chain, num, x):
    return atom(serial, 'C', res, chain, num, x, 0, 0, 'HETATM')

def make_pdb(path, ligand_names=('LIG',), ligand_positions=(2.0,)):
    lines=[atom(1,'N','ALA','A',1,0,0,0),atom(2,'CA','ALA','A',1,1,0,0),atom(3,'C','ALA','A',1,1,1,0),atom(4,'O','ALA','A',1,1,1,1),atom(5,'CB','ALA','A',1,0,0,1)]
    serial=6
    for name,x in zip(ligand_names,ligand_positions):
        for j in range(8):
            lines.append(ligand(serial,name,'A',106 if name=='LIG' else 107,x+j*0.1)); serial += 1
    path.write_text('\n'.join(lines)+'\n')

def test_single_candidate_discovery(tmp_path):
    p=tmp_path/'x.pdb'; make_pdb(p)
    rows=reference_candidates(p)
    assert len(rows)==1 and rows[0]['recommended']

def test_equivalent_candidates_are_groupable_by_resname(tmp_path):
    p=tmp_path/'x.pdb'; make_pdb(p,('LIG','LIG'),(2,3))
    rows=reference_candidates(p)
    assert len(rows)==1 and rows[0]['ligand_code']=='LIG'

def test_candidate_csv_written(tmp_path):
    p=tmp_path/'x.pdb'; make_pdb(p); out=tmp_path/'c.csv'; rows=write_reference_candidates(p,out)
    assert out.exists() and rows[0]['key'].startswith('LIG:')

def test_explicit_group_key_is_available(tmp_path):
    p=tmp_path/'x.pdb'; make_pdb(p)
    assert reference_candidates(p)[0]['key']=='LIG:A:106'

def test_complete_residue_classification(tmp_path):
    p=tmp_path/'x.pdb'; make_pdb(p); ref=HeteroGroup('LIG','A','106',(ligand(6,'LIG','A',106,2),))
    assert receptor_completeness(p,ref)[0]['classification']=='COMPLETE'

def test_remote_incomplete_residue_classification(tmp_path):
    p=tmp_path/'x.pdb'; p.write_text(atom(1,'N','LYS','A',1,20,0,0)+'\n'); ref=HeteroGroup('LIG','A','100',(ligand(2,'LIG','A',100,0),))
    assert receptor_completeness(p,ref)[0]['classification']=='INCOMPLETE_REMOTE'

def test_binding_site_incomplete_residue_classification(tmp_path):
    p=tmp_path/'x.pdb'; p.write_text(atom(1,'N','LYS','A',1,2,0,0)+'\n'); ref=HeteroGroup('LIG','A','100',(ligand(2,'LIG','A',100,0),))
    assert receptor_completeness(p,ref)[0]['classification']=='INCOMPLETE_BINDING_SITE'

def test_peripheral_incomplete_residue_classification(tmp_path):
    p=tmp_path/'x.pdb'; p.write_text(atom(1,'N','LYS','A',1,6,0,0)+'\n'); ref=HeteroGroup('LIG','A','100',(ligand(2,'LIG','A',100,0),))
    assert receptor_completeness(p,ref)[0]['classification']=='INCOMPLETE_PERIPHERAL'

def test_unknown_residue_is_not_silently_complete(tmp_path):
    p=tmp_path/'x.pdb'; p.write_text(atom(1,'N','UNK','A',1,20,0,0)+'\n'); ref=HeteroGroup('LIG','A','100',(ligand(2,'LIG','A',100,0),))
    assert receptor_completeness(p,ref)[0]['classification']=='UNKNOWN_RESIDUE'

def test_state_stale_can_be_represented(tmp_path):
    update_state(tmp_path,'DOCKING','STALE',errors=['hash mismatch'])
    assert project_status(tmp_path)['stages']['DOCKING']['status']=='STALE'

def test_compound_state_can_be_recorded(tmp_path):
    update_state(tmp_path,'DOCKING','PASS',outputs={'quercetin':'abc'})
    assert project_status(tmp_path)['stages']['DOCKING']['output_hashes']['quercetin']=='abc'

def test_compliance_missing_project_fails(tmp_path):
    result=evaluate(tmp_path)
    assert result['status']=='FAIL'

def test_compliance_requires_raw_hashes(tmp_path):
    (tmp_path/'analysis').mkdir(); (tmp_path/'analysis/evidence.json').write_text(json.dumps({'validation':{'sampling_result':'TARGET_VALIDATION_PASS','native_redocking_RMSD_A':1.0},'interaction_engine':'plip','sidock_internal_parser_status':'EXPERIMENTAL_SECONDARY','failures':[]}))
    assert any(c['rule_id']=='RAW_DLG_HASHES' and c['status']=='FAIL' for c in evaluate(tmp_path)['checks'])

def test_compliance_accepts_valid_core(tmp_path):
    (tmp_path/'analysis').mkdir(); (tmp_path/'provenance').mkdir()
    (tmp_path/'analysis/evidence.json').write_text(json.dumps({'validation':{'sampling_result':'TARGET_VALIDATION_PASS','native_redocking_RMSD_A':1.0},'interaction_engine':'plip','sidock_internal_parser_status':'EXPERIMENTAL_SECONDARY','failures':[]}))
    for f in ['input_immutability.json','raw_dlg_hashes.json','resolved_config.yaml']:
        (tmp_path/'provenance'/f).write_text('{}')
    assert evaluate(tmp_path)['status']=='PASS'

def test_candidate_ranking_is_deterministic(tmp_path):
    p=tmp_path/'x.pdb'; make_pdb(p,('LIG','LIG'),(2,3))
    a=[x.key for x in rank_reference_candidates(p)]; b=[x.key for x in rank_reference_candidates(p)]
    assert a==b
