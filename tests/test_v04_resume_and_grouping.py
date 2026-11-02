import json
from sidock_ubuntu.state import (COMPOUND_STAGES, config_hash, ensure_compounds, invalidate_for_config_change,
                                  invalidate_from, project_status, update_compound, verify_output_hashes)
from sidock_ubuntu.reference import confidence_for_candidates, group_equivalent_instances
from sidock_ubuntu.config import DEFAULT_CONFIG, deep_merge

def test_compound_template_has_all_stages(tmp_path):
    ensure_compounds(tmp_path,['q']); assert set(project_status(tmp_path)['compounds']['q'])==set(COMPOUND_STAGES)

def test_compound_plip_resume_state(tmp_path):
    ensure_compounds(tmp_path,['q']); update_compound(tmp_path,'q','preparation','PASS'); update_compound(tmp_path,'q','docking','PASS'); update_compound(tmp_path,'q','plip','FAIL',error_code='E701')
    assert project_status(tmp_path)['compounds']['q']['plip']['status']=='FAIL'

def test_compound_attempt_history_retained(tmp_path):
    ensure_compounds(tmp_path,['q']); update_compound(tmp_path,'q','docking','FAIL',error_code='E502'); update_compound(tmp_path,'q','docking','PASS')
    assert len(project_status(tmp_path)['compounds']['q']['docking']['attempts'])==2

def test_compound_output_hashes_recorded(tmp_path):
    ensure_compounds(tmp_path,['q']); update_compound(tmp_path,'q','docking','PASS',outputs={'x':'abc'})
    assert project_status(tmp_path)['compounds']['q']['docking']['output_hashes']['x']=='abc'

def test_stage_invalidation_propagates():
    state={'stages':{x:{'status':'PASS','warnings':[]} for x in ['DOCKING','POSE_EXTRACTION','CLUSTERING','PLIP','FIGURES','FINALIZATION']}}
    invalidate_from(state,'DOCKING','OUTPUT_HASH_MISMATCH')
    assert all(v['status']=='STALE' for v in state['stages'].values())

def test_config_interaction_change_does_not_invalidate_docking():
    old=DEFAULT_CONFIG; new=deep_merge(old,{'interaction':{'primary_engine':'plip'}})
    assert config_hash(old,'DOCKING')==config_hash(new,'DOCKING')

def test_config_ga_change_invalidates_docking_hash():
    old=DEFAULT_CONFIG; new=deep_merge(old,{'docking':{'evaluations':999}})
    assert config_hash(old,'DOCKING')!=config_hash(new,'DOCKING')

def test_config_threshold_change_invalidates_validation_hash():
    old=DEFAULT_CONFIG; new=deep_merge(old,{'validation':{'rmsd_threshold':1.5}})
    assert config_hash(old,'VALIDATION')!=config_hash(new,'VALIDATION')

def test_config_figure_change_invalidates_figures_hash():
    old=DEFAULT_CONFIG; new=deep_merge(old,{'figures':{'dpi':300}})
    assert config_hash(old,'FIGURES')!=config_hash(new,'FIGURES')

def test_config_unrelated_change_preserves_docking_hash():
    old=DEFAULT_CONFIG; new=deep_merge(old,{'output':{'figures':False}})
    assert config_hash(old,'DOCKING')==config_hash(new,'DOCKING')

def test_equivalent_identity_grouping():
    rows=[{'ligand_code':'P43','heavy_atom_count':10,'formal_charge':0,'chemical_fingerprint':'abc','key':'A'}, {'ligand_code':'P43','heavy_atom_count':10,'formal_charge':0,'chemical_fingerprint':'abc','key':'B'}]
    assert len(group_equivalent_instances(rows))==1

def test_same_size_different_graph_not_grouped():
    rows=[{'ligand_code':'X','heavy_atom_count':10,'formal_charge':0,'chemical_fingerprint':'abc'}, {'ligand_code':'X','heavy_atom_count':10,'formal_charge':0,'chemical_fingerprint':'def'}]
    assert len(group_equivalent_instances(rows))==2

def test_low_confidence_distinct_candidates():
    assert confidence_for_candidates([{'ligand_code':'A'},{'ligand_code':'B'}])=='LOW'

def test_moderate_confidence_equivalent_candidates():
    assert confidence_for_candidates([{'ligand_code':'A'},{'ligand_code':'A'}])=='MODERATE'

def test_high_confidence_single_candidate():
    assert confidence_for_candidates([{'ligand_code':'A'}])=='HIGH'

def test_hash_mismatch_is_detectable(tmp_path):
    f=tmp_path/'out'; f.write_text('one'); ensure_compounds(tmp_path,['q']); update_compound(tmp_path,'q','docking','PASS',outputs={str(f):'bad'})
    assert verify_output_hashes(tmp_path,project_status(tmp_path)['compounds']['q']['docking']) is False

def test_hash_match_is_detectable(tmp_path):
    f=tmp_path/'out'; f.write_text('one'); from sidock_ubuntu.state import sha256
    ensure_compounds(tmp_path,['q']); update_compound(tmp_path,'q','docking','PASS',outputs={str(f):sha256(f)})
    assert verify_output_hashes(tmp_path,project_status(tmp_path)['compounds']['q']['docking']) is True
