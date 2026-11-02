import json
from pathlib import Path

from sidock_ubuntu.finalization import load_finalization_input, snapshot_evidence
from sidock_ubuntu.state import aggregate_compound_stage_status, derived_project_status, ensure_compounds, update_compound, update_state


def _state_project(tmp_path, compounds):
    p = tmp_path / 'p'; p.mkdir()
    (p / 'project_state.json').write_text(json.dumps({'stages': {x: {'status': 'PENDING'} for x in ('REFERENCE','RECEPTOR','GRID','VALIDATION','FINALIZATION')}, 'compounds': {}}))
    ensure_compounds(p, list(compounds))
    return p


def test_all_compound_stages_aggregate_without_global_attempts(tmp_path):
    p = _state_project(tmp_path, ['A', 'B'])
    update_state(p, 'VALIDATION', 'PASS')
    for stage in ('REFERENCE','RECEPTOR','GRID'):
        update_state(p, stage, 'PASS')
    for cid in ('A', 'B'):
        for stage in ('preparation','docking','pose_extraction','clustering','plip','figures'):
            update_compound(p, cid, stage, 'PASS')
    update_state(p, 'FINALIZATION', 'PASS')
    state = json.loads((p/'project_state.json').read_text())
    assert aggregate_compound_stage_status(state, 'docking') == 'PASS'
    assert derived_project_status(state) == 'COMPLETE'


def test_mixed_compounds_are_partial_and_validation_blocks_downstream(tmp_path):
    p = _state_project(tmp_path, ['A', 'B'])
    update_state(p, 'VALIDATION', 'PASS')
    update_compound(p, 'A', 'docking', 'PASS'); update_compound(p, 'B', 'docking', 'FAIL')
    state = json.loads((p/'project_state.json').read_text())
    assert aggregate_compound_stage_status(state, 'docking') == 'PARTIAL'
    update_state(p, 'VALIDATION', 'FAIL')
    state = json.loads((p/'project_state.json').read_text())
    assert aggregate_compound_stage_status(state, 'docking') == 'BLOCKED'


def test_disk_finalization_preserves_zero_plip_counts_and_metrics(tmp_path):
    p = _state_project(tmp_path, ['A'])
    update_state(p, 'VALIDATION', 'PASS')
    for stage in ('preparation','docking','pose_extraction','clustering','plip','figures'):
        update_compound(p, 'A', stage, 'PASS')
    (p/'analysis').mkdir()
    (p/'analysis'/'ranking.csv').write_text('rank,ligand,energy_kcal_mol,pose_file,complex_file\n1,A,-7.25,analysis/poses/A_model_3.pdbqt,analysis/complexes/A_model_3_complex.pdbqt\n')
    (p/'analysis'/'clusters.csv').write_text('ligand,cluster_rank,cluster_size\nA,1,4\nA,2,2\n')
    (p/'analysis'/'interactions_plip_A.csv').write_text('compound_id,pose_id,interaction_type\nA,3,hydrogen_bond\nA,3,hydrophobic\n')
    e = snapshot_evidence(load_finalization_input(p))
    c = e['compounds'][0]
    assert (c['lowest_docking_energy'], c['cluster_count'], c['largest_cluster_population']) == (-7.25, 2, 4)
    assert (c['PLIP_H_bonds'], c['PLIP_hydrophobic'], c['PLIP_pi_stack'], c['PLIP_salt_bridge']) == (1, 1, 0, 0)
