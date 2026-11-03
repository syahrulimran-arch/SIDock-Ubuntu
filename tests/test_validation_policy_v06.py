import json

from sidock_ubuntu.config import DEFAULT_CONFIG, validate_config
from sidock_ubuntu.finalization import load_finalization_input, snapshot_evidence
from sidock_ubuntu.orchestrator import StageResult, StageRunner
from sidock_ubuntu.stage_outputs import ValidationStageOutput, load_validation_output, save_validation_output
from sidock_ubuntu.state import COMPOUND_STAGES, ensure_compounds, sha256, update_compound, update_state


def _failed_validation(project):
    value = ValidationStageOutput(
        native_dlg=None,
        minimum_rmsd=2.5,
        median_rmsd=3.0,
        lowest_energy_rmsd=2.7,
        rmsd_threshold=2.0,
        scientific_status='FAIL',
        execution_status='SUCCESS',
        native_like_pose_count=0,
    )
    meta = save_validation_output(project, value)
    update_state(
        project,
        'VALIDATION',
        'PASS',
        outputs={str(meta): sha256(meta)},
        result_metadata=str(meta.relative_to(project)),
    )


def test_skip_policy_does_not_execute_native_validation(tmp_path):
    update_state(tmp_path, 'GRID', 'PASS')
    called = []

    def should_not_run():
        called.append(True)
        raise AssertionError('native validation executor must not run in skip mode')

    config = {'validation': {'mode': 'skip', 'rmsd_threshold': 2.0}}
    runner = StageRunner(tmp_path, config, {'validation': should_not_run})
    result = runner.run_global(('VALIDATION',))

    assert result == [('VALIDATION', 'PASS')]
    assert called == []
    validation = load_validation_output(tmp_path)
    assert validation.scientific_status == 'UNVALIDATED'
    assert validation.execution_status == 'SKIPPED'
    policy = json.loads((tmp_path/'provenance'/'validation_policy.json').read_text())
    assert policy['mode'] == 'skip'


def _compound_project(project):
    ensure_compounds(project, ['A'])
    _failed_validation(project)
    return project


def test_required_policy_blocks_compounds_after_failed_native_validation(tmp_path):
    project = _compound_project(tmp_path)
    runner = StageRunner(
        project,
        {'validation': {'mode': 'required'}},
        {'preparation': lambda cid: StageResult('PASS', {})},
    )
    assert runner._run_compound_unlocked('A') == [('preparation', 'BLOCKED')]


def test_if_available_policy_allows_unvalidated_compound_docking(tmp_path):
    project = _compound_project(tmp_path)
    runner = StageRunner(
        project,
        {'validation': {'mode': 'if_available'}},
        {'preparation': lambda cid: StageResult('PASS', {})},
    )
    result = runner._run_compound_unlocked('A')
    assert result[0] == ('preparation', 'PASS')
    assert result[1] == ('docking', 'BLOCKED')  # no docking executor in this focused test


def test_finalization_labels_optional_failed_validation_as_unvalidated(tmp_path):
    ensure_compounds(tmp_path, ['A'])
    for stage in COMPOUND_STAGES:
        update_compound(tmp_path, 'A', stage, 'PASS')
    _failed_validation(tmp_path)
    state = json.loads((tmp_path/'project_state.json').read_text())
    state['stages']['VALIDATION']['warnings'] = [
        'UNVALIDATED_DOCKING',
        'NATIVE_VALIDATION_FAILED_CONTINUING_BY_POLICY',
    ]
    (tmp_path/'project_state.json').write_text(json.dumps(state))
    (tmp_path/'provenance').mkdir(exist_ok=True)
    (tmp_path/'provenance'/'validation_policy.json').write_text(json.dumps({'mode':'if_available'}))

    snapshot = load_finalization_input(tmp_path)
    assert snapshot.status == 'COMPLETE_WITH_WARNINGS'
    evidence = snapshot_evidence(snapshot)
    assert evidence['validation']['sampling_result'] == 'FAIL'
    assert evidence['validation']['docking_validation_status'] == 'UNVALIDATED'
    assert evidence['validation']['validation_policy'] == 'if_available'
    assert evidence['blocked_by'] is None


def test_validation_mode_configuration_is_explicit():
    assert DEFAULT_CONFIG['validation']['mode'] == 'required'
    config = json.loads(json.dumps(DEFAULT_CONFIG))
    config['validation']['mode'] = 'nonsense'
    assert any('validation.mode' in error for error in validate_config(config))
