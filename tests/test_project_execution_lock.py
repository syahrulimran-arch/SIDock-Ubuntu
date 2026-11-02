import json

from sidock_ubuntu.project_lock import lock_path, project_execution_lock


def test_project_lock_is_neutral_and_stale_metadata_is_not_authoritative(tmp_path):
    project = tmp_path / 'project'
    (project / 'provenance').mkdir(parents=True)
    path = lock_path(project)
    path.write_text(json.dumps({'pid': 999999, 'interface': 'CLI'}))
    with project_execution_lock(project, {'interface': 'GUI', 'action': 'RUN'}):
        metadata = json.loads(path.read_text())
        assert path == project / 'provenance' / 'project_execution.lock'
        assert metadata['interface'] == 'GUI'


def test_project_locks_are_independent_per_project(tmp_path):
    first = tmp_path / 'first'
    second = tmp_path / 'second'
    with project_execution_lock(first, {'interface': 'CLI'}):
        with project_execution_lock(second, {'interface': 'GUI'}):
            assert lock_path(first) != lock_path(second)

