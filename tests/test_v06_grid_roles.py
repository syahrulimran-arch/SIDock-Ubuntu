import json
from pathlib import Path

import pytest

from sidock_ubuntu.model import GridBox
from sidock_ubuntu.stages import run_grid_stage


def _receptor(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        'ATOM      1  C   LIG A   1       1.000   2.000   3.000  1.00  0.00      0.000 C\n',
        encoding='utf-8',
    )
    return path


@pytest.mark.parametrize('role', ['production', 'validation'])
def test_grid_stage_persists_role_geometry_maps_and_receptor(tmp_path, monkeypatch, role):
    receptor = _receptor(tmp_path/'prepared'/'source_receptor.pdbqt')
    grid_dir = tmp_path/'grid'/role
    grid = GridBox(11.1, 22.2, 33.3, 60, 62, 64, 0.375)

    def fake_run_command(cmd, cwd, log_path):
        assert Path(cwd) == grid_dir
        assert cmd[1:] == ['-p', 'receptor.gpf', '-l', 'receptor.glg']
        log_path.write_text('autogrid test log\n', encoding='utf-8')
        (grid_dir/'receptor.glg').write_text('GLG\n', encoding='utf-8')
        for name in ('receptor.C.map', 'receptor.e.map', 'receptor.d.map', 'receptor.maps.fld'):
            (grid_dir/name).write_text(name+'\n', encoding='utf-8')
        return 0

    monkeypatch.setattr('sidock_ubuntu.stages.run_command', fake_run_command)
    result = run_grid_stage(receptor, grid_dir, grid, ['C'], '/fake/autogrid4')

    assert result['grid_role'] == role
    assert result['center'] == [11.1, 22.2, 33.3]
    assert result['dimensions'] == [60, 62, 64]
    assert result['spacing'] == 0.375
    assert (grid_dir/'receptor.pdbqt').exists()
    assert {p.name for p in result['maps']} >= {'receptor.C.map', 'receptor.e.map', 'receptor.d.map'}

    provenance = json.loads((grid_dir/'grid_provenance.json').read_text(encoding='utf-8'))
    assert provenance['role'] == role
    assert provenance['center'] == result['center']
    assert provenance['dimensions'] == result['dimensions']
    assert provenance['spacing'] == result['spacing']
    assert provenance['ligand_types'] == ['C']
    assert Path(provenance['receptor_pdbqt']).name == 'receptor.pdbqt'
    assert len(provenance['receptor_sha256']) == 64


def test_grid_stage_legacy_role_remains_supported(tmp_path, monkeypatch):
    receptor = _receptor(tmp_path/'source_receptor.pdbqt')
    grid_dir = tmp_path/'grid'
    grid = GridBox(0.0, 0.0, 0.0, 40, 40, 40, 0.375)

    def fake_run_command(cmd, cwd, log_path):
        log_path.write_text('ok\n', encoding='utf-8')
        (grid_dir/'receptor.glg').write_text('GLG\n', encoding='utf-8')
        (grid_dir/'receptor.C.map').write_text('MAP\n', encoding='utf-8')
        return 0

    monkeypatch.setattr('sidock_ubuntu.stages.run_command', fake_run_command)
    result = run_grid_stage(receptor, grid_dir, grid, ['C'], '/fake/autogrid4')
    assert result['grid_role'] == 'legacy'
    assert json.loads((grid_dir/'grid_provenance.json').read_text())['role'] == 'legacy'
