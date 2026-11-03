import os
from pathlib import Path


def test_figures_module_forces_headless_backend(monkeypatch, tmp_path):
    # Importing the module must configure Matplotlib for a file-only backend
    # before pyplot can initialize Tk/Tcl inside a GUI worker process.
    monkeypatch.delenv('MPLBACKEND', raising=False)
    import importlib
    import sidock_ubuntu.figures as figures
    importlib.reload(figures)
    assert os.environ.get('MPLBACKEND') == 'Agg'


def test_activity_scatter_uses_agg_and_writes_png(tmp_path):
    from sidock_ubuntu.figures import activity_scatter
    import matplotlib

    out = tmp_path / 'scatter.png'
    result = activity_scatter(
        {'rows': [
            {'ligand': 'A', 'energy_kcal_mol': -7.0, 'pIC50': 6.0},
            {'ligand': 'B', 'energy_kcal_mol': -8.0, 'pIC50': 7.0},
        ]},
        out,
    )
    assert result == out
    assert out.is_file() and out.stat().st_size > 0
    assert matplotlib.get_backend().lower() == 'agg'
