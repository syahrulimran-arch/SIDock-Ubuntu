import os


def test_sidock_figures_force_headless_matplotlib_backend(tmp_path):
    from sidock_ubuntu.figures import activity_scatter

    assert os.environ.get('MPLBACKEND') == 'Agg'
    out=activity_scatter({'rows':[
        {'energy_kcal_mol':-7.0,'pIC50':5.0,'ligand':'A'},
        {'energy_kcal_mol':-8.0,'pIC50':6.0,'ligand':'B'},
    ]}, tmp_path/'scatter.png')
    assert out is not None and out.exists()
    import matplotlib
    assert str(matplotlib.get_backend()).lower() == 'agg'
