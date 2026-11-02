from pathlib import Path
import inspect
import sidock_ubuntu.analysis as analysis

def _dlg(tmp_path):
    d=tmp_path/'docking/dlg'; d.mkdir(parents=True); (tmp_path/'receptor').mkdir(); (tmp_path/'receptor/receptor.pdbqt').write_text('ATOM\n')
    p=d/'x.dlg'; p.write_text('Number of conformations = 1\nDOCKED: MODEL 1\nDOCKED: ATOM      1  C   LIG A   1       0.000   0.000   0.000\nEstimated Free Energy of Binding    = -1.0\nDOCKED: ENDMDL\nCLUSTERING HISTOGRAM\n  1 | -1.0 | 1 | -1.0 | 1 |\nNumber of multi-member\n')
    return d

def test_legacy_api_delegates_pose_and_cluster(tmp_path,monkeypatch):
    d=_dlg(tmp_path); calls=[]
    monkeypatch.setattr(analysis, 'detect_contacts', lambda *a, **k: [])
    original_pose=analysis.extract_docking_poses; original_cluster=analysis.analyze_autodock_clusters
    monkeypatch.setattr(analysis,'extract_docking_poses',lambda *a,**k: calls.append('pose') or original_pose(*a,**k))
    monkeypatch.setattr(analysis,'analyze_autodock_clusters',lambda *a,**k: calls.append('cluster') or original_cluster(*a,**k))
    analysis.analyze_docking_results(d,tmp_path/'analysis')
    assert calls==['pose','cluster']

def test_legacy_api_keeps_secondary_contact_output(tmp_path, monkeypatch):
    monkeypatch.setattr(analysis, 'detect_contacts', lambda *a, **k: [])
    out=analysis.analyze_docking_results(_dlg(tmp_path),tmp_path/'analysis')
    assert str(out['contacts']).endswith('interaction_crosscheck/interactions_sidock_experimental.csv')

def test_legacy_api_contains_no_duplicate_pose_or_cluster_parser():
    source = inspect.getsource(analysis.analyze_docking_results)
    assert '_write_pose_and_complexes' not in source
    assert 'parse_dlg_clusters(' not in source
    assert 'extract_docking_poses(' in source
    assert 'analyze_autodock_clusters(' in source
