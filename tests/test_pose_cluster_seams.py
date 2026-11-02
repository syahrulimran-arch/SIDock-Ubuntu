from pathlib import Path
import sidock_ubuntu.analysis as analysis
from sidock_ubuntu.analysis import extract_docking_poses, analyze_autodock_clusters, rank_poses_by_energy

def dlg_fixture(tmp_path):
    d=tmp_path/'x.dlg'; d.write_text('Number of conformations = 2\nDOCKED: MODEL 1\nDOCKED: ATOM      1  C   LIG A   1       0.000   0.000   0.000\nEstimated Free Energy of Binding    = -5.0\nDOCKED: ENDMDL\nDOCKED: MODEL 2\nDOCKED: ATOM      1  C   LIG A   1       1.000   0.000   0.000\nEstimated Free Energy of Binding    = -7.0\nDOCKED: ENDMDL\nCLUSTERING HISTOGRAM\n  1 | -7.0 | 2 | -6.0 | 1 |\nNumber of multi-member\n'); r=tmp_path/'receptor'; r.mkdir(); (r/'receptor.pdbqt').write_text('ATOM\n'); return d

def test_pose_extraction_does_not_call_legacy(tmp_path,monkeypatch):
    d=dlg_fixture(tmp_path); monkeypatch.setattr(analysis,'analyze_docking_results',lambda *a: (_ for _ in ()).throw(RuntimeError('LEGACY_ANALYSIS_CALLED')))
    out=extract_docking_poses(d,tmp_path/'poses',tmp_path/'complexes',tmp_path/'receptor/receptor.pdbqt')
    assert out.selected_model==2 and out.lowest_energy==-7.0

def test_cluster_parser_does_not_call_legacy(tmp_path,monkeypatch):
    d=dlg_fixture(tmp_path); monkeypatch.setattr(analysis,'analyze_docking_results',lambda *a: (_ for _ in ()).throw(RuntimeError('LEGACY_ANALYSIS_CALLED')))
    out=analyze_autodock_clusters(d,tmp_path/'clusters.csv'); assert out.cluster_count==1 and out.largest_cluster_population==1

def test_energy_ranking_preserves_lowest_first():
    rows=rank_poses_by_energy([{'model':1,'energy':-5},{'model':2,'energy':-7}]); assert [x['model'] for x in rows]==[2,1]

def test_pose_and_cluster_do_not_need_contacts_or_figures(tmp_path,monkeypatch):
    d=dlg_fixture(tmp_path); monkeypatch.setattr(analysis,'detect_contacts',lambda *a: (_ for _ in ()).throw(RuntimeError('CONTACTS_CALLED')))
    monkeypatch.setattr(analysis,'write_contacts',lambda *a: (_ for _ in ()).throw(RuntimeError('CONTACTS_CALLED')))
    extract_docking_poses(d,tmp_path/'poses',tmp_path/'complexes',tmp_path/'receptor/receptor.pdbqt'); analyze_autodock_clusters(d)
