import json
from pathlib import Path
from types import SimpleNamespace

from sidock_ubuntu.finalization import load_finalization_input, snapshot_evidence
from sidock_ubuntu.plip_runner import materialize_plip_complex
from sidock_ubuntu.stage_outputs import ReferenceStageOutput, save_reference_output
from sidock_ubuntu.state import ensure_compounds, update_compound


def _complete_compound(project, cid='A'):
    ensure_compounds(project,[cid])
    for stage in ('preparation','docking','pose_extraction','clustering','plip','figures'):
        update_compound(project,cid,stage,'PASS')


def test_finalization_recovers_clusters_directly_from_each_dlg(tmp_path):
    project=tmp_path/'p'; project.mkdir(); (project/'docking'/'dlg').mkdir(parents=True)
    _complete_compound(project,'A')
    (project/'docking'/'dlg'/'A.dlg').write_text(
        'Number of conformations = 10\n'
        'CLUSTERING HISTOGRAM\n'
        '  1 | -8.10 | 2 | -7.90 | 6 |\n'
        '  2 | -7.20 | 7 | -7.00 | 4 |\n'
        'Number of multi-member\n'
    )
    evidence=snapshot_evidence(load_finalization_input(project))
    row=evidence['compounds'][0]
    assert row['cluster_count']==2
    assert row['largest_cluster_population']==6


def test_finalization_records_target_and_reference(tmp_path):
    project=tmp_path/'p'; project.mkdir(); (project/'validation').mkdir(); (project/'input').mkdir()
    receptor=project/'input'/'protein.pdb'; receptor.write_text('ATOM\n')
    (project/'project.json').write_text(json.dumps({'receptor_source':str(receptor)}))
    (project/'gui_inputs.json').write_text(json.dumps({'pdb_id':'3A4A'}))
    save_reference_output(project,ReferenceStageOutput('GLC:A:801','GLC','A','801','receptor/reference_ligand.pdb'))
    evidence=snapshot_evidence(load_finalization_input(project))
    assert evidence['target']['pdb']=='3A4A'
    assert evidence['target']['reference_ligand']=='GLC:A:801'


def test_pdbqt_complex_materialization_marks_docked_ligand_as_hetatm(tmp_path):
    source=tmp_path/'complex.pdbqt'
    source.write_text(
        'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00      0.000 C\n'
        'MODEL 1\n'
        'ATOM      1  C1  UNL     1       2.000   0.000   0.000  1.00  0.00      0.000 C\n'
        'ATOM      2  O1  UNL     1       3.000   0.000   0.000  1.00  0.00     -0.200 OA\n'
        'ENDMDL\n'
    )

    def fake_runner(cmd,**kwargs):
        src=Path(cmd[1]); dst=Path(cmd[3])
        if src.name=='receptor.pdbqt':
            dst.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n')
        else:
            dst.write_text(
                'ATOM      1  C1  UNL     1       2.000   0.000   0.000  1.00 20.00           C\n'
                'ATOM      2  O1  UNL     1       3.000   0.000   0.000  1.00 20.00           O\n'
                'CONECT    1    2\nEND\n'
            )
        return SimpleNamespace(returncode=0)

    out=materialize_plip_complex(source,tmp_path/'plip_input.pdb',runner=fake_runner)
    text=out.read_text()
    assert 'HETATM' in text
    assert 'LIG Z 999' in text
    assert 'CONECT' in text
