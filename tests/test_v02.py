from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from sidock_ubuntu.analysis import parse_dlg_clusters
from sidock_ubuntu.contacts import detect_contacts
from sidock_ubuntu.pipeline import enforce_redocking_validation
from sidock_ubuntu.validation import heavy_atom_rmsd, fixed_frame_symmetry_corrected_heavy_atom_rmsd
from sidock_ubuntu.core import write_native_start_dpf
from sidock_ubuntu.export import export_ligand_complex
from sidock_ubuntu.model import DockingParameters


def _mol(distance: float):
    mol = Chem.AddHs(Chem.MolFromSmiles('CC'))
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.SetAtomPosition(0, (0, 0, 0)); conf.SetAtomPosition(1, (distance, 0, 0))
    for i in range(2, mol.GetNumAtoms()): conf.SetAtomPosition(i, (0, 0, 0))
    mol.AddConformer(conf)
    return mol


def test_rmsd_mapping_pass_and_fail():
    crystal = _mol(1.5)
    assert heavy_atom_rmsd(crystal, _mol(1.5)) <= 2.0
    assert heavy_atom_rmsd(crystal, _mol(6.0)) > 2.0


def test_fixed_frame_rmsd_translation_and_rotation():
    crystal = _mol(1.5)
    translated = _mol(1.5)
    conf = translated.GetConformer()
    for i in range(translated.GetNumAtoms()):
        p = conf.GetAtomPosition(i); conf.SetAtomPosition(i, (p.x + 10.0, p.y, p.z))
    assert fixed_frame_symmetry_corrected_heavy_atom_rmsd(crystal, translated) == pytest.approx(10.0)
    rotated = _mol(1.5)
    conf = rotated.GetConformer()
    for i in range(rotated.GetNumAtoms()):
        p = conf.GetAtomPosition(i); conf.SetAtomPosition(i, (-p.x, -p.y, p.z))
    assert fixed_frame_symmetry_corrected_heavy_atom_rmsd(crystal, rotated) > 0.0


def test_native_start_dpf_uses_absolute_ligand_center(tmp_path: Path):
    ligand = tmp_path / 'ligand.pdbqt'
    ligand.write_text('ATOM      1  C   UNL     1      10.000  20.000  30.000  1.00  0.00    +0.000 C  \nATOM      2  C   UNL     1      12.000  20.000  30.000  1.00  0.00    +0.000 C  \nTORSDOF 0\n')
    out = tmp_path / 'native.dpf'; write_native_start_dpf('r', ligand, ['C'], DockingParameters(), out)
    text = out.read_text(); assert 'about 11.000 20.000 30.000' in text; assert 'tran0 11.000 20.000 30.000' in text; assert 'quaternion0 0.0 0.0 0.0 1.0' in text; assert 'tran0 0 0 0' not in text


def test_epdb_dpf_generation(tmp_path: Path):
    ligand = tmp_path / 'ligand.pdbqt'
    ligand.write_text('ATOM      1  C   UNL     1      10.000  20.000  30.000  1.00  0.00    +0.000 C  \nTORSDOF 0\n')
    out = tmp_path / 'epdb.dpf'; write_native_start_dpf('r', ligand, ['C'], DockingParameters(), out, epdb=True)
    text = out.read_text(); assert '\nepdb\n' in text; assert 'do_local_only' not in text


def test_plip_complex_export_preserves_hetatm_identity(tmp_path: Path):
    receptor = tmp_path / 'r.pdb'; template = tmp_path / 't.pdb'; pose = tmp_path / 'p.pdb'; out = tmp_path / 'x.pdb'
    receptor.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n')
    template.write_text('HETATM    2  C1  1N1 A 501       1.000   1.000   1.000  1.00 20.00           C  \nEND\n')
    pose.write_text('ATOM      9  C   UNL     1       2.000   2.000   2.000  1.00  0.00           C  \n')
    export_ligand_complex(receptor, template, pose, out)
    text = out.read_text(); assert 'HETATM' in text and ' 1N1 A 501' in text


def test_dlg_cluster_parser(tmp_path: Path):
    dlg = tmp_path / 'x.dlg'
    dlg.write_text('Number of conformations = 4\nCLUSTERING HISTOGRAM\n   1 |     -8.20 |   3 |     -7.50 |   3 |###\n   2 |     -6.00 |   1 |     -6.00 |   1 |#\nNumber of multi-member conformational clusters found = 1\n')
    rows = parse_dlg_clusters(dlg)
    assert rows[0]['cluster_size'] == 3
    assert rows[0]['run_number'] == 3
    assert rows[0]['cluster_percentage'] == 75.0


def test_contact_parser_known_geometry(tmp_path: Path):
    receptor = tmp_path / 'r.pdbqt'; ligand = tmp_path / 'l.pdbqt'; out = tmp_path / 'contacts.csv'
    receptor.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00    +0.000 C \n')
    ligand.write_text('ATOM      1  C   UNL     1       3.000   0.000   0.000  1.00  0.00    +0.000 A \n')
    rows = detect_contacts(receptor, ligand, 'L1', out)
    assert rows and rows[0]['interaction_type'] == 'hydrophobic'


def test_failed_validation_blocks_pipeline():
    with pytest.raises(RuntimeError):
        enforce_redocking_validation({'status': 'FAIL'})
    enforce_redocking_validation({'status': 'FAIL'}, allow_failed=True)
