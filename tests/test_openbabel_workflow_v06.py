"""Preparation and GUI regressions run with no Meeko installation required."""
import json
import shutil
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from sidock_ubuntu import core, openbabel_receptor as ob
from sidock_ubuntu.guided_workflow import create_guided_project, select_protein
from sidock_ubuntu.manual_docking import root, execute


def molecule(smiles="CCO"):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=1234) == 0
    return mol


def write_sdf(path):
    with Chem.SDWriter(str(path)) as writer:
        writer.write(molecule())


@pytest.mark.skipif(not shutil.which("obabel"), reason="Open Babel executable required")
@pytest.mark.parametrize("smiles", ["CCO", "CC(=O)Nc1ccccc1"])
def test_real_openbabel_ligand_preserves_heavy_geometry(tmp_path, smiles):
    source = tmp_path / "ligand.sdf"
    with Chem.SDWriter(str(source)) as writer:
        writer.write(molecule(smiles))
    target = tmp_path / "first-run/prepared/ligand.pdbqt"
    before = source.read_bytes()
    assert core.prepare_ligand(source, target, tmp_path / "logs/preparation.log") == target
    assert source.read_bytes() == before
    assert "TORSDOF" in target.read_text()


def test_empty_ligand_output_rejected(tmp_path, monkeypatch):
    source = tmp_path / "ligand.sdf"
    write_sdf(source)
    target = tmp_path / "new/ligand.pdbqt"
    monkeypatch.setattr(core, "resolve_executable", lambda *_: "/fake/obabel")

    def fake_run(command, cwd, log):
        assert cwd.is_dir()
        target.write_text("ROOT\nENDROOT\nTORSDOF 0\n")
        return 0

    monkeypatch.setattr(core, "run_command", fake_run)
    with pytest.raises(RuntimeError, match="no atoms"):
        core.prepare_ligand(source, target, tmp_path / "logs/prepare.log")


def test_prepare_chain_calls_all_checks_without_docking_or_approval(tmp_path, monkeypatch):
    project = create_guided_project(tmp_path, "review")
    home = root(project)
    (home / "custom_parameters.dat").write_text("fixture-only\n")
    calls = []
    for name in ("generate_candidate", "generate_protonated_charge_candidate",
                 "audit_and_correct_hydrogen_types", "prepare_exploratory_receptor"):
        def stage(p, name=name):
            assert p == project
            calls.append(name)
            return {"candidate": "fixture.pdbqt", "scientific_approval": False}
        monkeypatch.setattr(ob, name, stage)
    result = execute(project, "receptor_openbabel_pipeline")
    assert calls == ["generate_candidate", "generate_protonated_charge_candidate",
                     "audit_and_correct_hydrogen_types", "prepare_exploratory_receptor"]
    assert result["scientific_approval"] is False
    assert not (home / "receptor.pdbqt").exists()
    assert not (home / "exploratory_smoke").exists()
    assert not (home / "runs").exists()
    assert not json.loads((home / "workflow.json").read_text())["approvals"]


def project_with_protein(tmp_path):
    project = create_guided_project(tmp_path, "protein")
    pdb = project / "input/protein.pdb"
    # Known chemistry with finite 3D coordinates plus an isolated Ni site.
    mol = Chem.RemoveHs(molecule("NCC(=O)O"))
    block = Chem.MolToPDBBlock(mol)
    metal = ("HETATM  100   NI  NI Z 999      20.000  20.000  20.000"
             "  1.00 20.00          NI\n")
    pdb.write_text(block.replace("END\n", "") + metal + "END\n")
    select_protein(project, ligand_key=None, remove_groups=[], center=(0, 0, 0))
    return project


@pytest.mark.skipif(not shutil.which("obabel"), reason="Open Babel executable required")
def test_real_openbabel_receptor_chain_retains_metal_and_requires_review(tmp_path):
    project = project_with_protein(tmp_path)
    home = root(project)
    before = (home / "inputs/protein_clean.pdb").read_bytes()
    # These synthetic values test type coverage only; no docking is run.
    types = ("C", "A", "N", "NA", "OA", "O", "HD", "H", "Ni")
    (home / "custom_parameters.dat").write_text("".join(
        f"atom_par {typ} 1 1 1 1 1 1 1 1 1 1\n" for typ in types))
    result = ob.prepare_openbabel_receptor(project)
    assert result["scientific_approval"] is False
    assert (home / "inputs/protein_clean.pdb").read_bytes() == before
    assert not (home / "receptor.pdbqt").exists()
    audit = result["stages"]["exploratory"]
    assert audit["heavy_atoms_preserved_exactly"] is True
    assert audit["ni_count"] == 1
    assert tuple(audit["ni_coordinates"][0]) == (20, 20, 20)
    assert audit["ni_charges_e"] == [0.0]
    assert all(stage["status"] == "REUSED" for stage in
               ob.prepare_openbabel_receptor(project)["stages"].values())


def test_gui_can_launch_preparation_and_separate_technical_check(tmp_path, monkeypatch):
    pytest.importorskip("httpx")
    from fastapi.testclient import TestClient
    from sidock_ubuntu.gui import create_app
    from sidock_ubuntu import manual_docking
    project = project_with_protein(tmp_path)
    calls = []
    monkeypatch.setattr(manual_docking, "launch", lambda p, action: calls.append(action))
    with TestClient(create_app(tmp_path)) as client:
        for action in ("receptor_openbabel_pipeline", "exploratory_smoke"):
            response = client.post(f"/projects/{project.name}/guided/action",
                                   data={"action": action}, follow_redirects=False)
            assert response.status_code == 303
            assert response.headers["location"].endswith("/guided/2")
    assert calls == ["receptor_openbabel_pipeline", "exploratory_smoke"]


def test_packaging_does_not_require_meeko():
    config = (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text()
    assert "meeko" not in config.lower()
