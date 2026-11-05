"""Step 3 regression: ligands can be prepared without a prepared receptor.

The Meeko executable is mocked in tests; actual Ubuntu Meeko checks are separate.
"""
import json
import re

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from sidock_ubuntu import manual_docking as manual
from sidock_ubuntu.guided_ui import step3
from sidock_ubuntu.guided_workflow import (
    create_guided_project, compound_profile, approved,
)


def new_compound_only_project(tmp_path, count=16):
    project = create_guided_project(tmp_path, "independent-ligand-step")
    source = project / "input/compounds.sdf"
    writer = Chem.SDWriter(str(source))
    for i in range(1, count + 1):
        mol = Chem.MolFromSmiles("CCO")
        mol.SetProp("_Name", f"AA-{i}")
        mol.SetProp("Compound_ID", f"AA-{i}")
        writer.write(mol)
    writer.close()
    config_file = project / "project.json"
    config = json.loads(config_file.read_text())
    config["ligand_source"] = str(source.resolve())
    config_file.write_text(json.dumps(config))
    return project


def fake_split(source, out_dir, *, reuse_existing=False):
    out_dir.mkdir(parents=True, exist_ok=True)
    output = []
    supplier = Chem.SDMolSupplier(str(source), removeHs=False)
    for mol in supplier:
        name = mol.GetProp("Compound_ID")
        out = out_dir / f"{name}.sdf"
        if not out.exists():
            molecule = Chem.AddHs(mol)
            assert AllChem.EmbedMolecule(molecule, randomSeed=42) == 0
            w = Chem.SDWriter(str(out))
            w.write(molecule)
            w.close()
        output.append(out)
    return output


def test_step3_button_enabled_before_step2_approval(tmp_path):
    project = new_compound_only_project(tmp_path)
    assert not approved(project, 2)
    html = step3(project)
    assert "Molecules:</b> 16" in html
    assert "Prepared:</b> 0" in html
    assert "Preparation is independent of Step 2" in html
    assert re.search(
        r"value='prepare'>\s*<button type='submit'\s*>Minimize and prepare compounds",
        html,
    )
    assert "Step 3 approval unavailable" in html


def test_step3_prepares_all_16_with_no_receptor_and_reuses(tmp_path, monkeypatch):
    project = new_compound_only_project(tmp_path)
    calls = []
    monkeypatch.setattr(manual, "split_and_minimize_sdf", fake_split)

    def fake_meeko(source, destination, log):
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(
            "ROOT\n"
            "HETATM    1  C1  UNL A   1      10.000  20.000  30.000  1.00  0.00    +0.100 C\n"
            "ENDROOT\nTORSDOF 0\n"
        )
        calls.append(destination.name)
        return destination

    monkeypatch.setattr(manual, "prepare_ligand", fake_meeko)
    result = manual.prepare(project)
    assert result["prepared"] == 16
    assert len(calls) == 16
    assert not (manual.root(project) / "receptor.pdbqt").exists()
    profile = compound_profile(project)
    assert profile["prepared"]
    assert profile["prepared_count"] == 16
    assert not approved(project, 2)
    assert not approved(project, 3)
    html = step3(project)
    assert "Step 3 approval is pending" in html
    assert "no rerun needed" in html

    def no_rerun(*args, **kwargs):
        raise AssertionError("Reprocessing a complete library must not happen")

    monkeypatch.setattr(manual, "split_and_minimize_sdf", no_rerun)
    monkeypatch.setattr(manual, "prepare_ligand", no_rerun)
    assert manual.prepare(project)["status"] == "REUSED"


def test_partial_ligand_preparation_resumes_without_overwriting(tmp_path, monkeypatch):
    project = new_compound_only_project(tmp_path)
    monkeypatch.setattr(manual, "split_and_minimize_sdf", fake_split)
    calls = []
    fail = {"enabled": True}

    def occasional_failure(source, destination, log):
        calls.append(destination.name)
        if destination.name == "AA-7.pdbqt" and fail["enabled"]:
            raise RuntimeError("simulated Meeko error for AA-7")
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(
            "ROOT\nHETATM    1  C1  UNL A   1      10.000  20.000  30.000  1.00  0.00    +0.100 C\nENDROOT\nTORSDOF 0\n"
        )
        return destination

    monkeypatch.setattr(manual, "prepare_ligand", occasional_failure)
    with pytest.raises(RuntimeError, match="AA-7"):
        manual.prepare(project)
    completed = manual.root(project) / "ligands/prepared"
    snapshot = {p.name: p.read_bytes() for p in completed.glob("*.pdbqt")}
    assert len(snapshot) == 6
    assert not (manual.root(project) / "ligands/manifest.json").exists()

    fail["enabled"] = False
    result = manual.prepare(project)
    assert result["prepared"] == 16
    assert len(calls) == 17  # first seven attempted; only ten new on retry
    assert {p.name: p.read_bytes() for p in completed.glob("*.pdbqt")
            if p.name in snapshot} == snapshot
    assert compound_profile(project)["prepared"]


def test_incomplete_preparation_cannot_be_marked_verified(tmp_path):
    project = new_compound_only_project(tmp_path)
    folder = manual.root(project) / "ligands/prepared"
    folder.mkdir(parents=True)
    (folder / "AA-1.pdbqt").write_text("ROOT\nENDROOT\n")
    assert compound_profile(project)["prepared_count"] == 0
    assert not compound_profile(project)["prepared"]



def test_meeko_preparation_creates_empty_destination_folder(tmp_path, monkeypatch):
    from sidock_ubuntu import core

    source = tmp_path / "source.sdf"
    source.write_text("fixture")
    destination = tmp_path / "first-run/ligands/prepared/AA-1.pdbqt"
    log_path = tmp_path / "first-run/logs/prepare_AA-1.log"
    monkeypatch.setattr(core, "resolve_meeko_script", lambda *_: "mk_prepare_ligand.py")

    def fake_run(command, cwd, log):
        assert cwd.is_dir()
        assert destination.parent.is_dir()
        destination.write_text(
            "ROOT\nHETATM    1  C1  UNL A   1      10.000  20.000  30.000  1.00  0.00    +0.100 C\nENDROOT\n"
        )
        return 0

    monkeypatch.setattr(core, "run_command", fake_run)
    assert core.prepare_ligand(source, destination, log_path) == destination
    assert destination.stat().st_size > 0


def test_meeko_empty_output_rejected(tmp_path, monkeypatch):
    from sidock_ubuntu import core

    source = tmp_path / "sample.sdf"
    source.write_text("fixture")
    target = tmp_path / "prepared/sample.pdbqt"
    monkeypatch.setattr(core, "resolve_meeko_script", lambda *_: "mk_prepare_ligand.py")

    def invalid_run(command, cwd, log):
        target.write_text("ROOT\nENDROOT\n")
        return 0

    monkeypatch.setattr(core, "run_command", invalid_run)
    with pytest.raises(RuntimeError, match="no atoms"):
        core.prepare_ligand(source, target, tmp_path / "logs/ligand.log")
