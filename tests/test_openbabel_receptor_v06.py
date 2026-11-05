"""Step 2 rigid PDBQT export without Meeko: no AutoDock execution."""
import json
from pathlib import Path

import pytest

from sidock_ubuntu.guided_workflow import create_guided_project, select_protein
from sidock_ubuntu.guided_ui import step2
from sidock_ubuntu.manual_docking import root
from sidock_ubuntu import openbabel_receptor as ob


def pdb_atom(serial, name, residue, chain, number, xyz, element):
    x, y, z = xyz
    return (
        f"HETATM{serial:5d} {name:>4s} {residue:>3s} {chain}{number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2s}\n"
    )


def source_project(tmp_path, *, metal=True):
    project = create_guided_project(tmp_path, "sample")
    pdb = project / "input/protein.pdb"
    pdb.write_text(
        pdb_atom(1, "CA", "ALA", "A", 1, (4.0, 5.0, 6.0), "C")
        + pdb_atom(2, "O", "ALA", "A", 1, (5.0, 6.0, 7.0), "O")
        + (pdb_atom(3, "NI", "NI", "A", 201, (7.2, 8.0, 9.4), "NI")
           if metal else "")
        + "END\n"
    )
    project_json = project / "project.json"
    data = json.loads(project_json.read_text())
    data["receptor_source"] = str(pdb)
    project_json.write_text(json.dumps(data))
    select_protein(
        project, ligand_key=None, remove_groups=[],
        center=(5.0, 6.0, 7.0), remove_waters=False
    )
    return project


def install_fake_obabel(monkeypatch, *, missing_ni=False):
    monkeypatch.setattr(ob.shutil, "which", lambda executable: "/usr/bin/obabel" if executable == "obabel" else None)

    def run(command, cwd, log):
        assert "-xr" in command and "-xc" in command
        assert "-xn" in command and "-xp" in command and "-xh" in command
        assert command[-2:] == ["--partialcharge", "gasteiger"]
        original = Path(command[1]).read_text().splitlines()
        converted = []
        for line in original:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            element = line[76:78].strip()
            if missing_ni and element == "NI":
                continue
            typ = "Ni" if element == "NI" else "OA" if element == "O" else "C"
            converted.append(line[:66].ljust(70) + f"{0.125:6.3f}" + f" {typ:>2s}")
        Path(command[3]).write_text("\n".join(converted) + "\nEND\n")
        Path(log).parent.mkdir(parents=True, exist_ok=True)
        Path(log).write_text("1 molecule converted\n")
        return 0

    from sidock_ubuntu import core
    monkeypatch.setattr(core, "run_command", run)


def test_step2_openbabel_generates_downloadable_rigid_receptor_without_meeko(tmp_path, monkeypatch):
    project = source_project(tmp_path)
    original = (project / "input/protein.pdb").read_bytes()
    reviewed = (root(project) / "inputs/protein_clean.pdb").read_bytes()
    install_fake_obabel(monkeypatch)

    output = ob.generate_candidate(project)
    assert output["status"] == "CANDIDATE_REQUIRES_REVIEW"
    assert output["input_heavy_atom_count"] == 3
    assert output["metal_atoms_missing"] == 0
    assert output["missing_input_heavy_atoms"] == []
    assert output["metal_parameterization"] == "NOT_VALIDATED"
    candidate = Path(output["candidate"])
    assert candidate.is_file()
    assert candidate.suffix == ".pdbqt"
    assert "NI A 201" in candidate.read_text()
    assert not (root(project) / "receptor.pdbqt").exists()
    assert (project / "input/protein.pdb").read_bytes() == original
    assert (root(project) / "inputs/protein_clean.pdb").read_bytes() == reviewed
    assert (root(project) / "inputs/receptor_openbabel_audit.json").is_file()
    assert "receptor_openbabel_candidate.pdbqt" in step2(project)
    assert ob.generate_candidate(project)["status"] == "REUSED"


def test_step2_openbabel_flags_lost_metal_without_discarding_candidate(tmp_path, monkeypatch):
    project = source_project(tmp_path)
    install_fake_obabel(monkeypatch, missing_ni=True)
    result = ob.generate_candidate(project)
    assert result["status"] == "CANDIDATE_MISSING_ATOMS"
    assert result["metal_atoms_missing"] == 1
    assert len(result["missing_input_heavy_atoms"]) == 1
    assert Path(result["candidate"]).is_file()
    assert not (root(project) / "receptor.pdbqt").exists()


def test_step2_openbabel_nonmetal_input_is_supported(tmp_path, monkeypatch):
    project = source_project(tmp_path, metal=False)
    install_fake_obabel(monkeypatch)
    result = ob.generate_candidate(project)
    assert result["input_heavy_atom_count"] == 2
    assert result["metal_atoms_input"] == 0
    assert result["missing_input_heavy_atoms"] == []
    assert "Generate rigid receptor PDBQT with Open Babel" in step2(project)


def test_step2_openbabel_missing_executable_reports_actionable_error(tmp_path, monkeypatch):
    project = source_project(tmp_path)
    monkeypatch.setattr(ob.shutil, "which", lambda _: None)
    with pytest.raises(RuntimeError, match="obabel"):
        ob.generate_candidate(project)
    assert not (root(project) / "inputs/receptor_openbabel_candidate.pdbqt").exists()
