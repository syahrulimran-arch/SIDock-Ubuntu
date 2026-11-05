"""Ni reconstruction tests: no actual Meeko, AutoGrid or AutoDock executions.

These use fixture-only force-field numbers to verify file wiring and safeguards.
No fixture parameter values are scientific recommendations.
"""
import json
from pathlib import Path

import pytest

from sidock_ubuntu import core
from sidock_ubuntu.guided_workflow import (
    create_guided_project, select_protein, approve, approved,
)
from sidock_ubuntu.manual_docking import install_parameter_file, root
from sidock_ubuntu.nickel_reconstruction import (
    prepare_intermediate, reconstruct, parse_charges, sites,
)


def pdb_line(serial, name, code, chain, resnum, xyz, element="C"):
    x,y,z = xyz
    return (
        f"HETATM{serial:5d} {name:>4s} {code:>3s} {chain}{resnum:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2s}\n"
    )


def build_project(tmp_path):
    project = create_guided_project(tmp_path, "urease-nickel-test")
    pdb = project / "input/protein.pdb"
    pdb.write_text(
        "HEADER    HYDROLASE                               25-FEB-99   4UBP\n"
        + pdb_line(1, "N", "KCX", "C", 220, (20., 70., 62.), "N")
        + pdb_line(2, "CA", "KCX", "C", 220, (21., 70., 62.), "C")
        + pdb_line(3, "CX", "KCX", "C", 220, (21., 72., 70.), "C")
        + pdb_line(4, "OQ1", "KCX", "C", 220, (22., 72., 71.), "O")
        + pdb_line(5, "OQ2", "KCX", "C", 220, (20., 72., 71.), "O")
        + pdb_line(6, "C1", "HAE", "C", 800, (29., 72., 75.), "C")
        + pdb_line(7, "C2", "HAE", "C", 800, (30., 71., 75.), "C")
        + pdb_line(8, "N", "HAE", "C", 800, (31., 72., 75.), "N")
        + pdb_line(9, "O", "HAE", "C", 800, (30., 72., 74.), "O")
        + pdb_line(10, "NI", "NI", "C", 798, (28.500, 73.200, 72.000), "NI")
        + pdb_line(11, "NI", "NI", "C", 799, (31.250, 73.400, 73.500), "NI")
        + "END\n"
    )
    choice = select_protein(
        project, ligand_key="HAE:C:800", remove_groups=[], remove_waters=True
    )
    assert choice["mode"] == "CRYSTAL_LIGAND"
    return project


def fake_meeko(source, destination, log, *, allow_bad_res):
    assert allow_bad_res is False
    content=Path(source).read_text()
    assert "HAE C 800" not in content
    assert " NI C 798" not in content
    assert " NI C 799" not in content
    assert "KCX C 220" in content
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        "ATOM      1  C   GLY C  40      30.000  70.000  72.000  1.00  0.00    +0.100 C\n"
        "ATOM      2  OQ1 KCX C 220      28.237  72.536  71.006  1.00  0.00    -0.528 OA\n"
        "END\n"
    )
    log = Path(log)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text("Files written.\n")
    return destination


def example_ad4_full(path):
    # Synthetic parameters for fixture only: MUST NOT use in scientific docking.
    atom_values = {
        "C": "4.00 0.15 33.51 -0.00143 0 0 0 -1 -1 0",
        "N": "3.50 0.16 22.44 -0.00162 0 0 0 -1 -1 1",
        "OA": "3.20 0.20 17.15 -0.00251 1.9 5.0 5 -1 -1 2",
        "HD": "2.00 0.02 0.00 0.00051 0 0 2 -1 -1 3",
        "Ni": "1.00 0.01 1.00 -0.001 0 0 0 -1 -1 4",
    }
    path.write_text(
        "\n".join(f"atom_par {name} {values}" for name,values in atom_values.items())
        + "\n", encoding="utf-8"
    )
    return path


def prepare(project, monkeypatch):
    monkeypatch.setattr(core, "prepare_receptor", fake_meeko)
    first = prepare_intermediate(project)
    assert first["status"] == "PREPARED"
    assert first["production_receptor_installed"] is False
    assert not (root(project)/"receptor.pdbqt").exists()
    assert (root(project)/"inputs/protein_openbabel_Ni_intermediate.pdbqt").is_file()
    assert prepare_intermediate(project)["status"] == "REUSED"
    return first


def test_metal_intermediate_never_installs_metal_free_receptor(tmp_path, monkeypatch):
    project = build_project(tmp_path)
    original = (project/"input/protein.pdb").read_bytes()
    cleaned = (root(project)/"inputs/protein_clean.pdb").read_bytes()
    prepare(project, monkeypatch)
    assert original == (project/"input/protein.pdb").read_bytes()
    assert cleaned == (root(project)/"inputs/protein_clean.pdb").read_bytes()
    assert not (root(project)/"runs").exists()


def test_charges_are_explicit_and_complete():
    assert parse_charges("C:798=1.234,C:799=1.111") == {"C:798": 1.234, "C:799": 1.111}
    for invalid in ("", "C:798=...", "C:798=nan", "C:798=10", "C:798=1;C:798=2"):
        with pytest.raises(ValueError):
            parse_charges(invalid)


def test_reconstruction_requires_user_parameters_and_provenance(tmp_path, monkeypatch):
    project = build_project(tmp_path)
    prepare(project, monkeypatch)
    charges = parse_charges("C:798=1.100,C:799=1.200")
    with pytest.raises(RuntimeError, match="acknowledge"):
        reconstruct(project, atom_type="Ni", charges=charges,
                    parameter_source="fixture", charge_source="fixture")
    with pytest.raises(RuntimeError, match="parameter library"):
        reconstruct(project, atom_type="Ni", charges=charges,
                    parameter_source="fixture", charge_source="fixture",
                    nonbonded_model_acknowledged=True)
    assert not (root(project)/"receptor.pdbqt").exists()


def test_crystal_ni_coordinates_and_ad4_type_restored_and_audited(tmp_path, monkeypatch):
    project = build_project(tmp_path)
    original = (project/"input/protein.pdb").read_bytes()
    clean = (root(project)/"inputs/protein_clean.pdb").read_bytes()
    prepare(project, monkeypatch)
    parameter_file = example_ad4_full(project/"input/test_model.dat")
    install_parameter_file(project, parameter_file)

    result = reconstruct(
        project, atom_type="Ni",
        charges=parse_charges("C:798=1.100,C:799=1.200"),
        parameter_source="fixture-only force-field test, not scientifically valid",
        charge_source="fixture-only charges, not scientifically valid",
        nonbonded_model_acknowledged=True,
    )
    assert result["status"] == "RECONSTRUCTED_NOT_VALIDATED"
    receptor = (root(project)/"receptor.pdbqt").read_text()
    assert receptor.index(" NI C 798") < receptor.index("END")
    assert receptor.index(" NI C 799") < receptor.index("END")
    metal_lines = [l for l in receptor.splitlines()
                   if l.startswith("HETATM") and l[17:20].strip() == "NI"]
    assert len(metal_lines) == 2
    assert [l[77:79].strip() for l in metal_lines] == ["Ni", "Ni"]
    assert [float(l[70:76]) for l in metal_lines] == [1.1, 1.2]
    assert [tuple(float(l[i:i+8]) for i in (30,38,46)) for l in metal_lines] == [
        (28.500, 73.200, 72.000), (31.250, 73.400, 73.500)
    ]
    audit = json.loads((root(project)/"inputs/Ni_reconstruction_audit.json").read_text())
    assert audit["native_redocking"] == "NOT_PERFORMED"
    assert audit["NI_count"] == 2
    assert audit["scientific_review_required"] is True
    assert original == (project/"input/protein.pdb").read_bytes()
    assert clean == (root(project)/"inputs/protein_clean.pdb").read_bytes()
    assert not (root(project)/"runs").exists()
    from sidock_ubuntu.manual_docking import settings
    assert settings(project)["allow_missing_metal"] is False
    approve(project,2)
    assert approved(project,2)

    # Changing the custom force field invalidates protein approval.
    (root(project)/"custom_parameters.dat").write_text(
        (root(project)/"custom_parameters.dat").read_text() + "# mutation\n"
    )
    assert not approved(project,2)


def test_reconstruction_cannot_use_missing_or_extra_site_charges(tmp_path, monkeypatch):
    project = build_project(tmp_path)
    prepare(project, monkeypatch)
    install_parameter_file(project, example_ad4_full(project/"input/test_model.dat"))
    with pytest.raises(ValueError, match="exactly one charge"):
        reconstruct(
            project, atom_type="Ni", charges={"C:798":1.2},
            parameter_source="fixture", charge_source="fixture",
            nonbonded_model_acknowledged=True,
        )
