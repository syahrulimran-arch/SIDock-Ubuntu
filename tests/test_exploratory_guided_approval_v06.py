"""In-memory fixture verifies Step 2 -> Step 3 -> Step 4 integration.

The fixture is a synthetic receptor and audited *technical* smoke result.
These tests do not run AutoGrid, AutoDock, Open Babel or scientific validation.
"""
from __future__ import annotations

import hashlib
import json

import pytest
from rdkit import Chem

from sidock_ubuntu.exploratory_approval import (
    ACKNOWLEDGEMENT, MODE, readiness,
)
from sidock_ubuntu.guided_workflow import (
    approve, approve_exploratory, approved, approval_mode,
    create_guided_project, load_workflow, select_protein,
)
from sidock_ubuntu.guided_ui import step2, step3, step4, step5
from sidock_ubuntu.manual_docking import (
    _run_identity, root, selected_receptor_for_run,
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True))


def pdb_atom(serial, atom, residue, chain, number, xyz, element):
    x, y, z = xyz
    return (f"HETATM{serial:5d} {atom:>4s} {residue:>3s} {chain}{number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2s}\n")


def pdbqt_atom(serial, atom, residue, chain, number, xyz, element, ad_type):
    line = pdb_atom(serial, atom, residue, chain, number, xyz, element)
    return line[:66].ljust(70) + f"{0.0:6.3f}" + f" {ad_type:>2s}\n"


def project_with_exploratory_evidence(tmp_path):
    p = create_guided_project(tmp_path, "other-protein")
    original = p / "input/protein.pdb"
    original.write_text(
        "HEADER    SYNTHETIC TEST RECEPTOR\n"
        + pdb_atom(1, "CA", "ALA", "A", 1, (4, 5, 6), "C")
        + pdb_atom(2, "NI", "NI", "A", 201, (7, 8, 9), "NI")
        + "END\n"
    )
    select_protein(
        p, ligand_key=None, remove_groups=[],
        center=(5.0, 6.0, 7.0), remove_waters=False
    )
    home = root(p)
    inp = home / "inputs"
    clean = inp / "protein_clean.pdb"
    # This file is deliberately different; exploratory approval must not
    # overwrite an unrelated scientifically reviewed manual receptor.
    scientific = home / "receptor.pdbqt"
    scientific.write_text(
        pdbqt_atom(1, "CA", "ALA", "A", 1, (1, 1, 1), "C", "C")
    )
    scientific_bytes = scientific.read_bytes()

    pdbqt = (
        pdbqt_atom(1, "CA", "ALA", "A", 1, (4, 5, 6), "C", "C")
        + pdbqt_atom(2, "NI", "NI", "A", 201, (7, 8, 9), "NI", "Ni")
    )
    structural = inp / "receptor_openbabel_candidate.pdbqt"
    charged = inp / "receptor_openbabel_protonated_charge_candidate.pdbqt"
    typed = inp / "receptor_openbabel_hydrogen_typed_candidate.pdbqt"
    exploratory = inp / "receptor_openbabel_exploratory.pdbqt"
    for path in (structural, charged, typed, exploratory):
        path.write_text(pdbqt)
    params = home / "custom_parameters.dat"
    params.write_text(
        "atom_par C 1 1 1 1 1 1 1 1 1 1\n"
        "atom_par Ni 1 1 1 1 1 1 1 1 1 1\n"
    )
    write_json(inp / "receptor_openbabel_audit.json", {
        "input_sha256": sha(clean),
        "candidate_sha256": sha(structural),
        "missing_input_heavy_atoms": [],
        "metal_atoms_missing": 0,
    })
    write_json(inp / "receptor_openbabel_protonated_charge_audit.json", {
        "source_clean_pdb_sha256": sha(clean),
        "source_structural_candidate_sha256": sha(structural),
        "candidate_sha256": sha(charged),
    })
    write_json(inp / "receptor_openbabel_hydrogen_typed_audit.json", {
        "source_candidate_sha256": sha(charged),
        "candidate_sha256": sha(typed),
        "heavy_atoms_preserved_exactly": True,
    })
    write_json(inp / "receptor_openbabel_exploratory_audit.json", {
        "source_candidate_sha256": sha(typed),
        "candidate_sha256": sha(exploratory),
        "parameter_file_sha256": sha(params),
        "heavy_atoms_preserved_exactly": True,
        "missing_parameter_types": [],
        "ni_count": 1,
        "scientific_approval": False,
        "technical_exploratory_mode": True,
    })
    smoke = home / "exploratory_smoke/exploratory_smoke_audit.json"
    gridlog = home / "exploratory_smoke/grid/production/autogrid.glg"
    dlg = home / "exploratory_smoke/docking/dlg/ligand.dlg"
    gridmap = home / "exploratory_smoke/grid/production/receptor.C.map"
    for path, content in (
        (gridlog, "autogrid4: Successful Completion\n"),
        (dlg, "autodock4: Successful Completion\n"),
        (gridmap, "synthetic map\n"),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    write_json(smoke, {
        "status": "PASSED",
        "scientific_approval": False,
        "identity": {
            "receptor_sha256": sha(exploratory),
            "parameter_sha256": sha(params),
            "mode": "TECHNICAL_EXPLORATORY_ONLY",
        },
        "technical_exploratory_approval": {
            "status": "APPROVED_FOR_TECHNICAL_SMOKE_ONLY",
            "receptor_sha256": sha(exploratory),
            "parameter_sha256": sha(params),
        },
        "grid": {"maps_exist": True, "glg": str(gridlog), "maps": [str(gridmap)]},
        "docking": {"successful_completion": True, "dlg": str(dlg)},
        "atom_type_errors": [],
    })
    # Add already-prepared library without asking Step 3 to re-minimize.
    sdf = p / "input/compounds.sdf"
    molecule = Chem.MolFromSmiles("CCO")
    molecule.SetProp("_Name", "Z-1")
    molecule.SetProp("Compound_ID", "Z-1")
    with Chem.SDWriter(str(sdf)) as writer:
        writer.write(molecule)
    ligand = home / "ligands/prepared/Z-1.pdbqt"
    ligand.parent.mkdir(parents=True, exist_ok=True)
    ligand.write_text(pdbqt_atom(1, "C", "UNL", "X", 1, (3, 3, 3), "C", "C"))
    write_json(home / "ligands/manifest.json", {
        "source_sha256": sha(sdf),
        "compounds": ["Z-1.pdbqt"],
    })
    return p, scientific_bytes, exploratory, params, smoke


def test_full_guided_exploratory_approval_and_ligand_reuse(tmp_path):
    p, scientific_bytes, exploratory, params, smoke = project_with_exploratory_evidence(tmp_path)
    assert readiness(p)["ready"]
    assert "Approve for exploratory docking and continue" in step2(p)
    with pytest.raises(RuntimeError, match="Explicitly acknowledge"):
        approve_exploratory(p, acknowledgement="")

    approval = approve_exploratory(p, acknowledgement=ACKNOWLEDGEMENT)
    assert approval["mode"] == MODE
    assert approval["scientific_approval"] is False
    assert approved(p, 2)
    assert approval_mode(p) == MODE
    assert (root(p) / "receptor.pdbqt").read_bytes() == scientific_bytes

    assert "EXPLORATORY_UNVALIDATED_NI_MODEL" in step3(p)
    assert "Approve and continue" in step3(p)
    assert "Minimize and prepare compounds" not in step3(p)
    approve(p, 3)
    assert approved(p, 3)

    receptor_for_run, mode = selected_receptor_for_run(p)
    assert receptor_for_run == exploratory
    assert mode == MODE
    token, ligands = _run_identity(p)
    assert len(token) == 16
    assert [x.name for x in ligands] == ["Z-1.pdbqt"]

    html = step4(p)
    assert "EXPLORATORY_UNVALIDATED_NI_MODEL" in html
    assert "Generate AutoGrid4 maps" in html
    assert "disabled>Generate AutoGrid4 maps" not in html
    assert "EXPLORATORY_UNVALIDATED_NI_MODEL" in step5(p)
    assert (root(p) / "receptor.pdbqt").read_bytes() == scientific_bytes
    assert "4UBP" not in html


def test_parameter_change_invalidates_technical_approval(tmp_path):
    p, _, exploratory, params, smoke = project_with_exploratory_evidence(tmp_path)
    approve_exploratory(p, acknowledgement=ACKNOWLEDGEMENT)
    assert approved(p, 2)
    params.write_text(params.read_text() + "# changed\n")
    assert not approved(p, 2)
    assert approval_mode(p) == "UNAPPROVED"
    assert not readiness(p)["ready"]
    with pytest.raises(RuntimeError, match="approval expired"):
        _run_identity(p)
    assert exploratory.exists()
    assert smoke.exists()


def test_failed_smoke_audit_cannot_unlock_workflow(tmp_path):
    p, _, _, _, smoke = project_with_exploratory_evidence(tmp_path)
    report = json.loads(smoke.read_text())
    report["status"] = "FAILED"
    write_json(smoke, report)
    assert not readiness(p)["ready"]
    with pytest.raises(RuntimeError, match="smoke"):
        approve_exploratory(p, acknowledgement=ACKNOWLEDGEMENT)
    assert not approved(p, 2)


def test_changed_receptor_invalidates_approval_and_does_not_mutate_scientific_model(tmp_path):
    p, scientific_bytes, exploratory, _, _ = project_with_exploratory_evidence(tmp_path)
    approve_exploratory(p, acknowledgement=ACKNOWLEDGEMENT)
    exploratory.write_text(exploratory.read_text() + "REMARK changed\n")
    assert not approved(p, 2)
    assert (root(p) / "receptor.pdbqt").read_bytes() == scientific_bytes
    with pytest.raises(RuntimeError, match="approval expired"):
        selected_receptor_for_run(p)
