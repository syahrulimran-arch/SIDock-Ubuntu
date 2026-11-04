"""Manual AutoDock setup tests (never execute AutoGrid or AutoDock)."""
import hashlib
import json

import pytest

from sidock_ubuntu.manual_docking import (
    _run_identity, install_receptor, prepare, root, save_settings, settings,
)


def atom(serial, res, chain, number, x, y, z, element="C"):
    return (
        f"HETATM{serial:5d}  {element:<3s} {res:>3s} {chain}{number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2s}\n"
    )


def setup_project(tmp_path):
    p = tmp_path / "project"
    (p / "input").mkdir(parents=True)
    pdb = p / "input/protein.pdb"
    pdb.write_text(
        "HEADER    HYDROLASE\n"
        + atom(1, "HAE", "C", 800, 28.000, 70.000, 72.000)
        + atom(2, "HAE", "C", 800, 30.000, 72.000, 76.000)
        + atom(3, "NI", "C", 798, 30.000, 72.000, 73.000, "NI")
        + atom(4, "NI", "C", 799, 31.000, 73.000, 74.000, "NI")
    )
    sdf = p / "input/compounds.sdf"
    sdf.write_text("synthetic-fixture-only\n", encoding="utf-8")
    (p / "project.json").write_text(json.dumps({
        "receptor_source": str(pdb), "ligand_source": str(sdf),
    }))
    return p


def setup_manual_receptor_and_ligand(p):
    pdbqt = p / "input/experimental.pdbqt"
    pdbqt.write_text(
        "ATOM      1  CA  GLY C  40      30.000  70.000  72.000  1.00  0.00    +0.100 C\n",
        encoding="utf-8",
    )
    install_receptor(p, pdbqt)
    dst = root(p) / "ligands/prepared/L1.pdbqt"
    dst.parent.mkdir(parents=True)
    dst.write_text(
        "ROOT\nHETATM    1  C1  UNL A   1      29.000  72.000  74.000  1.00  0.00    +0.100 C\nENDROOT\nTORSDOF 0\n"
    )
    sdf = p / "input/compounds.sdf"
    (root(p) / "ligands/manifest.json").write_text(json.dumps({
        "source_sha256": hashlib.sha256(sdf.read_bytes()).hexdigest(),
        "compounds": ["L1.pdbqt"],
    }))


def test_grid_suggested_from_crystallographic_hae(tmp_path):
    p = setup_project(tmp_path)
    t = settings(p)
    assert (t["center_x"],t["center_y"],t["center_z"]) == (29.0,71.0,74.0)
    assert (t["npts_x"],t["npts_y"],t["npts_z"]) == (80,80,80)


def test_invalid_grid_points_cannot_be_saved(tmp_path):
    p = setup_project(tmp_path)
    t = settings(p)
    t["npts_x"] = 79
    with pytest.raises(ValueError, match="npts_x"):
        save_settings(p, t)
    assert not (root(p) / "settings.json").exists()


def test_missing_catalytic_metals_require_explicit_acknowledgment(tmp_path):
    p = setup_project(tmp_path)
    setup_manual_receptor_and_ligand(p)
    with pytest.raises(RuntimeError, match="nickel"):
        _run_identity(p)
    t = settings(p)
    t["allow_missing_metal"] = True
    save_settings(p, t)
    token, ligands = _run_identity(p)
    assert len(token) == 16
    assert [x.name for x in ligands] == ["L1.pdbqt"]


def test_manual_setup_reuses_prepared_ligands_without_reprocessing(tmp_path):
    p = setup_project(tmp_path)
    setup_manual_receptor_and_ligand(p)
    result = prepare(p)
    assert result["status"] == "REUSED"
    assert result["prepared"] == 1


def test_grid_changes_create_new_run_identity_without_overwriting(tmp_path):
    p = setup_project(tmp_path)
    setup_manual_receptor_and_ligand(p)
    t = settings(p)
    t["allow_missing_metal"] = True
    save_settings(p, t)
    original, _ = _run_identity(p)
    t["center_x"] += 3
    save_settings(p, t)
    changed, _ = _run_identity(p)
    assert original != changed
