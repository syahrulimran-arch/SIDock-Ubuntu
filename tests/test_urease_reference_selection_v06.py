"""4UBP regression: modified KCX must never outrank the true HAE inhibitor."""
from pathlib import Path

from sidock_ubuntu.intake import (
    choose_reference_ligand, eligible_reference_groups, reference_candidates,
    semantic_reference_groups,
)
from sidock_ubuntu.core import clean_receptor_pdb


def atom(serial, name, residue, chain, number, element):
    return (
        f"HETATM{serial:5d} {name:<4s} {residue:>3s} {chain}{number:4d}    "
        f"{28.0 + serial * .01:8.3f}{70.0:8.3f}{63.0:8.3f}"
        f"{1.00:6.2f}{19.00:6.2f}          {element:>2s}\n"
    )


def _write_sample(path: Path):
    lines = ["HEADER    HYDROLASE                               25-FEB-99   4UBP\n",
             "MODRES 4UBP KCX C  220  LYS  LYSINE NZ-CARBOXYLIC ACID\n"]
    backbone = ["N", "CA", "CB", "CG", "CD", "CE", "NZ", "C", "O", "CX", "OQ1", "OQ2"]
    for i, name in enumerate(backbone, 1):
        lines.append(atom(i, name, "KCX", "C", 220, "O" if name.startswith("O") else ("N" if name in ("N", "NZ") else "C")))
    for j, name in enumerate(("C1", "C2", "N1", "O1", "O2"), 20):
        lines.append(atom(j, name, "HAE", "C", 800, name[0]))
    lines.extend((atom(41, "NI", "NI", "C", 798, "NI"),
                  atom(42, "NI", "NI", "C", 799, "NI"),
                  atom(43, "C", "ACE", "A", 0, "C"),
                  atom(44, "O", "ACE", "A", 0, "O")))
    path.write_text("".join(lines), encoding="utf-8")


def test_true_hydroxamate_selected_while_modified_lysine_excluded(tmp_path):
    pdb = tmp_path / "4ubp.pdb"
    _write_sample(pdb)
    groups = eligible_reference_groups(pdb)
    assert [g.key for g in groups] == ["HAE:C:800"]
    assert choose_reference_ligand(pdb).key == "HAE:C:800"
    assert [g.key for g in semantic_reference_groups(pdb)] == ["HAE:C:800"]
    assert [r["key"] for r in reference_candidates(pdb)] == ["HAE:C:800"]


def test_receptor_cleaning_preserves_urease_active_centre(tmp_path):
    pdb = tmp_path / "4ubp.pdb"
    _write_sample(pdb)
    clean = tmp_path / "clean.pdb"
    clean_receptor_pdb(pdb, clean, "HAE", "C", "800")
    result = clean.read_text()
    assert "HAE C 800" not in result
    assert "KCX C 220" in result
    assert " NI C 798" in result
    assert " NI C 799" in result


def test_backbone_excluded_even_without_modres(tmp_path):
    pdb = tmp_path / "4ubp.pdb"
    _write_sample(pdb)
    pdb.write_text("\n".join(line for line in pdb.read_text().splitlines() if not line.startswith("MODRES")) + "\n")
    assert choose_reference_ligand(pdb).key == "HAE:C:800"


def test_reference_correction_dry_run_does_not_mutate_files(tmp_path):
    import json
    from sidock_ubuntu.reference_repair import repair_reference
    project = tmp_path / "case"
    (project / "input").mkdir(parents=True)
    pdb = project / "input/protein.pdb"
    _write_sample(pdb)
    (project / "project.json").write_text(json.dumps({"receptor_source": str(pdb)}))
    (project / "project_state.json").write_text(json.dumps({"stages": {"REFERENCE": {"status": "PASS"}}, "compounds": {}}))
    old_bytes = (project / "project_state.json").read_bytes()
    data = repair_reference(project, expected="HAE:C:800")
    assert data["proposed_reference"] == "HAE:C:800"
    assert data["applied"] is False
    assert (project / "project_state.json").read_bytes() == old_bytes
    assert not (project / "provenance").exists()


def test_explicit_reference_override_accepts_unknown_ccd_mapping():
    # Missing optional CCD is not evidence of an invalid crystallographic inhibitor.
    from sidock_ubuntu.reference import select_reference_instance
    candidate = {
        "key": "HAE:C:800", "ligand_code": "HAE", "chain": "C",
        "residue_number": "800", "mapping_fraction": None,
        "atom_count": 5, "heavy_atom_count": 5,
    }
    selected = select_reference_instance([candidate], override="HAE:C:800")
    assert selected["selected"]["key"] == "HAE:C:800"
    assert selected["selected"]["mapping_fraction"] is None
    assert selected["selected"]["selection_method"] == "EXPLICIT_OVERRIDE"


def test_explicit_reference_override_rejects_confirmed_zero_mapping():
    import pytest
    from sidock_ubuntu.reference import select_reference_instance
    candidate = {
        "key": "HAE:C:800", "ligand_code": "HAE", "chain": "C",
        "residue_number": "800", "mapping_fraction": 0.0,
    }
    with pytest.raises(RuntimeError, match="INVALID_REFERENCE_OVERRIDE"):
        select_reference_instance([candidate], override="HAE:C:800")
