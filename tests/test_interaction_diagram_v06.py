"""Synthetic contact fixtures validate figure wiring, not binding predictions."""
import csv
import json

import pytest
from rdkit import Chem

from sidock_ubuntu.interaction_diagram import draw_interaction_diagrams
from sidock_ubuntu.plip_runner import normalize_plip_interactions


def row_with_report(tmp_path, *, count=1, smiles="CCO", mapping="1:101,2:102,3:103"):
    report = tmp_path / "report.xml"
    report.write_text("<report><bindingsite><identifiers><hetid>LIG</hetid><chain>Z</chain>"
                      "<position>999</position><smiles>" + smiles + "</smiles></identifiers>"
                      "<mappings><smiles_to_pdb>" + mapping + "</smiles_to_pdb></mappings>"
                      "</bindingsite></report>")
    return {"compound_id": "test", "plip_status": "COMPLETE", "plip_reports": [str(report)],
            "interactions": [{"protein_residue_name": "ASP", "protein_residue_number": str(100 + i),
                              "protein_chain": "B", "interaction_type": "hydrogen_bond",
                              "ligand_atom": "103", "distance_A": "2.91"} for i in range(count)]}


def test_png_has_source_record_and_verified_atom_endpoints(tmp_path):
    row = row_with_report(tmp_path)
    figures = draw_interaction_diagrams(row, Chem.MolFromSmiles("CCO"), tmp_path)
    assert figures[0]["path"] == "figures/test_interactions_2d.png"
    assert (tmp_path / "test_interactions_2d.png").read_bytes().startswith(b"\x89PNG")
    data = json.loads((tmp_path / "test_interactions_2d.json").read_text())
    assert data["scientific_approval"] is False
    assert data["contacts"][0]["endpoint"] == "PLIP_ATOM_OR_GROUP"
    assert data["contacts"][0]["distance_A"] == "2.91"
    assert "Exploratory" in figures[0]["caption"]


def test_graph_mismatch_uses_explicit_ligand_level_lines(tmp_path):
    row = row_with_report(tmp_path, smiles="CCN")
    figures = draw_interaction_diagrams(row, Chem.MolFromSmiles("CCO"), tmp_path)
    data = json.loads((tmp_path / "test_interactions_2d.json").read_text())
    assert data["contacts"][0]["endpoint"] == "LIGAND_LEVEL"
    assert "no verified atom endpoint" in figures[0]["caption"]


def test_dense_contact_list_is_paginated_without_dropping_contacts(tmp_path):
    row = row_with_report(tmp_path, count=15)
    figures = draw_interaction_diagrams(row, Chem.MolFromSmiles("CCO"), tmp_path)
    assert len(figures) == 2
    assert (tmp_path / "test_interactions_2d_p2.png").is_file()
    data = json.loads((tmp_path / "test_interactions_2d.json").read_text())
    assert len(data["contacts"]) == 15
    assert {c["page"] for c in data["contacts"]} == {1, 2}


def test_successful_zero_contacts_still_gets_a_truthful_figure(tmp_path):
    row = row_with_report(tmp_path, count=0)
    assert len(draw_interaction_diagrams(row, Chem.MolFromSmiles("CCO"), tmp_path)) == 1
    data = json.loads((tmp_path / "test_interactions_2d.json").read_text())
    assert data["contacts"] == []
    row["plip_status"] = "UNAVAILABLE"
    with pytest.raises(RuntimeError, match="did not complete"):
        draw_interaction_diagrams(row, Chem.MolFromSmiles("CCO"), tmp_path)


def test_new_plip_pi_stacks_and_group_indices_are_preserved(tmp_path):
    report = tmp_path / "report.xml"
    report.write_text("<report><bindingsite><pi_stacks><pi_stack><reschain>A</reschain>"
                      "<restype>PHE</restype><resnr>12</resnr><centdist>4.17</centdist>"
                      "<prot_idx_list><idx>7</idx><idx>8</idx></prot_idx_list>"
                      "<lig_idx_list><idx>101</idx><idx>102</idx></lig_idx_list>"
                      "</pi_stack></pi_stacks></bindingsite></report>")
    path = normalize_plip_interactions([report], tmp_path / "contacts.csv")
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["interaction_type"] == "pi-stack"
    assert rows[0]["distance_A"] == "4.17"
    assert rows[0]["ligand_atom"] == "101,102"
