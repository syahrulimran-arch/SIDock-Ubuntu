"""Receptor errors must be actionable and must not alter the prepared PDB."""
import json
from pathlib import Path

from sidock_ubuntu.receptor_diagnostics import analyze_meeko_log, read_receptor_diagnostics

MEEKO_EXCERPT = """
/.venv/lib/python3.13/site-packages/meeko/polymer.py:1237: RuntimeWarning:
Input residues {'C:798': 'NI', 'C:799': 'NI'} not in residue_templates
Molecule contains unexpected elements -> template for NI will be None.
Error: Creation of data structure for receptor failed.
Template generation failed for unknown residues: {'NI'}
No template matched for residue_key='B:111'
No template matched for residue_key='C:403'
Error: Creation of data structure for receptor failed.
Details:
- Template matching failed for: ['A:0', 'B:111', 'B:119', 'C:326', 'C:403']
"""


def test_parser_distinguishes_nickel_and_unmatched_protein_residues():
    result = analyze_meeko_log(MEEKO_EXCERPT)
    assert result["unknown_components"] == ["NI"]
    assert result["template_mismatches"] == [
        "A:0", "B:111", "B:119", "C:326", "C:403"
    ]
    assert any("Nickel" in note for note in result["warnings"])
    assert any("ACE" in note for note in result["warnings"])
    assert any("--allow_bad_res" in note for note in result["warnings"])


def test_current_project_log_read_only(tmp_path):
    project = tmp_path / "test_project"
    path = project / "manual_docking/logs/receptor_meeko.log"
    path.parent.mkdir(parents=True)
    path.write_text(MEEKO_EXCERPT, encoding="utf-8")
    before = path.read_bytes()
    result = read_receptor_diagnostics(project)
    assert result["available"]
    assert "NI" in result["unknown_components"]
    assert before == path.read_bytes()


def test_missing_meeko_log_not_misdiagnosed(tmp_path):
    result = read_receptor_diagnostics(tmp_path)
    assert result["available"] is False
    assert "unknown_components" not in result



def test_exploratory_failure_reads_its_own_log_not_old_standard_failure(tmp_path):
    folder = tmp_path / "manual_docking/logs"
    folder.mkdir(parents=True)
    (folder/"receptor_meeko.log").write_text(
        "Template generation failed for unknown residues: {'NI'}\n",
        encoding="utf-8",
    )
    (folder/"receptor_meeko_metal_free.log").write_text(
        "Template matching failed for: ['A:0', 'C:326']\n",
        encoding="utf-8",
    )
    first = read_receptor_diagnostics(tmp_path)
    second = read_receptor_diagnostics(tmp_path, exploratory_metal_free=True)
    assert first["unknown_components"] == ["NI"]
    assert second["unknown_components"] == []
    assert second["template_mismatches"] == ["A:0", "C:326"]
