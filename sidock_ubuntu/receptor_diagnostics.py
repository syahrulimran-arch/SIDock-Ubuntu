"""Readable diagnostics for Meeko receptor failures.

This module only reads saved logs and reviewed protein components. It NEVER
alters receptor coordinates, deletes residues, or starts docking.
"""
from __future__ import annotations

import ast
import re
from pathlib import Path


def analyze_meeko_log(text: str) -> dict:
    text = str(text or "")
    mismatches = set()
    unknown_templates = set()
    for match in re.finditer(r"No template matched for residue_key=['\"]([^'\"]+)", text):
        mismatches.add(match.group(1))
    for match in re.finditer(r"Template matching failed for:\s*(\[[^\]\n]*\])", text):
        try:
            values = ast.literal_eval(match.group(1))
            mismatches.update(str(value) for value in values)
        except (ValueError, SyntaxError):
            pass
    for match in re.finditer(r"Template generation failed for unknown residues:\s*(\{[^}\n]*\})", text):
        try:
            values = ast.literal_eval(match.group(1))
            unknown_templates.update(str(value) for value in values)
        except (ValueError, SyntaxError):
            pass

    warnings = []
    if "NI" in unknown_templates or re.search(r"template for NI will be None", text):
        warnings.append(
            "Nickel ions are not parameterized by the standard Meeko receptor templates. "
            "Do not remove catalytic nickel merely to make the job pass. "
            "Use a metal-aware prepared receptor PDBQT with scientifically justified parameters, "
            "or explicitly select an exploratory metal-free model."
        )
    if "KCX" in unknown_templates or re.search(r"['\"]KCX['\"]\s*not in residue_templates", text):
        warnings.append(
            "KCX is a modified lysine. Verify its carbamyl atoms and partial charges after "
            "preparation; successful PDBQT output alone does not validate metal coordination."
        )
    if "ACE" in unknown_templates or any(x.startswith("A:0") for x in mismatches):
        warnings.append(
            "The N-terminal ACE/protein cap may require a template or careful terminal chemistry. "
            "Do not treat a possible covalent cap as a freely bound inhibitor."
        )
    if "HAE" in unknown_templates:
        warnings.append(
            "HAE is still in the receptor preparation input. If it is the selected crystallographic "
            "inhibitor, extract it as the reference and remove it from the receptor first."
        )
    if mismatches:
        warnings.append(
            "Meeko cannot match the listed protein residues to its templates. "
            "Inspect backbone connectivity/terminal chemistry; do not automatically delete internal "
            "residues or use --allow_bad_res without assessing the binding site."
        )
    if not warnings and text.strip():
        warnings.append(
            "Meeko failed or logged a warning; inspect the exact output below. "
            "No residues have been deleted automatically."
        )
    significant = [line for line in text.splitlines() if line.strip()]
    return {
        "template_mismatches": sorted(mismatches),
        "unknown_components": sorted(unknown_templates),
        "warnings": warnings,
        "last_log_lines": "\n".join(significant[-70:]),
    }


def read_receptor_diagnostics(project: Path, *, exploratory_metal_free: bool = False) -> dict:
    name = "receptor_meeko_metal_free.log" if exploratory_metal_free else "receptor_meeko.log"
    path = Path(project) / "manual_docking/logs" / name
    if not path.is_file():
        return {"available": False, "path": str(path)}
    text = path.read_text(encoding="utf-8", errors="replace")
    return {
        "available": True,
        "path": str(path),
        **analyze_meeko_log(text),
    }
