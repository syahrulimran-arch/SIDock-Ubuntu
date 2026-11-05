"""Rigid receptor PDBQT export without Meeko.

This module produces an inspectable candidate file. It NEVER silently approves,
installs, models missing atoms, deletes metals, or declares metal charges valid.
"""
from __future__ import annotations

import collections
import hashlib
import json
import math
import shutil
from pathlib import Path
from datetime import datetime, timezone

def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def _atoms(path: Path, *, input_pdb: bool) -> list[dict]:
    records = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 54:
            raise RuntimeError(f"Malformed atom coordinates in {path.name}: {line[:30]}")
        try:
            xyz = tuple(float(line[i:i + 8]) for i in (30, 38, 46))
        except ValueError as exc:
            raise RuntimeError(f"Invalid PDB/PDBQT coordinates in {path.name}") from exc
        if not all(math.isfinite(v) for v in xyz):
            raise RuntimeError(f"Nonfinite coordinates in {path.name}")
        chain = line[21:22].strip() or "-"
        residue = line[17:20].strip()
        number = line[22:27].strip()
        atom = line[12:16].strip()
        if input_pdb:
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not element:
                raise RuntimeError(
                    f"Missing element column for {residue}:{chain}:{number}:{atom}. "
                    "Correct the PDB element metadata before audit."
                )
            typ = element
            charge = None
        else:
            if len(line) < 79:
                raise RuntimeError(f"Truncated PDBQT atom line: {line[:50]}")
            typ = line[77:79].strip()
            try:
                charge = float(line[70:76])
            except ValueError as exc:
                raise RuntimeError(f"Unreadable PDBQT charge for {atom}") from exc
            if not typ or not math.isfinite(charge):
                raise RuntimeError(f"Invalid atom type/charge on PDBQT record {atom}")
            element = ""
        records.append({
            "identity": (chain, residue, number, atom),
            "xyz": xyz, "element": element, "atom_type": typ,
            "charge": charge,
        })
    return records

def generate_candidate(project: Path) -> dict:
    """Reviewable Open Babel export; no Meeko, no job approval or docking."""
    from .guided_workflow import _json, _sha_file, _save_json, _event
    from .manual_docking import root
    from .core import run_command

    project = Path(project).resolve()
    home = root(project)
    src = home / "inputs/protein_clean.pdb"
    choice = _json(home / "inputs/protein_choice.json")
    if not choice or not src.is_file() or _sha_file(src) != choice.get("clean_sha256"):
        raise RuntimeError("Save/review the Step 2 protein first. The cleaned PDB hash must match its review.")
    obabel = shutil.which("obabel")
    if not obabel:
        raise RuntimeError("Open Babel executable 'obabel' not found on Ubuntu PATH.")
    source_atoms = _atoms(src, input_pdb=True)
    if not source_atoms:
        raise RuntimeError("Reviewed PDB has no atoms.")
    output = home / "inputs/receptor_openbabel_candidate.pdbqt"
    audit_path = home / "inputs/receptor_openbabel_audit.json"
    log = home / "logs/receptor_openbabel.log"
    existing = _json(audit_path)
    if output.is_file():
        if existing.get("input_sha256") == _sha(src) and existing.get("candidate_sha256") == _sha(output):
            return {**existing, "status": "REUSED", "candidate": str(output), "audit": str(audit_path)}
        raise RuntimeError(
            f"Existing Open Babel receptor candidate differs: {output}. "
            "Preserved it; archive/review explicitly before regenerating."
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    # -xr rigid; -xc combine disconnected pieces; -xn preserve atom/residue names;
    # -xp atom indices; -xh preserve existing H. No automatic protonation or
    # deletion, and no biologically unjustified charges are filled by SIDock.
    cmd = [
        obabel, str(src.resolve()), "-O", str(output.resolve()),
        "-xr", "-xc", "-xn", "-xp", "-xh",
        "--partialcharge", "gasteiger",
    ]
    exit_code = run_command(cmd, output.parent, log)
    if exit_code != 0 or not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(
            f"Open Babel rigid receptor conversion failed (exit={exit_code}); "
            f"see {log}. A partial output, if present, is NOT approved."
        )
    generated = _atoms(output, input_pdb=False)
    if not generated:
        raise RuntimeError(f"Open Babel output contains no atoms: {output}")
    all_lines = output.read_text(encoding="utf-8", errors="replace").splitlines()
    if any(line.startswith(("ROOT", "BRANCH", "TORSDOF")) for line in all_lines):
        raise RuntimeError(f"Output has ligand torsion-tree records, not a rigid receptor: {output}")
    # Compare identifiers and coordinates. Hydrogen assignment is kept separate
    # because Open Babel may omit/add H; all original non-H atoms must be traced.
    def snap(atom: dict) -> tuple:
        return (*atom["identity"], *(round(v, 3) for v in atom["xyz"]))
    input_heavy = [a for a in source_atoms if a["element"] not in {"H", "D"}]
    output_snapshots = collections.Counter(snap(a) for a in generated)
    missing = [a for a in input_heavy if output_snapshots[snap(a)] == 0]
    for atom in input_heavy:
        key = snap(atom)
        if output_snapshots[key]:
            output_snapshots[key] -= 1
    original_metals = [a for a in input_heavy if a["element"] in {
        "NI","ZN","CU","FE","MN","CO","MG","CA","CD","HG","NA","K",
    }]
    missing_metal = [a for a in original_metals if any(
        a["identity"] == m["identity"] and a["xyz"] == m["xyz"] for m in missing
    )]
    types = sorted({a["atom_type"] for a in generated})
    audit = {
        "status": "CANDIDATE_REQUIRES_REVIEW" if not missing else "CANDIDATE_MISSING_ATOMS",
        "engine": "OPEN_BABEL",
        "command": cmd,
        "input_sha256": _sha(src),
        "candidate_sha256": _sha(output),
        "input_heavy_atom_count": len(input_heavy),
        "candidate_atom_count": len(generated),
        "missing_input_heavy_atoms": [
            {"chain": a["identity"][0], "resname": a["identity"][1],
             "resid": a["identity"][2], "atom": a["identity"][3], "xyz": a["xyz"]}
            for a in missing
        ],
        "metal_atoms_input": len(original_metals),
        "metal_atoms_missing": len(missing_metal),
        "atom_types": types,
        "total_pdbqt_partial_charge_e": round(sum(a["charge"] for a in generated), 4),
        "metal_parameterization": "NOT_VALIDATED",
        "protonation_review": "NOT_PERFORMED",
        "scientific_validation": "NOT_PERFORMED",
        "approved_for_docking": False,
        "log": str(log),
        "created": datetime.now(timezone.utc).isoformat(),
    }
    _save_json(audit_path, audit)
    _event(project, "RECEPTOR_OPENBABEL_CANDIDATE_GENERATED", {
        "candidate_sha256": audit["candidate_sha256"],
        "status": audit["status"], "missing_heavy_atoms": len(missing),
    })
    return {**audit, "candidate": str(output), "audit": str(audit_path)}
