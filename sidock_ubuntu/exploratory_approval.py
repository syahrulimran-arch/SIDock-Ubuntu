"""Hash-bound, explicitly acknowledged *technical* exploratory receptor approval.

The successful smoke test is technical evidence, not scientific metal-site
validation. The receptor is never installed over manual_docking/receptor.pdbqt.
No PDB ID, chain, residue index or number of nickel ions is assumed.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

MODE = "EXPLORATORY_UNVALIDATED_NI_MODEL"
ACKNOWLEDGEMENT = "I_ACKNOWLEDGE_UNVALIDATED_EXPLORATORY_DOCKING"


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _json(path: Path) -> dict:
    try:
        result = json.loads(path.read_text(encoding="utf-8"))
        return result if isinstance(result, dict) else {}
    except (OSError, ValueError):
        return {}


def readiness(project: Path) -> dict:
    """Check only recorded source provenance and a successful, matching smoke run.

    No source files, receptor PDBQTs or approval records are altered.
    """
    from .manual_docking import root, _input_metals, pdbqt_atom_types
    home = root(project)
    inp = home / "inputs"
    source = inp / "protein_clean.pdb"
    choice = inp / "protein_choice.json"
    structural = inp / "receptor_openbabel_candidate.pdbqt"
    charged = inp / "receptor_openbabel_protonated_charge_candidate.pdbqt"
    typed = inp / "receptor_openbabel_hydrogen_typed_candidate.pdbqt"
    candidate = inp / "receptor_openbabel_exploratory.pdbqt"
    params = home / "custom_parameters.dat"
    smoke_file = home / "exploratory_smoke/exploratory_smoke_audit.json"
    source_audit_file = inp / "receptor_openbabel_audit.json"
    charge_audit_file = inp / "receptor_openbabel_protonated_charge_audit.json"
    type_audit_file = inp / "receptor_openbabel_hydrogen_typed_audit.json"
    receptor_audit_file = inp / "receptor_openbabel_exploratory_audit.json"
    files = {
        "protein_choice": choice, "clean_pdb": source,
        "structural_pdbqt": structural, "charged_pdbqt": charged,
        "typed_pdbqt": typed, "exploratory_receptor": candidate,
        "parameters": params, "smoke": smoke_file,
        "structural_audit": source_audit_file,
        "charge_audit": charge_audit_file,
        "hydrogen_audit": type_audit_file,
        "receptor_audit": receptor_audit_file,
    }
    reasons = []
    for key, path in files.items():
        if not path.is_file() or path.stat().st_size == 0:
            reasons.append(f"Missing {key}: {path.name}")
    if reasons:
        return {"ready": False, "reasons": reasons,
                "candidate": str(candidate), "mode": MODE}
    source_hash = _sha(source)
    candidate_hash = _sha(candidate)
    parameter_hash = _sha(params)
    saved_choice = _json(choice)
    structural_audit = _json(source_audit_file)
    charged_audit = _json(charge_audit_file)
    hydrogen_audit = _json(type_audit_file)
    receptor_audit = _json(receptor_audit_file)
    smoke = _json(smoke_file)
    if saved_choice.get("clean_sha256") != source_hash:
        reasons.append("Current Step 2 cleaned PDB no longer matches its reviewed hash.")
    if structural_audit.get("input_sha256") != source_hash or \
            structural_audit.get("candidate_sha256") != _sha(structural):
        reasons.append("Structural PDBQT source/export hashes do not match.")
    if structural_audit.get("missing_input_heavy_atoms") or structural_audit.get("metal_atoms_missing"):
        reasons.append("Structural PDBQT audit reports omitted original heavy/metal atoms.")
    if charged_audit.get("source_clean_pdb_sha256") != source_hash or \
            charged_audit.get("source_structural_candidate_sha256") != _sha(structural) or \
            charged_audit.get("candidate_sha256") != _sha(charged):
        reasons.append("Charged PDBQT provenance no longer matches the source.")
    if hydrogen_audit.get("source_candidate_sha256") != _sha(charged) or \
            hydrogen_audit.get("candidate_sha256") != _sha(typed) or \
            hydrogen_audit.get("heavy_atoms_preserved_exactly") is not True:
        reasons.append("Hydrogen-typed PDBQT audit is missing/stale or lost heavy atoms.")
    if receptor_audit.get("source_candidate_sha256") != _sha(typed) or \
            receptor_audit.get("candidate_sha256") != candidate_hash or \
            receptor_audit.get("parameter_file_sha256") != parameter_hash or \
            receptor_audit.get("heavy_atoms_preserved_exactly") is not True or \
            receptor_audit.get("missing_parameter_types"):
        reasons.append("Exploratory receptor provenance, heavy atoms or AD4 atom types failed review.")
    original_metals = _input_metals(project)
    if original_metals and receptor_audit.get("ni_count") != original_metals:
        reasons.append("Original nickel-ion count differs from the exploratory receptor audit.")
    if receptor_audit.get("scientific_approval") is not False or \
            receptor_audit.get("technical_exploratory_mode") is not True:
        reasons.append("Receptor audit does not identify an unvalidated exploratory model.")
    if smoke.get("status") != "PASSED" or smoke.get("scientific_approval") is not False:
        reasons.append("A successful, explicitly unvalidated AutoGrid/AutoDock smoke audit is required.")
    smoke_identity = smoke.get("identity", {})
    smoke_approval = smoke.get("technical_exploratory_approval", {})
    if smoke_identity.get("receptor_sha256") != candidate_hash or \
            smoke_identity.get("parameter_sha256") != parameter_hash or \
            smoke_identity.get("mode") != "TECHNICAL_EXPLORATORY_ONLY":
        reasons.append("Smoke-test receptor/force-field hashes or run mode do not match.")
    if smoke_approval.get("status") != "APPROVED_FOR_TECHNICAL_SMOKE_ONLY" or \
            smoke_approval.get("receptor_sha256") != candidate_hash or \
            smoke_approval.get("parameter_sha256") != parameter_hash:
        reasons.append("Technical smoke-test approval does not match this receptor/library.")
    if smoke.get("atom_type_errors") or not smoke.get("grid", {}).get("maps_exist") or \
            not smoke.get("docking", {}).get("successful_completion"):
        reasons.append("Smoke test has atom-type, grid-map or docking completion errors.")
    # The atom type list must be representable by the *same* user-supplied AD4 library.
    if not reasons:
        declared = {
            fields[1] for line in params.read_text(encoding="utf-8", errors="replace").splitlines()
            if (fields := line.split()) and len(fields) >= 12 and fields[0] == "atom_par"
        }
        needed = set(pdbqt_atom_types(candidate))
        if not needed or needed - declared:
            reasons.append("Candidate atom types are absent from the bound custom AD4 library.")
    hashes = {name: _sha(path) for name, path in files.items()}
    evidence = {
        "mode": MODE,
        "clean_sha256": source_hash,
        "receptor_sha256": candidate_hash,
        "parameter_sha256": parameter_hash,
        "evidence_hashes": hashes,
        "scientific_approval": False,
    }
    digest = hashlib.sha256(json.dumps(evidence, sort_keys=True).encode()).hexdigest()
    return {
        "ready": not reasons, "reasons": reasons,
        "mode": MODE, "candidate": str(candidate),
        "receptor_sha256": candidate_hash,
        "parameter_sha256": parameter_hash,
        "digest": digest if not reasons else None,
        "scientific_approval": False,
    }
