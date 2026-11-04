"""Audited fixed-centre Ni preparation for manually vetted AutoDock4 studies.

Metal coordinates are copied from the reviewed protein PDB; they are never
optimized or invented. Meeko receives a separate *temporary* Ni-free protein,
but the intermediate is never approved as a production receptor. Reconstitution
requires a user-supplied complete AutoDock4 force-field library and individually
specified metal partial charges with provenance. This does not establish that a
nonbonded AutoDock4 model faithfully describes Ni coordination.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

from .manual_docking import root, parameter_file_path, install_parameter_file, pdbqt_atom_types
from .guided_workflow import _json, _sha_file, _save_json, install_reviewed_receptor, _event


def _reviewed_source(project: Path) -> tuple[Path, dict]:
    project = Path(project).resolve()
    clean = root(project) / "inputs/protein_clean.pdb"
    choice = _json(root(project) / "inputs/protein_choice.json")
    if not choice or not clean.is_file() or _sha_file(clean) != choice.get("clean_sha256"):
        raise RuntimeError("Step 2 reviewed protein was changed or not saved; repeat protein review.")
    return clean, choice


def _atoms(pdb: Path) -> list[dict]:
    sites = []
    seen = set()
    for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("HETATM") or line[17:20].strip().upper() != "NI":
            continue
        element = (line[76:78].strip().upper() or "NI")
        if element != "NI":
            raise ValueError(f"Residue NI with unexpected element {element}; manual review required")
        chain = line[21:22].strip() or "-"
        residue_number = line[22:27].strip()
        key = f"{chain}:{residue_number}"
        if key in seen:
            raise ValueError(f"Multiple NI atoms in {key}; manual review required")
        seen.add(key)
        try:
            xyz = [float(line[n:n+8]) for n in (30, 38, 46)]
        except ValueError as exc:
            raise ValueError(f"Invalid Ni coordinates for {key}") from exc
        if not all(math.isfinite(v) for v in xyz):
            raise ValueError(f"Nonfinite Ni coordinates at {key}")
        sites.append({
            "key": key, "chain": chain, "residue_number": residue_number,
            "atom_name": line[12:16].strip(), "xyz": xyz,
            "pdb_line": line,
        })
    return sites


def sites(project: Path) -> list[dict]:
    clean, _ = _reviewed_source(project)
    return [
        {"key": site["key"], "xyz": site["xyz"], "atom_name": site["atom_name"]}
        for site in _atoms(clean)
    ]


def prepare_intermediate(project: Path) -> dict:
    """Run Meeko without Ni; do NOT install an incomplete production receptor."""
    from .core import prepare_receptor
    project = Path(project).resolve()
    clean, choice = _reviewed_source(project)
    metals = _atoms(clean)
    if not metals:
        raise RuntimeError("Reviewed protein has no Ni ions. Use ordinary receptor preparation.")

    folder = root(project) / "inputs"
    temporary = folder / "protein_for_meeko_without_Ni.pdb"
    intermediate = folder / "protein_meeko_Ni_intermediate.pdbqt"
    log = root(project) / "logs/receptor_meeko_Ni_intermediate.log"
    metadata = folder / "Ni_intermediate.json"

    # Require an exact current-clean-hash match before reusing previous work.
    old = _json(metadata)
    if (intermediate.is_file() and intermediate.stat().st_size
            and old.get("reviewed_clean_sha256") == _sha_file(clean)
            and old.get("intermediate_sha256") == _sha_file(intermediate)
            and old.get("excluded_Ni_sites") == [x["key"] for x in metals]):
        return {"status": "REUSED", "intermediate": str(intermediate),
                "excluded_Ni_sites": old["excluded_Ni_sites"], "log": str(log),
                "production_receptor_installed": False}

    excluded = {site["pdb_line"] for site in metals}
    lines = clean.read_text(encoding="utf-8", errors="replace").splitlines()
    temporary.write_text("\n".join(line for line in lines if line not in excluded)+"\n",
                         encoding="utf-8")
    if _atoms(temporary):
        raise RuntimeError("Temporary Meeko input still contains Ni. Preparation aborted.")
    prepared = prepare_receptor(temporary, intermediate, log, allow_bad_res=False)
    if any(x.upper() == "NI" for x in pdbqt_atom_types(prepared)):
        raise RuntimeError("Unexpected Ni in Meeko intermediate; aborted.")
    record = {
        "reviewed_clean_sha256": _sha_file(clean),
        "temporary_input_sha256": _sha_file(temporary),
        "intermediate_sha256": _sha_file(prepared),
        "excluded_Ni_sites": [x["key"] for x in metals],
        "Ni_coordinates_A": {x["key"]: x["xyz"] for x in metals},
        "method": "Meeko preparation with Ni omitted only from temporary input",
        "status": "UNVALIDATED_INTERMEDIATE_NO_METAL",
        "meeko_log": str(log),
    }
    _save_json(metadata, record)
    _event(project, "NICKEL_INTERMEDIATE_PREPARED", record)
    return {"status": "PREPARED", "intermediate": str(prepared),
            "excluded_Ni_sites": record["excluded_Ni_sites"], "log": str(log),
            "production_receptor_installed": False}


def parse_charges(raw: str) -> dict[str, float]:
    """Exact format C:798=1.0,C:799=1.0; no implicit guesses or default."""
    parsed = {}
    for part in re.split(r"[;,]", raw):
        if not part.strip():
            continue
        match = re.fullmatch(r"\s*([A-Za-z0-9-]+):([0-9]+[A-Za-z]?)\s*=\s*"
                             r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*", part)
        if not match:
            raise ValueError("Metal charges must use C:798=<charge>,C:799=<charge> (e units)")
        key = f"{match.group(1)}:{match.group(2)}"
        value = float(match.group(3))
        if key in parsed or not math.isfinite(value) or abs(value) > 4:
            raise ValueError(f"Duplicate or invalid Ni partial charge for {key}")
        parsed[key] = value
    if not parsed:
        raise ValueError("All site-specific Ni partial charges are required")
    return parsed


def _atom_parameter_types(path: Path) -> set[str]:
    found = set()
    for line in path.read_text(encoding="utf-8", errors="strict").splitlines():
        columns = line.split("#", 1)[0].split()
        if columns and columns[0] == "atom_par":
            if len(columns) != 12:
                raise ValueError(f"Malformed atom_par entry in {path.name}: {line[:100]}")
            found.add(columns[1])
    return found


def _metal_pdbqt_line(original: str, serial: int, charge: float, atom_type: str) -> str:
    if len(atom_type) > 2 or not atom_type.isascii() or not atom_type.isalnum():
        raise ValueError("AutoDock4 metal atom type must be 1–2 ASCII letters/digits")
    if serial > 99999:
        raise ValueError("PDBQT atom serial overflow")
    line = "HETATM" + f"{serial:5d}" + original[11:66]
    # Fixed PDBQT columns: partial charge 70:76; AutoDock atom type 77:79.
    line = line.ljust(70) + f"{charge:+6.3f}" + " " + f"{atom_type:>2s}"
    return line + "\n"


def reconstruct(project: Path, *, atom_type: str, charges: dict[str, float],
                parameter_source: str, charge_source: str,
                nonbonded_model_acknowledged: bool = False) -> dict:
    project = Path(project).resolve()
    clean, _ = _reviewed_source(project)
    metals = _atoms(clean)
    if not metals:
        raise RuntimeError("No reviewed nickel ions to reconstruct")
    if not nonbonded_model_acknowledged:
        raise RuntimeError("Explicitly acknowledge the limitations of fixed-centre nonbonded Ni treatment")
    if not parameter_source.strip() or not charge_source.strip():
        raise ValueError("Record scientific sources for both AutoDock parameters and Ni charges")

    intermediate = root(project) / "inputs/protein_meeko_Ni_intermediate.pdbqt"
    metadata = _json(root(project) / "inputs/Ni_intermediate.json")
    if (not intermediate.is_file()
            or metadata.get("reviewed_clean_sha256") != _sha_file(clean)
            or metadata.get("intermediate_sha256") != _sha_file(intermediate)
            or metadata.get("excluded_Ni_sites") != [x["key"] for x in metals]):
        raise RuntimeError("No current Ni-free Meeko intermediate. Prepare it in Step 2 first.")
    if set(charges) != {x["key"] for x in metals}:
        raise ValueError("Specify exactly one charge for each crystallographic Ni site: "
                         + ", ".join(x["key"] for x in metals))
    for site in metals:
        charge = charges[site["key"]]
        if not isinstance(charge, (float, int)) or not math.isfinite(charge) or abs(charge) > 4:
            raise ValueError(f"Invalid explicit charge for {site['key']}")

    param_file = parameter_file_path(project)
    if param_file is None:
        raise RuntimeError(
            "Complete AutoDock4 parameter library with Ni atom_par record is required; "
            "upload the scientifically selected library before reconstruction."
        )
    declared = _atom_parameter_types(param_file)
    present = set(pdbqt_atom_types(intermediate))
    missing = present | {atom_type}
    missing -= declared
    if missing:
        raise RuntimeError("Custom AD4 parameter library lacks atom types: " + ", ".join(sorted(missing)))
    if not atom_type or atom_type.upper() not in {"NI"}:
        raise ValueError("Nickel-specific AD4 atom type Ni or NI is required")

    intermediate_text = intermediate.read_text(encoding="utf-8", errors="replace")
    atom_lines = [line for line in intermediate_text.splitlines()
                  if line.startswith(("ATOM", "HETATM"))]
    if not atom_lines:
        raise RuntimeError("Meeko intermediate PDBQT has no atoms")
    if any(line[17:20].strip() == "NI" for line in atom_lines):
        raise RuntimeError("Intermediate already contains nickel; refusing to duplicate ions")
    serials = [int(line[6:11]) for line in atom_lines]
    result_lines = [
        _metal_pdbqt_line(site["pdb_line"], max(serials)+index+1,
                          float(charges[site["key"]]), atom_type)
        for index, site in enumerate(metals)
    ]
    final = root(project) / "inputs/protein_Ni_reconstructed.pdbqt"
    # Maintain Meeko protein atoms exactly; append Ni at rigid receptor tail.
    # Keep the metal records inside the rigid PDBQT body, before any END
    # marker. Meeko versions differ in whether they include terminal records.
    base_lines = intermediate_text.splitlines()
    end_marker = next((i for i, line in enumerate(base_lines)
                       if line.strip() in {"END", "ENDMDL"}), len(base_lines))
    base_lines[end_marker:end_marker] = [line.rstrip("\n") for line in result_lines]
    final.write_text("\n".join(base_lines) + "\n", encoding="utf-8")
    output_metal = [x for x in final.read_text().splitlines()
                    if x.startswith(("ATOM", "HETATM")) and x[17:20].strip() == "NI"]
    if len(output_metal) != len(metals):
        final.unlink(missing_ok=True)
        raise RuntimeError("Reconstruction atom count mismatch")
    for source_site, output in zip(metals, output_metal):
        observed = [float(output[i:i+8]) for i in (30, 38, 46)]
        if any(abs(a-b) > 0.001 for a,b in zip(observed,source_site["xyz"])):
            raise RuntimeError("Reconstructed nickel coordinates shifted")
        if output[77:79].strip() != atom_type:
            raise RuntimeError("Reconstructed nickel atom type mismatch")
    metal_separation = None
    if len(metals) == 2:
        metal_separation = math.dist(metals[0]["xyz"], metals[1]["xyz"])
    record = {
        "model": "fixed_centres_nonbonded_custom_AD4_unvalidated",
        "scientific_review_required": True,
        "native_redocking": "NOT_PERFORMED",
        "NI_count": len(metals),
        "NI_sites": [{"key": s["key"],"xyz_A": s["xyz"],
                      "partial_charge_e": charges[s["key"]],"AD4_type":atom_type}
                     for s in metals],
        "Ni_Ni_distance_A": metal_separation,
        "intermediate_sha256": _sha_file(intermediate),
        "reviewed_clean_sha256": _sha_file(clean),
        "AD4_parameters_sha256": _sha_file(param_file),
        "parameter_source": parameter_source.strip(),
        "charge_source": charge_source.strip(),
        "reconstructed_sha256": _sha_file(final),
        "limitations": [
            "Unvalidated fixed-centre nonbonded metal model; no explicit Ni-ligand coordination bonds.",
            "Meeko assigned protein charges without explicit metal electronic polarization.",
            "The nickel and KCX metal environment require experimental redocking/geometry review.",
        ],
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    _save_json(root(project) / "inputs/Ni_reconstruction_audit.json", record)
    install_reviewed_receptor(project, final)
    # Record the actual accepted audit fingerprint alongside the receptor.
    review = _json(root(project) / "inputs/receptor_review.json")
    review.update({
        "metal_model": record["model"], "metal_audit_sha256":
        _sha_file(root(project) / "inputs/Ni_reconstruction_audit.json"),
    })
    _save_json(root(project) / "inputs/receptor_review.json", review)
    _event(project, "NICKEL_RECONSTRUCTED", record)
    return {
        "status": "RECONSTRUCTED_NOT_VALIDATED",
        "receptor": str(root(project) / "receptor.pdbqt"),
        "reconstructed_sites": len(metals),
        "Ni_Ni_distance_A": metal_separation,
        "audit": str(root(project) / "inputs/Ni_reconstruction_audit.json"),
        "requires": "Human parameter review and crystallographic HAE redocking",
    }
