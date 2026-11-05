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
import tempfile
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
        if output.stat().st_size == 0 and not existing:
            # A prior Open Babel failure may leave a zero-byte candidate. It is
            # not a scientific artifact and may be replaced by a new attempt.
            output.unlink()
        else:
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
    gasteiger_log = log.with_name(".receptor_openbabel_gasteiger.log")
    exit_code = run_command(cmd, output.parent, gasteiger_log)
    commands = [cmd]
    charge_assignment = "GASTEIGER"
    # Open Babel 3.1.0 can return success while its Gasteiger charge
    # perception rejects this metal-containing protein and writes 0 bytes.
    # Retry with Open Babel's neutral-charge PDBQT export. This is still an
    # inspectable candidate only; it is not a validated AutoDock receptor.
    if exit_code != 0 or not output.is_file() or output.stat().st_size == 0:
        fallback_cmd = [
            obabel, str(src.resolve()), "-O", str(output.resolve()),
            "-xr", "-xc", "-xn", "-xp", "-xh",
            "--partialcharge", "none",
        ]
        fallback_log = log.with_name(".receptor_openbabel_neutral.log")
        exit_code = run_command(fallback_cmd, output.parent, fallback_log)
        commands.append(fallback_cmd)
        charge_assignment = "NONE_FALLBACK_AFTER_GASTEIGER_EMPTY_OUTPUT"
        gasteiger_text = gasteiger_log.read_text(encoding="utf-8", errors="replace") if gasteiger_log.is_file() else ""
        fallback_text = fallback_log.read_text(encoding="utf-8", errors="replace") if fallback_log.is_file() else ""
        log.write_text(
            "Open Babel Gasteiger attempt:\n" + gasteiger_text
            + "\nOpen Babel neutral-charge fallback:\n" + fallback_text,
            encoding="utf-8",
        )
        for temporary_log in (gasteiger_log, fallback_log):
            if temporary_log.exists():
                temporary_log.unlink()
    else:
        log.write_text(
            gasteiger_log.read_text(encoding="utf-8", errors="replace"),
            encoding="utf-8",
        )
        if gasteiger_log.exists():
            gasteiger_log.unlink()
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
        "command": commands[-1],
        "commands_attempted": commands,
        "charge_assignment": charge_assignment,
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


_METAL_ELEMENTS = {
    "LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "AL",
    "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG", "AG", "AU",
}


def _ad4_atom_types(path: Path) -> set[str]:
    types = set()
    for line in path.read_text(encoding="utf-8", errors="strict").splitlines():
        parts = line.split()
        if parts and parts[0] == "atom_par" and len(parts) >= 12:
            types.add(parts[1])
    return types


def _metal_line(line: str, serial: int) -> str:
    """Copy a structural metal record, changing only serial and charge.

    The zero charge is a deliberate unassigned sentinel. It is not a claimed
    Ni/metal charge and therefore cannot be mistaken for a parameterization.
    """
    updated = line[:6] + f"{serial:5d}" + line[11:]
    if len(updated) < 79:
        updated = updated.ljust(79)
    return updated[:70] + f"{0.0:6.3f}" + updated[76:]


def generate_protonated_charge_candidate(project: Path) -> dict:
    """Create an audited, non-approved protonated/charged receptor candidate.

    Open Babel is applied only to a temporary copy with crystallographic metals
    removed. The existing structural candidate supplies the metal records,
    which are restored with zero/unassigned charges and never parameterized by
    this function. No original PDB, structural candidate, or docking artifact
    is overwritten.
    """
    from .guided_workflow import _json, _save_json, _event
    from .manual_docking import root
    from .core import run_command

    project = Path(project).resolve()
    home = root(project)
    src = home / "inputs/protein_clean.pdb"
    structural = home / "inputs/receptor_openbabel_candidate.pdbqt"
    structural_audit_path = home / "inputs/receptor_openbabel_audit.json"
    if not src.is_file() or not structural.is_file() or not structural_audit_path.is_file():
        raise RuntimeError("Generate and audit the existing Open Babel structural candidate first")
    structural_audit = _json(structural_audit_path)
    if structural_audit.get("candidate_sha256") != _sha(structural):
        raise RuntimeError("Structural Open Babel candidate hash does not match its audit")
    if structural_audit.get("missing_input_heavy_atoms") or structural_audit.get("metal_atoms_missing"):
        raise RuntimeError("Structural candidate is incomplete; protonation/charge workflow is blocked")
    obabel = shutil.which("obabel")
    if not obabel:
        raise RuntimeError("Open Babel executable 'obabel' not found on Ubuntu PATH.")
    source_atoms = _atoms(src, input_pdb=True)
    structural_atoms = _atoms(structural, input_pdb=False)
    metals = [a for a in source_atoms if a["element"].upper() in _METAL_ELEMENTS]
    metal_keys = {(a["identity"], a["xyz"]) for a in metals}
    structural_metal_lines = []
    for line in structural.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 79:
            continue
        atom = next((a for a in structural_atoms
                     if a["identity"] == (line[21:22].strip() or "-", line[17:20].strip(),
                                           line[22:27].strip(), line[12:16].strip())
                     and a["xyz"] == tuple(float(line[i:i + 8]) for i in (30, 38, 46))), None)
        # Atom type alone is insufficient: AD4 uses NA for some non-metal
        # nitrogen atoms. Metal identity is therefore anchored to the reviewed
        # PDB coordinates and residue/atom identity.
        if atom and (atom["identity"], atom["xyz"]) in metal_keys:
            structural_metal_lines.append(line)
    if len(structural_metal_lines) != len(metals):
        raise RuntimeError("Existing structural candidate does not retain every crystallographic metal")

    library = home / "custom_parameters.dat"
    if not library.is_file():
        raise RuntimeError("Select a project AutoDock4 parameter library in Step 2 before preparing receptor charges")
    library_types = _ad4_atom_types(library)
    output = home / "inputs/receptor_openbabel_protonated_charge_candidate.pdbqt"
    audit_path = home / "inputs/receptor_openbabel_protonated_charge_audit.json"
    log = home / "logs/receptor_openbabel_protonated_charge.log"
    existing = _json(audit_path)
    if output.is_file() and existing.get("source_structural_candidate_sha256") == _sha(structural) \
            and existing.get("parameter_library_sha256") == _sha(library) \
            and existing.get("candidate_sha256") == _sha(output):
        return {**existing, "status": "REUSED", "candidate": str(output), "audit": str(audit_path)}
    if output.exists():
        raise RuntimeError(f"Existing protonated/charged candidate differs; preserved it: {output}")

    with tempfile.TemporaryDirectory(prefix="sidock-openbabel-charge-") as temp_name:
        temp = Path(temp_name)
        nonmetal_pdb = temp / "protein_nonmetal.pdb"
        probe = temp / "protein_nonmetal_protonated.pdbqt"
        kept = []
        for line in src.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith(("ATOM  ", "HETATM")) and line[76:78].strip().upper() in _METAL_ELEMENTS:
                continue
            kept.append(line)
        nonmetal_pdb.write_text("\n".join(kept) + "\n", encoding="utf-8")
        cmd = [obabel, str(nonmetal_pdb), "-O", str(probe), "-h",
               "-xr", "-xc", "-xn", "-xp", "-xh", "--partialcharge", "gasteiger"]
        probe_log = temp / "openbabel.log"
        exit_code = run_command(cmd, temp, probe_log)
        probe_text = probe_log.read_text(encoding="utf-8", errors="replace") if probe_log.is_file() else ""
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text("Open Babel non-metal protonation/Gasteiger probe:\n" + probe_text, encoding="utf-8")
        if exit_code != 0 or not probe.is_file() or probe.stat().st_size == 0:
            raise RuntimeError(f"Open Babel protonation/charge probe failed (exit={exit_code}); see {log}")
        probe_atoms = _atoms(probe, input_pdb=False)
        if not probe_atoms:
            raise RuntimeError("Open Babel protonation/charge probe produced no atoms")
        probe_lines = probe.read_text(encoding="utf-8", errors="replace").splitlines()

        input_nonmetal_heavy = [a for a in source_atoms
                                if a["element"] not in {"H", "D"}
                                and a["element"].upper() not in _METAL_ELEMENTS]
        def snap(atom: dict) -> tuple:
            return (*atom["identity"], *(round(v, 3) for v in atom["xyz"]))
        probe_snapshots = collections.Counter(snap(a) for a in probe_atoms)
        missing = [a for a in input_nonmetal_heavy if probe_snapshots[snap(a)] == 0]
        for atom in input_nonmetal_heavy:
            probe_snapshots[snap(atom)] -= 1
        if missing:
            raise RuntimeError(f"Open Babel protonation probe lost {len(missing)} non-metal heavy atoms")
        coordinate_deltas = []
        for expected in input_nonmetal_heavy:
            matches = [a for a in probe_atoms if a["identity"] == expected["identity"]]
            if matches:
                coordinate_deltas.append(max(abs(x - y) for x, y in zip(expected["xyz"], matches[0]["xyz"])))

        max_serial = max((int(line[6:11]) for line in probe_lines
                          if line.startswith(("ATOM  ", "HETATM")) and line[6:11].strip().isdigit()), default=0)
        metal_lines = [_metal_line(line, max_serial + i + 1)
                       for i, line in enumerate(structural_metal_lines)]
        final_lines = []
        inserted = False
        for line in probe_lines:
            if line in {"TER", "END", "ENDMDL"} and not inserted:
                final_lines.extend(metal_lines)
                inserted = True
            if line not in {"END", "ENDMDL"}:
                final_lines.append(line)
        if not inserted:
            final_lines.extend(metal_lines)
        final_lines.append("END")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("\n".join(final_lines) + "\n", encoding="utf-8")

    final_atoms = _atoms(output, input_pdb=False)
    all_types = sorted({a["atom_type"] for a in final_atoms})
    missing_types = sorted(set(all_types) - library_types)
    if missing_types:
        output.unlink()
        raise RuntimeError("AD4 library lacks generated atom types: " + ", ".join(missing_types))
    final_metal_keys = {(a["identity"], a["xyz"]) for a in metals}
    final_metals = [a for a in final_atoms if (a["identity"], a["xyz"]) in final_metal_keys]
    nonmetal_charges = [a["charge"] for a in final_atoms
                        if (a["identity"], a["xyz"]) not in final_metal_keys]
    polar_h = [a for a in final_atoms if a["atom_type"].upper() in {"HD", "HS"}]
    audit = {
        "status": "CANDIDATE_REQUIRES_REVIEW",
        "engine": "OPEN_BABEL",
        "command": cmd,
        "source_clean_pdb_sha256": _sha(src),
        "source_structural_candidate": str(structural),
        "source_structural_candidate_sha256": _sha(structural),
        "candidate_sha256": _sha(output),
        "parameter_library": str(library),
        "parameter_library_sha256": _sha(library),
        "input_heavy_atom_count": len([a for a in source_atoms if a["element"] not in {"H", "D"}]),
        "final_atom_count": len(final_atoms),
        "final_heavy_atom_count": len([a for a in final_atoms if a["atom_type"].upper() not in {"H", "HD", "HS"}]),
        "added_hydrogen_count": len(final_atoms) - len(source_atoms),
        "polar_hydrogen_count": len(polar_h),
        "nonzero_nonmetal_charge_count": sum(abs(q) > 1e-8 for q in nonmetal_charges),
        "nonmetal_charge_sum_e": round(sum(nonmetal_charges), 4),
        "nonmetal_heavy_coordinate_max_abs_delta_angstrom": round(max(coordinate_deltas, default=0.0), 6),
        "atom_types": all_types,
        "ad4_atom_types_missing": missing_types,
        "metal_count": len(final_metals),
        "metal_elements": sorted({a["atom_type"] for a in final_metals}),
        "metal_charge_status": "UNASSIGNED_NOT_INVENTED",
        "metal_charges_e": [a["charge"] for a in final_metals],
        "metal_parameterization": "NOT_VALIDATED",
        "protonation_method": "OPEN_BABEL_ADD_H_GASTEIGER_NONMETAL_PROBE",
        "protonation_review": "NOT_PERFORMED",
        "scientific_validation": "NOT_PERFORMED",
        "approved_for_docking": False,
        "limitations": [
            "Gasteiger charges are a heuristic diagnostic, not a validated protein force field.",
            "No pH-specific protonation, tautomer, or histidine-state review was performed.",
            "Metal polarization, coordination bonding, and metal charge assignment were not parameterized.",
            "AD4 atom-type coverage does not constitute scientific approval for docking.",
        ],
        "log": str(log),
        "created": datetime.now(timezone.utc).isoformat(),
    }
    _save_json(audit_path, audit)
    _event(project, "RECEPTOR_OPENBABEL_PROTONATED_CHARGE_CANDIDATE_GENERATED", {
        "candidate_sha256": audit["candidate_sha256"], "status": audit["status"],
        "metal_charge_status": audit["metal_charge_status"],
    })
    return {**audit, "candidate": str(output), "audit": str(audit_path)}


def audit_and_correct_hydrogen_types(project: Path) -> dict:
    """Audit hydrogen connectivity and correct only demonstrable H typing.

    The input candidate is copied line-for-line except for the PDBQT atom type
    of hydrogens. Bonding is assigned from same-residue coordinate proximity,
    with element-specific covalent-distance windows. No coordinates, charges,
    heavy atoms, or metal records are changed.
    """
    from .guided_workflow import _json, _save_json, _event
    from .manual_docking import root

    project = Path(project).resolve()
    home = root(project)
    source = home / "inputs/receptor_openbabel_protonated_charge_candidate.pdbqt"
    source_audit_path = home / "inputs/receptor_openbabel_protonated_charge_audit.json"
    if not source.is_file() or not source_audit_path.is_file():
        raise RuntimeError("Generate the existing Open Babel protonated/charged candidate first")
    source_audit = _json(source_audit_path)
    if source_audit.get("candidate_sha256") != _sha(source):
        raise RuntimeError("Protonated/charged candidate hash does not match its audit")
    output = home / "inputs/receptor_openbabel_hydrogen_typed_candidate.pdbqt"
    audit_path = home / "inputs/receptor_openbabel_hydrogen_typed_audit.json"
    if output.is_file():
        existing = _json(audit_path)
        if existing.get("source_candidate_sha256") == _sha(source) and existing.get("candidate_sha256") == _sha(output):
            return {**existing, "status": "REUSED", "candidate": str(output), "audit": str(audit_path)}
        raise RuntimeError(f"Existing hydrogen-typed candidate differs; preserved it: {output}")

    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    records = []
    for index, line in enumerate(lines):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        records.append({
            "index": index, "line": line, "identity": (line[21:22].strip() or "-", line[17:20].strip(),
            line[22:27].strip(), line[12:16].strip()),
            "xyz": tuple(float(line[i:i + 8]) for i in (30, 38, 46)),
            "type": line[77:79].strip(), "charge": float(line[70:76]),
        })
    hydrogen = [r for r in records if r["line"][12:16].strip().startswith("H")]
    heavy = [r for r in records if r not in hydrogen]
    polar_types = {"N", "NA", "OA", "O", "S", "SA", "P", "PA"}
    connectivity = []
    corrections = []
    unresolved = []
    for h in hydrogen:
        same_residue = [a for a in heavy if a["identity"][:3] == h["identity"][:3]]
        candidates = sorted((math.dist(h["xyz"], a["xyz"]), a) for a in same_residue)
        if not candidates:
            unresolved.append({"hydrogen": h["identity"], "reason": "no_same_residue_heavy_atom"})
            continue
        distance, bonded = candidates[0]
        # These windows deliberately include crystallographic S-H distances
        # while rejecting nonbonded contacts in ordinary protein geometry.
        element_type = bonded["type"].upper()
        max_distance = {"C": 1.25, "A": 1.25, "N": 1.30, "NA": 1.30,
                        "OA": 1.25, "O": 1.25, "S": 1.60, "SA": 1.60,
                        "P": 1.65, "PA": 1.65}.get(element_type, 1.45)
        if distance > max_distance:
            unresolved.append({"hydrogen": h["identity"], "nearest_heavy": bonded["identity"],
                               "nearest_type": bonded["type"], "distance_angstrom": round(distance, 4),
                               "reason": "outside_element_specific_bond_window"})
            continue
        new_type = "HD" if element_type in polar_types else "H"
        connectivity.append({
            "hydrogen": h["identity"], "hydrogen_serial": int(h["line"][6:11]),
            "hydrogen_type_before": h["type"], "hydrogen_type_after": new_type,
            "bonded_heavy_atom": bonded["identity"], "bonded_heavy_type": bonded["type"],
            "distance_angstrom": round(distance, 4), "hydrogen_charge_e": h["charge"],
        })
        if h["type"] != new_type:
            updated = h["line"]
            if len(updated) < 79:
                updated = updated.ljust(79)
            lines[h["index"]] = updated[:77] + f"{new_type:>2}" + updated[79:]
            corrections.append(connectivity[-1])

    if unresolved:
        raise RuntimeError(f"Could not validate connectivity for {len(unresolved)} hydrogens")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    corrected_records = _atoms(output, input_pdb=False)
    original_records = _atoms(source, input_pdb=False)
    heavy_key = lambda a: (a["identity"], a["xyz"], a["charge"])
    original_heavy = [a for a in original_records if not a["identity"][3].startswith("H")]
    corrected_heavy = [a for a in corrected_records if not a["identity"][3].startswith("H")]
    heavy_preserved = collections.Counter(heavy_key(a) for a in original_heavy) == collections.Counter(heavy_key(a) for a in corrected_heavy)
    if not heavy_preserved:
        output.unlink()
        raise RuntimeError("Hydrogen correction altered a heavy atom coordinate or charge")

    residue_charge = collections.defaultdict(float)
    for atom in corrected_records:
        if atom["atom_type"] != "Ni":
            residue_charge[atom["identity"][:3]] += atom["charge"]
    audit = {
        "status": "CANDIDATE_REQUIRES_REVIEW",
        "source_candidate": str(source),
        "source_candidate_sha256": _sha(source),
        "candidate_sha256": _sha(output),
        "original_hydrogen_count": len(hydrogen),
        "original_hydrogen_type_distribution": dict(collections.Counter(r["type"] for r in hydrogen)),
        "corrected_hydrogen_type_distribution": dict(collections.Counter(r["hydrogen_type_after"] for r in connectivity)),
        "hydrogen_bonded_heavy_type_distribution": dict(collections.Counter(r["bonded_heavy_type"] for r in connectivity)),
        "polar_hydrogen_count": sum(r["hydrogen_type_after"] == "HD" for r in connectivity),
        "carbon_bound_hydrogen_count": sum(r["hydrogen_type_after"] == "H" for r in connectivity),
        "hydrogen_charge_nonzero_count": sum(abs(r["hydrogen_charge_e"]) > 1e-8 for r in connectivity),
        "heavy_atoms_preserved_exactly": heavy_preserved,
        "hydrogen_type_correction": "CARBON_BOUND_H_TO_H; HETEROATOM_BOUND_H_TO_HD",
        "charge_model": "UNCHANGED_OPEN_BABEL_GASTEIGER_NONMETAL_PROBE",
        "nonmetal_charge_sum_e": round(sum(a["charge"] for a in corrected_records if a["atom_type"] != "Ni"), 4),
        "residues_with_abs_charge_at_least_0_05_e": [
            {"chain": k[0], "resname": k[1], "resid": k[2], "charge_e": round(v, 4)}
            for k, v in sorted(residue_charge.items()) if abs(v) >= 0.05
        ],
        "metal_charge_status": "UNASSIGNED_NOT_INVENTED",
        "approved_for_docking": False,
        "scientific_validation": "NOT_PERFORMED",
        "limitations": [
            "Connectivity was assigned from coordinates and residue-local proximity; no explicit bond graph was present in PDBQT.",
            "Hydrogen coordinates and Gasteiger heavy-atom charges came from Open Babel and were not pH-validated.",
            "All generated hydrogen charges are zero in the source candidate and were not invented by this correction.",
            "Residue charge sums are diagnostic distributions, not validated protonation-state assignments.",
            "Metal coordination and metal charges remain unparameterized.",
        ],
        "created": datetime.now(timezone.utc).isoformat(),
    }
    _save_json(audit_path, audit)
    _event(project, "RECEPTOR_OPENBABEL_HYDROGEN_TYPES_AUDITED", {
        "candidate_sha256": audit["candidate_sha256"], "corrections": len(corrections),
        "polar_hydrogens": audit["polar_hydrogen_count"],
    })
    return {**audit, "candidate": str(output), "audit": str(audit_path)}


def prepare_exploratory_receptor(project: Path) -> dict:
    """Strip only nonpolar H for a separate, explicitly exploratory AD4 copy."""
    from .guided_workflow import _json, _save_json, _event
    from .manual_docking import root

    project = Path(project).resolve()
    home = root(project)
    source = home / "inputs/receptor_openbabel_hydrogen_typed_candidate.pdbqt"
    source_audit_path = home / "inputs/receptor_openbabel_hydrogen_typed_audit.json"
    parameter_file = home / "custom_parameters.dat"
    if not source.is_file() or not source_audit_path.is_file():
        raise RuntimeError("Generate the corrected hydrogen-typed receptor first")
    if not parameter_file.is_file():
        raise RuntimeError("Exploratory receptor requires manual_docking/custom_parameters.dat")
    source_audit = _json(source_audit_path)
    if source_audit.get("candidate_sha256") != _sha(source):
        raise RuntimeError("Hydrogen-typed candidate hash does not match its audit")
    output = home / "inputs/receptor_openbabel_exploratory.pdbqt"
    audit_path = home / "inputs/receptor_openbabel_exploratory_audit.json"
    if output.is_file():
        existing = _json(audit_path)
        if (existing.get("source_candidate_sha256") == _sha(source)
                and existing.get("parameter_file_sha256") == _sha(parameter_file)
                and existing.get("candidate_sha256") == _sha(output)):
            return {**existing, "status": "REUSED", "candidate": str(output), "audit": str(audit_path)}
        raise RuntimeError(f"Existing exploratory receptor differs; preserved it: {output}")

    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    records = []
    for index, line in enumerate(lines):
        if line.startswith(("ATOM  ", "HETATM")):
            records.append({"index": index, "line": line,
                            "identity": (line[21:22].strip() or "-", line[17:20].strip(),
                                         line[22:27].strip(), line[12:16].strip()),
                            "xyz": tuple(float(line[i:i + 8]) for i in (30, 38, 46)),
                            "type": line[77:79].strip(), "charge": float(line[70:76])})
    hydrogens = [r for r in records if r["identity"][3].startswith("H")]
    heavy = [r for r in records if r not in hydrogens]
    remove = set()
    merge_targets = collections.defaultdict(list)
    unresolved = []
    for h in hydrogens:
        if h["type"] != "H":
            continue
        same_residue = [a for a in heavy if a["identity"][:3] == h["identity"][:3]]
        candidates = sorted((math.dist(h["xyz"], a["xyz"]), a) for a in same_residue)
        if not candidates:
            unresolved.append(h["identity"])
            continue
        distance, bonded = candidates[0]
        if bonded["type"].upper() not in {"C", "A"} or distance > 1.25:
            unresolved.append({"hydrogen": h["identity"], "nearest": bonded["identity"],
                               "nearest_type": bonded["type"], "distance_angstrom": round(distance, 4)})
            continue
        remove.add(h["index"])
        merge_targets[bonded["index"]].append(h["charge"])
    expected_nonpolar = source_audit.get("carbon_bound_hydrogen_count")
    expected_polar = source_audit.get("polar_hydrogen_count")
    if not isinstance(expected_nonpolar, int) or not isinstance(expected_polar, int):
        raise RuntimeError("Hydrogen-type audit lacks the required expected counts")
    if unresolved or len(remove) != expected_nonpolar:
        raise RuntimeError(
            f"Expected {expected_nonpolar} audited nonpolar H; "
            f"found {len(remove)}; unresolved={len(unresolved)}"
        )
    for target_index, charges in merge_targets.items():
        line = lines[target_index]
        current = float(line[70:76])
        updated = line[:70] + f"{current + sum(charges):6.3f}" + line[76:]
        lines[target_index] = updated
    retained_lines = [line for i, line in enumerate(lines) if i not in remove]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(retained_lines) + "\n", encoding="utf-8")
    final_atoms = _atoms(output, input_pdb=False)
    final_heavy = [a for a in final_atoms if not a["identity"][3].startswith("H")]
    final_polar_h = [a for a in final_atoms if a["identity"][3].startswith("H")]
    original_heavy = [a for a in _atoms(source, input_pdb=False) if not a["identity"][3].startswith("H")]
    heavy_preserved = collections.Counter((a["identity"], a["xyz"]) for a in original_heavy) == collections.Counter((a["identity"], a["xyz"]) for a in final_heavy)
    types = sorted({a["atom_type"] for a in final_atoms})
    parameter_types = _ad4_atom_types(parameter_file)
    missing_types = sorted(set(types) - parameter_types)
    if not heavy_preserved or len(final_polar_h) != expected_polar or missing_types:
        output.unlink()
        raise RuntimeError(
            f"Exploratory receptor validation failed: heavy={heavy_preserved}, "
            f"polar_H={len(final_polar_h)}/{expected_polar}, missing_types={missing_types}"
        )
    ni = [a for a in final_atoms if a["atom_type"] == "Ni"]
    audit = {
        "status": "TECHNICAL_EXPLORATORY_CANDIDATE",
        "source_candidate": str(source), "source_candidate_sha256": _sha(source),
        "candidate": str(output), "candidate_sha256": _sha(output),
        "parameter_file": str(parameter_file), "parameter_file_sha256": _sha(parameter_file),
        "removed_carbon_bound_nonpolar_hydrogens": len(remove),
        "retained_heavy_atom_count": len(final_heavy), "retained_polar_hydrogen_count": len(final_polar_h),
        "heavy_atoms_preserved_exactly": heavy_preserved, "atom_types": types,
        "missing_parameter_types": missing_types,
        "ni_count": len(ni), "ni_coordinates": [a["xyz"] for a in ni],
        "ni_charges_e": [a["charge"] for a in ni],
        "ni_charge_status": "ZERO_UNASSIGNED_NOT_INVENTED",
        "scientific_approval": False, "technical_exploratory_mode": True,
        "limitations": ["Provisional Open Babel charges retained.", "Ni charges and coordination are not scientifically validated.",
                        "This artifact is eligible only for a hash-bound technical smoke test."],
        "created": datetime.now(timezone.utc).isoformat(),
    }
    _save_json(audit_path, audit)
    _event(project, "RECEPTOR_OPENBABEL_EXPLORATORY_CANDIDATE_GENERATED", {
        "candidate_sha256": audit["candidate_sha256"], "parameter_file_sha256": audit["parameter_file_sha256"]})
    return {**audit, "candidate": str(output), "audit": str(audit_path)}


def prepare_openbabel_receptor(project: Path) -> dict:
    """Complete the reviewed Open Babel chain; no smoke run or approval."""
    from .manual_docking import root
    if not (root(project) / "custom_parameters.dat").is_file():
        raise RuntimeError("Select the AutoDock4 parameter library in Step 2 first")
    stages = {}
    for name, prepare in (
        ("structural", generate_candidate),
        ("charges", generate_protonated_charge_candidate),
        ("hydrogen_types", audit_and_correct_hydrogen_types),
        ("exploratory", prepare_exploratory_receptor),
    ):
        stages[name] = prepare(project)
    return {"status": "PREPARED_REQUIRES_TECHNICAL_CHECK",
            "candidate": stages["exploratory"]["candidate"],
            "stages": stages, "scientific_approval": False,
            "next_action": "Run one-ligand technical check, review its outcome, then approve."}
