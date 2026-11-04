"""Human-vetted five-step SIDock workflow.

The guided workflow is independent of the legacy fully automated StageRunner.
Every approval is bound to input/artifact SHA256 hashes. Changing an upstream
artifact requires a fresh approval; it never silently deletes old docking runs.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

from .intake import (
    WATER, COMMON_IONS, COMMON_BUFFER,
    eligible_reference_groups, hetero_groups, _modified_polymer_keys,
)
from .manual_docking import (
    root as manual_root, settings, save_settings, install_receptor,
    prepared_ligands, _run_identity, active_run, _read_json, _save_json,
)
from .model import Project
from .state import STAGES
from .gui_service import safe_project_name, workspace_path

STEPS = (2, 3, 4)


def _sha_file(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _digest(payload: dict) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


def _json(path: Path) -> dict:
    return _read_json(path)


def _source(project: Path, kind: str) -> Path | None:
    data = _json(project / "project.json")
    value = data.get("receptor_source" if kind == "protein" else "ligand_source")
    if value and Path(value).is_file():
        return Path(value)
    fallback = project / "input" / ("protein.pdb" if kind == "protein" else "compounds.sdf")
    return fallback if fallback.is_file() else None


def create_guided_project(workspace: Path, name: str) -> Path:
    project = workspace_path(workspace, safe_project_name(name))
    if project.exists():
        raise FileExistsError("A project with that name already exists.")
    project.mkdir(parents=True)
    try:
        Project(str(project)).save()
        (project / "input").mkdir(exist_ok=True)
        state = {"stages": {stage: {"status": "PENDING"} for stage in STAGES},
                 "compounds": {}, "updated": None}
        _save_json(project / "project_state.json", state)
        _save_json(workflow_path(project), {"schema": 1, "approvals": {}, "history": []})
    except Exception:
        shutil.rmtree(project)
        raise
    return project


def workflow_path(project: Path) -> Path:
    return manual_root(project) / "workflow.json"


def load_workflow(project: Path) -> dict:
    return {"schema": 1, "approvals": {}, "history": [], **_json(workflow_path(project))}


def _update_workflow(project: Path, value: dict) -> dict:
    _save_json(workflow_path(project), value)
    return value


def _event(project: Path, action: str, detail: dict | None = None):
    workflow = load_workflow(project)
    workflow["history"].append({
        "at": datetime.now(timezone.utc).isoformat(),
        "action": action, "detail": detail or {},
    })
    _update_workflow(project, workflow)


def save_input(project: Path, source: Path, kind: str) -> Path:
    project = project.resolve()
    if kind not in ("protein", "compounds"):
        raise ValueError("Unknown input category")
    source = Path(source)
    suffixes = (".pdb", ".ent") if kind == "protein" else (".sdf",)
    if source.suffix.lower() not in suffixes or not source.is_file() or source.stat().st_size == 0:
        raise ValueError(f"Invalid {kind} upload")
    existing = _source(project, "protein" if kind == "protein" else "compounds")
    if existing is not None:
        raise RuntimeError(f"Existing {kind} input cannot be overwritten. Create a new project for a different input.")
    destination = project / "input" / ("protein.pdb" if kind == "protein" else "compounds.sdf")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.resolve() != destination.resolve():
        shutil.copy2(source, destination)
    data_path = project / "project.json"
    data = _json(data_path)
    data["receptor_source" if kind == "protein" else "ligand_source"] = str(destination.resolve())
    _save_json(data_path, data)
    _event(project, f"UPLOAD_{kind.upper()}", {"sha256": _sha_file(destination)})
    return destination


def protein_profile(project: Path) -> dict:
    pdb = _source(project, "protein")
    if pdb is None:
        return {"available": False, "groups": [], "ligands": []}
    modified = _modified_polymer_keys(pdb)
    heteros = hetero_groups(pdb)
    ligands = {g.key for g in eligible_reference_groups(pdb)}
    # Any nonpolymer component with coordinates can be manually proposed as a
    # grid marker, including a likely protein cap. Its structural role must be
    # recorded and disclosed; proposal is not native-ligand validation.
    for group in heteros:
        atom_names = {atom[12:16].strip().upper() for atom in group.atoms}
        is_modified = ((group.resname, group.chain, group.resseq) in modified
                       or {'N', 'CA', 'C', 'O'}.issubset(atom_names))
        is_cap = group.resname == "ACE" and group.resseq in {"0", "1"}
        heavy = sum((atom[76:78].strip().upper() or atom[12:16].strip()[0:1].upper())
                    not in {"H", "D"} for atom in group.atoms)
        if (heavy >= 2 and group.resname not in WATER | COMMON_IONS
                and not is_modified):
            ligands.add(group.key)
    groups = []
    for group in heteros:
        key = group.key
        atom_names = {atom[12:16].strip().upper() for atom in group.atoms}
        is_modified = ((group.resname, group.chain, group.resseq) in modified
                       or {'N', 'CA', 'C', 'O'}.issubset(atom_names))
        is_cap = group.resname == "ACE" and group.resseq in {"0", "1"}
        if group.resname in WATER:
            role = "WATER"
        elif is_modified:
            role = "MODIFIED_PROTEIN"
        elif is_cap:
            role = "POTENTIAL_PROTEIN_CAP"
        elif group.resname in COMMON_IONS:
            role = "ION"
        elif group.resname in COMMON_BUFFER:
            role = "BUFFER_OR_SOLVENT"
        elif key in ligands:
            role = "LIGAND_CANDIDATE"
        else:
            role = "NONSTANDARD_OR_OTHER"
        groups.append({
            "key": key, "code": group.resname, "chain": group.chain,
            "number": group.resseq, "atoms": group.atom_count, "role": role,
            "default_remove": role in {"WATER", "BUFFER_OR_SOLVENT", "LIGAND_CANDIDATE"},
            "protected": is_modified,
        })
    headers = []
    for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("HEADER", "TITLE ", "COMPND")):
            headers.append(line.strip())
        if len(headers) == 10: break
    return {"available": True, "protein": str(pdb), "protein_sha256": _sha_file(pdb),
            "headers": headers, "groups": groups,
            "ligands": [x for x in groups if x["key"] in ligands]}


def _coords(lines: list[str]) -> list[tuple[float, float, float]]:
    out = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            try:
                out.append(tuple(float(line[idx:idx+8]) for idx in (30, 38, 46)))
            except ValueError:
                pass
    return out


def select_protein(
    project: Path, *,
    ligand_key: str | None,
    remove_groups: list[str],
    center: tuple[float, float, float] | None = None,
    remove_waters: bool = True,
    padding: float = 8.0,
    acknowledge_non_ligand_component: bool = False,
) -> dict:
    project = project.resolve()
    profile = protein_profile(project)
    if not profile["available"]:
        raise ValueError("Upload a protein PDB in Step 2.")
    valid_groups = {g["key"]: g for g in profile["groups"]}
    candidates = {g["key"] for g in profile["ligands"]}
    if ligand_key and ligand_key not in candidates:
        raise ValueError("Selected reference component must be available in the uploaded protein.")
    selected_role = valid_groups[ligand_key]["role"] if ligand_key else None
    unusual_site = selected_role in {"POTENTIAL_PROTEIN_CAP", "BUFFER_OR_SOLVENT", "NONSTANDARD_OR_OTHER"}
    if unusual_site and not acknowledge_non_ligand_component:
        raise ValueError(
            f"{ligand_key} is classified as {selected_role}, not a verified native inhibitor. "
            "Tick the scientific acknowledgment if you still want to use its coordinates."
        )
    removal = set(remove_groups)
    if removal - set(valid_groups):
        raise ValueError("Unknown residue specified for removal")
    if any(valid_groups[k]["protected"] for k in removal):
        raise ValueError("Modified protein residues cannot be removed by this workflow")
    if any(valid_groups[k]["role"] == "POTENTIAL_PROTEIN_CAP" for k in removal) and not acknowledge_non_ligand_component:
        raise ValueError("Removing a possible protein cap requires explicit scientific acknowledgment.")
    if ligand_key is not None and not unusual_site:
        removal.add(ligand_key)
    # A likely covalent cap should not be stripped automatically. The reviewer
    # may elect to remove it, but must acknowledge that this changes the protein.
    if remove_waters:
        removal.update(g["key"] for g in profile["groups"] if g["role"] == "WATER")
    else:
        removal.difference_update(g["key"] for g in profile["groups"] if g["role"] == "WATER")
    if not 1.0 <= padding <= 20.0:
        raise ValueError("Ligand box padding must be between 1 and 20 Å")

    pdb = _source(project, "protein")
    assert pdb is not None
    selected_lines = []
    kept = []
    removed_atoms = {}
    for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("CONECT"):
            continue  # broken inter-residue CONECT after deleting a ligand
        if line.startswith(("ATOM", "HETATM")):
            key = f"{line[17:20].strip().upper()}:{line[21:22].strip() or '-'}:{line[22:27].strip()}"
            if key == ligand_key:
                selected_lines.append(line)
            if line.startswith("HETATM") and key in removal:
                removed_atoms[key] = removed_atoms.get(key, 0) + 1
                continue
        kept.append(line)
    if ligand_key and not selected_lines:
        raise RuntimeError("Selected ligand has no coordinates in protein PDB")

    t = settings(project)
    if ligand_key:
        coords = _coords(selected_lines)
        if not coords:
            raise RuntimeError("Selected crystal ligand contains no numeric coordinates")
        mins = [min(v[idx] for v in coords) for idx in range(3)]
        maxs = [max(v[idx] for v in coords) for idx in range(3)]
        center = tuple(round((mins[idx]+maxs[idx])/2, 3) for idx in range(3))
        spans = [(maxs[idx] - mins[idx]) + 2*padding for idx in range(3)]
        npts = [int(math.ceil(span/t["spacing"]/2)*2) for span in spans]
        if any(v > 126 for v in npts):
            raise ValueError("Selected ligand plus padding exceeds AutoDock4 126-point grid limit; reduce padding or increase spacing.")
        npts = [max(20, n) for n in npts]
        mode = "USER_SELECTED_COMPONENT" if unusual_site else "CRYSTAL_LIGAND"
    else:
        if center is None or len(center) != 3 or not all(math.isfinite(float(x)) for x in center):
            raise ValueError("Without a crystal ligand, enter all three grid center coordinates manually.")
        center = tuple(float(x) for x in center)
        npts = [t["npts_x"], t["npts_y"], t["npts_z"]]
        mode = "MANUAL_CENTER"

    new_settings = {**t, "center_x": center[0], "center_y": center[1], "center_z": center[2],
                    "npts_x": npts[0], "npts_y": npts[1], "npts_z": npts[2]}
    new_settings = save_settings(project, new_settings)
    folder = manual_root(project) / "inputs"
    folder.mkdir(parents=True, exist_ok=True)
    clean = folder / "protein_clean.pdb"
    clean.write_text("\n".join(kept) + "\n", encoding="utf-8")
    crystal = folder / "crystal_ligand.pdb"
    if ligand_key:
        crystal.write_text("\n".join(selected_lines) + "\nEND\n", encoding="utf-8")
    elif crystal.exists():
        crystal.unlink()
    choice = {
        "mode": mode, "ligand_key": ligand_key,
        "component_role": selected_role,
        "user_selected_without_native_validation": bool(unusual_site),
        "scientific_acknowledgment": bool(acknowledge_non_ligand_component),
        "site_warning": (
            f"{ligand_key} is classified as {selected_role}; using its coordinates "
            "to place a grid does not establish a native inhibitor binding site. "
            "A possible N-terminal acetyl cap may be covalent and is retained "
            "unless specifically removed."
        ) if unusual_site else None,
        "protein_sha256": _sha_file(pdb),
        "clean_sha256": _sha_file(clean),
        "removed_residues": sorted(removal), "removed_atom_counts": removed_atoms,
        "retained_residues": [g["key"] for g in profile["groups"] if g["key"] not in removal],
        "remove_waters": remove_waters, "padding": padding,
        "center": list(center), "npts": npts, "spacing": new_settings["spacing"],
        "crystal_ligand_sha256": _sha_file(crystal),
    }
    _save_json(folder / "protein_choice.json", choice)
    _event(project, "PROTEIN_REVIEW_SAVED", choice)
    return choice


def prepare_protein(project: Path, *, exploratory_metal_free: bool = False) -> dict:
    """Prepare exactly the reviewed receptor; optional explicit metal-free branch.

    The scientific clean PDB and original input are NEVER modified. The
    exploratory branch removes nickel only from a distinct Meeko input file,
    keeps a separate log, and records which ions were excluded.
    """
    from .core import prepare_receptor
    from .manual_docking import settings as manual_settings
    project = project.resolve()
    choice = _json(manual_root(project) / "inputs/protein_choice.json")
    if not choice:
        raise RuntimeError("Review and save protein components before preparing receptor")
    src = manual_root(project) / "inputs/protein_clean.pdb"
    if _sha_file(src) != choice["clean_sha256"]:
        raise RuntimeError("The cleaned PDB changed after review; review the protein again")
    folder = manual_root(project) / "inputs"
    nickel_residues = sorted({
        f"{line[17:20].strip().upper()}:{line[21:22].strip() or '-'}:{line[22:27].strip()}"
        for line in src.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith("HETATM") and line[17:20].strip().upper() == "NI"
    })
    if nickel_residues and not exploratory_metal_free:
        raise RuntimeError(
            "Standard Meeko receptor preparation blocked before execution: "
            "the reviewed protein still contains unsupported nickel residues "
            + ", ".join(nickel_residues) + ". "
            "Choose a scientifically parameterized metal-aware receptor PDBQT "
            "for urease docking, or explicitly approve the exploratory metal-free "
            "preparation button in Step 2. No protein or docking files were changed."
        )
    removal_audit = None
    if exploratory_metal_free:
        if not manual_settings(project).get("allow_missing_metal"):
            raise RuntimeError(
                "Explicit exploratory metal-free acknowledgment required in Step 2"
            )
        lines = src.read_text(encoding="utf-8", errors="replace").splitlines()
        metals = [
            line for line in lines
            if line.startswith("HETATM") and line[17:20].strip().upper() == "NI"
        ]
        if not metals:
            raise RuntimeError(
                "Reviewed receptor has no crystallographic NI atoms. "
                "Use standard receptor preparation instead."
            )
        temp = folder / "protein_meeko_metal_free_input.pdb"
        temp.write_text(
            "\n".join(line for line in lines if line not in metals) + "\n",
            encoding="utf-8",
        )
        removal_audit = {
            "reason": "EXPLORATORY_METAL_FREE_MEEKO_PREPARATION",
            "excluded_metal_atoms": len(metals),
            "excluded_residues": sorted({
                f"{line[17:20].strip()}:{line[21:22].strip() or '-'}:{line[22:27].strip()}"
                for line in metals
            }),
            "reviewed_clean_pdb_sha256": _sha_file(src),
            "temporary_meeko_pdb_sha256": _sha_file(temp),
            "warning": (
                "Exploratory metal-free receptor: nickel coordination is absent; "
                "not a validated urease metal-site docking model."
            ),
        }
        src_for_meeko = temp
        target = folder / "protein_meeko_metal_free.pdbqt"
        log = manual_root(project) / "logs/receptor_meeko_metal_free.log"
        _save_json(folder / "metal_free_preparation_audit.json", removal_audit)
    else:
        src_for_meeko = src
        target = folder / "protein_meeko.pdbqt"
        log = manual_root(project) / "logs/receptor_meeko.log"
    prepared = prepare_receptor(src_for_meeko, target, log, allow_bad_res=False)
    if exploratory_metal_free:
        from .manual_docking import pdbqt_atom_types
        if any(t.upper() == "NI" for t in pdbqt_atom_types(prepared)):
            raise RuntimeError(
                "Metal-free Meeko output unexpectedly contains nickel; refusing to install"
            )
    install_reviewed_receptor(project, prepared)
    result = {
        "prepared_receptor": str(prepared),
        "installed_receptor": str(manual_root(project) / "receptor.pdbqt"),
        "log": str(log),
        "exploratory_metal_free": bool(exploratory_metal_free),
        "metal_removal_audit": removal_audit,
        "note": "Metal-free preparation cannot validate native Ni coordination."
                if exploratory_metal_free else "No metal treatment or force-field validation performed.",
    }
    _event(project, "RECEPTOR_PREPARED", {
        "exploratory_metal_free": bool(exploratory_metal_free),
        "receptor_sha256": _sha_file(manual_root(project) / "receptor.pdbqt"),
        "audit": removal_audit,
    })
    return result


def install_reviewed_receptor(project: Path, source: Path) -> Path:
    """Bind an imported/prepared PDBQT to the explicitly reviewed cleaned protein."""
    project = Path(project).resolve()
    choice = _json(manual_root(project) / "inputs/protein_choice.json")
    clean = manual_root(project) / "inputs/protein_clean.pdb"
    if not choice or choice.get("clean_sha256") != _sha_file(clean):
        raise RuntimeError("Save the current protein component selection before importing a receptor")
    dest = install_receptor(project, source)
    _save_json(manual_root(project) / "inputs/receptor_review.json", {
        "clean_sha256": _sha_file(clean),
        "receptor_sha256": _sha_file(dest),
        "source_path": str(Path(source).resolve()),
        "at": datetime.now(timezone.utc).isoformat(),
    })
    _event(project, "RECEPTOR_IMPORTED", {"sha256": _sha_file(dest)})
    return dest


def compound_profile(project: Path) -> dict:
    sdf = _source(project, "compounds")
    if sdf is None:
        return {"available": False, "total": 0, "errors": []}
    from rdkit import Chem
    supplier = Chem.SDMolSupplier(str(sdf), sanitize=True, removeHs=False)
    labels = []
    errors = []
    for i, mol in enumerate(supplier, 1):
        if mol is None:
            errors.append(f"Compound {i}: failed RDKit parsing")
            continue
        name = (mol.GetProp("Compound_ID") if mol.HasProp("Compound_ID") else
                mol.GetProp("_Name") if mol.HasProp("_Name") else "").strip()
        if not name:
            errors.append(f"Compound {i}: missing name/Compound_ID")
        labels.append(name)
    if len(set(labels)) != len(labels):
        errors.append("Duplicate compound identifiers")
    normalized = [re.sub(r"[^A-Za-z0-9_.-]+", "_", name) for name in labels]
    if len(set(normalized)) != len(normalized):
        errors.append("Compound identifiers collide after filename normalization")
    if not labels:
        errors.append("No valid SDF compounds")
    manifest = _json(manual_root(project) / "ligands/manifest.json")
    ligands = prepared_ligands(project)
    fingerprint = _sha_file(sdf)
    prepared = manifest.get("source_sha256") == fingerprint and len(ligands) == len(labels) and (
        set(manifest.get("compounds", [])) == {x.name for x in ligands}
    )
    return {"available": True, "source": str(sdf), "source_sha256": fingerprint,
            "total": len(labels), "names": labels, "errors": errors,
            "prepared": bool(prepared), "prepared_count": len(ligands)}


def import_existing_compounds(project: Path) -> dict:
    """Copy already prepared same-project PDBQTs without running minimization.

    Reuse requires an exact match of the expected sanitized SDF identifiers.
    Original PDBQT files and any docking outputs are never altered.
    """
    project = Path(project).resolve()
    profile = compound_profile(project)
    if not profile.get("available") or profile.get("errors"):
        raise RuntimeError("Existing compound PDBQTs can only be reused for a valid, uniquely named SDF")
    names = [re.sub(r"[^A-Za-z0-9_.-]+", "_", name) for name in profile["names"]]
    originals = project / "ligands" / "prepared"
    missing = [name for name in names
               if not (originals / (name + ".pdbqt")).is_file()
               or (originals / (name + ".pdbqt")).stat().st_size == 0]
    if missing:
        raise RuntimeError("Existing prepared PDBQT files missing for: " + ", ".join(missing))
    target = manual_root(project) / "ligands"
    manifest = target / "manifest.json"
    previous = _json(manifest)
    prepared = prepared_ligands(project)
    if prepared and (previous.get("source_sha256") != profile["source_sha256"]
                     or {x.name for x in prepared} != {name+".pdbqt" for name in names}):
        raise RuntimeError("Manual ligand preparation already exists and differs; refused to overwrite")
    (target / "prepared").mkdir(parents=True, exist_ok=True)
    (target / "source").mkdir(parents=True, exist_ok=True)
    for name in names:
        destination = target / "prepared" / (name + ".pdbqt")
        if not destination.is_file():
            shutil.copy2(originals / (name+".pdbqt"), destination)
        minimized = project / "ligands" / "source" / (name + ".sdf")
        if minimized.is_file() and not (target/"source"/(name+".sdf")).exists():
            shutil.copy2(minimized, target/"source"/(name+".sdf"))
    _save_json(manifest, {
        "source_sha256": profile["source_sha256"],
        "compounds": [name+".pdbqt" for name in names],
        "preparation_method": "IMPORTED_EXISTING_PROJECT_PDBQT_NO_REMINIMIZATION",
    })
    _event(project, "REUSE_COMPOUND_PDBQT", {"count": len(names)})
    return {"count": len(names), "method": "imported_existing_PDBQT", "source_minimization_rerun": False}


def _fingerprint(project: Path, step: int) -> str | None:
    if step == 2:
        src = _source(project, "protein")
        choice = manual_root(project) / "inputs/protein_choice.json"
        clean = manual_root(project) / "inputs/protein_clean.pdb"
        receptor = manual_root(project) / "receptor.pdbqt"
        if not all(x and x.is_file() for x in (src, choice, clean, receptor)):
            return None
        receptor_review = _json(manual_root(project) / "inputs/receptor_review.json")
        if (receptor_review.get("clean_sha256") != _sha_file(clean) or
                receptor_review.get("receptor_sha256") != _sha_file(receptor)):
            return None
        metal_audit = None
        metal_parameter_hash = None
        if receptor_review.get("metal_model"):
            audit = manual_root(project) / "inputs/Ni_reconstruction_audit.json"
            param_file = manual_root(project) / "custom_parameters.dat"
            if (not audit.is_file() or not param_file.is_file()
                    or receptor_review.get("metal_audit_sha256") != _sha_file(audit)):
                return None
            metal_data = _json(audit)
            if (metal_data.get("reconstructed_sha256") != _sha_file(receptor)
                    or metal_data.get("reviewed_clean_sha256") != _sha_file(clean)
                    or metal_data.get("AD4_parameters_sha256") != _sha_file(param_file)):
                return None
            metal_audit = _sha_file(audit)
            metal_parameter_hash = _sha_file(param_file)
        return _digest({
            "source": _sha_file(src), "choice": _sha_file(choice),
            "clean": _sha_file(clean), "receptor": _sha_file(receptor),
            "metal_audit": metal_audit, "metal_parameters": metal_parameter_hash,
            "exploratory_metal_free": settings(project)["allow_missing_metal"],
        })
    if step == 3:
        source = _source(project, "compounds")
        manifest = manual_root(project) / "ligands/manifest.json"
        if source is None or not manifest.is_file():
            return None
        files = prepared_ligands(project)
        if not files:
            return None
        return _digest({
            "source": _sha_file(source), "manifest": _sha_file(manifest),
            "prepared": [(p.name, _sha_file(p)) for p in files],
        })
    if step == 4:
        run = active_run(project)
        if run is None:
            return None
        provenance = run / "grid/production/grid_provenance.json"
        if not provenance.is_file():
            return None
        settings_path = manual_root(project) / "settings.json"
        return _digest({
            "settings": _sha_file(settings_path), "grid_provenance": _sha_file(provenance),
            "receptor": _fingerprint(project, 2), "compounds": _fingerprint(project, 3),
        })
    raise ValueError("Unexpected step")


def approved(project: Path, step: int) -> bool:
    entry = load_workflow(project).get("approvals", {}).get(str(step), {})
    current = _fingerprint(project, step)
    if not current or entry.get("digest") != current:
        return False
    if step >= 3 and not approved(project, 2):
        return False
    if step == 4 and not approved(project, 3):
        return False
    return True


def approve(project: Path, step: int) -> dict:
    if step not in STEPS:
        raise ValueError("Step must be 2, 3 or 4")
    if step >= 3 and not approved(project, step - 1):
        raise RuntimeError(f"Step {step-1} must be approved first.")
    if step == 2:
        choice = _json(manual_root(project) / "inputs/protein_choice.json")
        if not choice:
            raise RuntimeError("Select a crystal ligand or manual coordinates and save protein review")
        if not (manual_root(project) / "receptor.pdbqt").is_file():
            raise RuntimeError("Prepare or import a receptor PDBQT before approving protein")
        receptor_review = _json(manual_root(project) / "inputs/receptor_review.json")
        if (receptor_review.get("clean_sha256") != _sha_file(manual_root(project) / "inputs/protein_clean.pdb")
                or receptor_review.get("receptor_sha256") != _sha_file(manual_root(project) / "receptor.pdbqt")):
            raise RuntimeError("Receptor PDBQT is not associated with the current reviewed protein. Reimport or reprepare it.")
        receptor = manual_root(project) / "receptor.pdbqt"
        if _sha_file(receptor) is None:
            raise RuntimeError("Invalid receptor")
        # Restored nickel remains an unvalidated model, but must be backed by
        # the same audited custom force-field library that will reach AutoGrid.
        if receptor_review.get("metal_model") and _fingerprint(project, 2) is None:
            raise RuntimeError(
                "Nickel reconstruction audit, receptor or custom parameters changed; "
                "reconstruct or reimport and review before approving Step 2."
            )
        # Do not disguise a metal-stripped receptor as a validated metalloprotein.
        from .manual_docking import _input_metals, _pdbqt_metals
        if _input_metals(project) and _pdbqt_metals(project) < _input_metals(project) and not settings(project)["allow_missing_metal"]:
            raise RuntimeError("Receptor is missing crystallographic nickel. Configure explicitly as exploratory metal-free, or import a chemically parameterized metal receptor.")
    if step == 3:
        profile = compound_profile(project)
        if profile.get("errors"):
            raise RuntimeError("Compound issues: " + "; ".join(profile["errors"]))
        if not profile.get("prepared"):
            raise RuntimeError("Ligands must be prepared and verified before approval")
    if step == 4:
        if not active_run(project) or not (active_run(project) / "grid/production/grid_provenance.json").is_file():
            raise RuntimeError("Generate and inspect AutoGrid maps before approving docking setup")
    digest = _fingerprint(project, step)
    if digest is None:
        raise RuntimeError(f"Step {step} inputs or outputs not ready")
    current = load_workflow(project)
    current["approvals"][str(step)] = {
        "digest": digest, "at": datetime.now(timezone.utc).isoformat(),
    }
    for later in STEPS:
        if later > step:
            current["approvals"].pop(str(later), None)
    _update_workflow(project, current)
    _event(project, "APPROVE_STEP", {"step": step, "digest": digest})
    return current["approvals"][str(step)]


def ready_to_dock(project: Path) -> None:
    for step in STEPS:
        if not approved(project, step):
            raise RuntimeError(f"Step {step} must be approved with current inputs before docking.")


def launch_docking(project: Path) -> dict:
    ready_to_dock(project)
    from .manual_docking import launch
    return launch(project, "dock")
