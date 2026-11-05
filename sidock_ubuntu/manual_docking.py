"""Explicit AutoDock4 workflow separate from SIDock's automated scientific stages.

Manual protocol: import a PDBQT receptor, prepare SDF ligands, specify grid
center and grid points, run AutoGrid4, run AutoDock4. The automated project
state and all prior DLGs remain unchanged. This is not a native validation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import fcntl
from contextlib import contextmanager
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

from .core import (
    pdbqt_atom_types, prepare_ligand, resolve_executable,
    split_and_minimize_sdf,
)
from .model import DockingParameters, GridBox
from .stages import run_compound_docking, run_grid_stage


def root(project: Path) -> Path:
    return Path(project).resolve() / "manual_docking"


def _read_json(path: Path, default=None):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return default if default is not None else {}


def _save_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")
    temp.replace(path)


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _crystal_center(project: Path) -> tuple[float, float, float] | None:
    # This is only a suggested grid center, not automatic ligand selection
    # or an assertion that a modified protein residue is a native ligand.
    pdb = Path(project) / "input" / "protein.pdb"
    reference = Path(project) / "receptor" / "reference_ligand.pdb"
    coords = []
    if pdb.exists():
        for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("HETATM") and line[17:20].strip() == "HAE":
                try:
                    coords.append(tuple(float(line[i:i+8]) for i in (30, 38, 46)))
                except ValueError:
                    pass
    if not coords and reference.exists():
        for line in reference.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith(("ATOM", "HETATM")):
                try:
                    coords.append(tuple(float(line[i:i+8]) for i in (30, 38, 46)))
                except ValueError:
                    pass
    if not coords:
        return None
    return tuple(round((min(v[i] for v in coords) + max(v[i] for v in coords)) / 2, 3) for i in range(3))


def settings(project: Path) -> dict:
    existing = _read_json(root(project) / "settings.json")
    suggested = _crystal_center(project) or (0., 0., 0.)
    defaults = {
        "center_x": suggested[0], "center_y": suggested[1], "center_z": suggested[2],
        "npts_x": 80, "npts_y": 80, "npts_z": 80,
        "spacing": 0.375, "ga_runs": 150, "parallel_jobs": 10,
        "ga_num_evals": 2500000, "ga_pop_size": 150,
        "ga_num_generations": 27000, "allow_missing_metal": False,
    }
    return {**defaults, **existing}


def validate_settings(value: dict) -> dict:
    val = {}
    for key in ("center_x", "center_y", "center_z", "spacing"):
        number = float(value[key])
        if not math.isfinite(number):
            raise ValueError(f"{key} must be finite")
        val[key] = number
    if not 0.1 <= val["spacing"] <= 1.0:
        raise ValueError("spacing must be between 0.1 and 1.0 Å")
    for key in ("npts_x", "npts_y", "npts_z"):
        val[key] = int(value[key])
        if not 20 <= val[key] <= 126 or val[key] % 2:
            raise ValueError(f"{key} must be even, 20 to 126")
    for key, limit in (("ga_runs", 10000), ("parallel_jobs", 32),
                       ("ga_num_evals", 50000000), ("ga_pop_size", 1000),
                       ("ga_num_generations", 1000000)):
        val[key] = int(value[key])
        if not 1 <= val[key] <= limit:
            raise ValueError(f"{key} must be between 1 and {limit}")
    val["allow_missing_metal"] = value.get("allow_missing_metal") in (True, "true", "on", "yes", "1")
    return val


def save_settings(project: Path, new: dict) -> dict:
    value = validate_settings(new)
    _save_json(root(project) / "settings.json", value)
    return value


def _input_metals(project: Path) -> int:
    pdb = Path(project) / "input" / "protein.pdb"
    if not pdb.exists():
        return 0
    return sum(1 for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines()
               if line.startswith("HETATM") and line[17:20].strip() in {"NI"})


def _pdbqt_metals(project: Path, receptor: Path | None = None) -> int:
    path = receptor if receptor is not None else root(project) / "receptor.pdbqt"
    if not path.exists():
        return 0
    return sum(1 for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
               if line.startswith(("ATOM", "HETATM"))
               and (line[77:79].strip().upper() == "NI"
                    or line[17:20].strip() == "NI"))


def parameter_file_path(project: Path) -> Path | None:
    file = root(project) / "custom_parameters.dat"
    return file if file.is_file() else None


def install_parameter_file(project: Path, source: Path) -> Path:
    """Accept only a complete user-supplied AD4 force field library.

    Nickel parameters are never guessed. A partial patch containing only Ni
    is insufficient because AD4's parameter_file replaces the default table.
    """
    source = Path(source)
    if not source.is_file() or source.stat().st_size == 0:
        raise ValueError("Select an existing AutoDock4 parameter library file")
    if source.stat().st_size > 2_000_000:
        raise ValueError("Custom parameter file is unreasonably large")
    lines = source.read_text(encoding="utf-8", errors="strict").splitlines()
    entries = set()
    for line in lines:
        parts = line.split()
        if parts and parts[0] == "atom_par" and len(parts) >= 12:
            entries.add(parts[1])
    if not {"C", "N", "OA", "HD"}.issubset(entries):
        raise ValueError("Supply a full AD4 parameter file with standard C, N, OA and HD atom_par records, not a Ni-only patch")
    receptor = root(project) / "receptor.pdbqt"
    present_types = set(pdbqt_atom_types(receptor)) if receptor.is_file() else set()
    missing = {t for t in present_types if t not in entries}
    if missing:
        raise ValueError("Parameter library lacks receptor atom types: " + ", ".join(sorted(missing)))
    target = root(project) / "custom_parameters.dat"
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.resolve() != target.resolve():
        shutil.copy2(source, target)
    return target


def install_receptor(project: Path, source: Path) -> Path:
    project = Path(project).resolve()
    source = Path(source).resolve()
    if source.suffix.lower() != ".pdbqt" or not source.is_file():
        raise ValueError("A prepared receptor .pdbqt file is required")
    if not any(line.startswith(("ATOM", "HETATM")) for line in source.read_text(encoding="utf-8", errors="replace").splitlines()):
        raise ValueError("The uploaded PDBQT has no receptor atoms")
    target = root(project) / "receptor.pdbqt"
    target.parent.mkdir(parents=True, exist_ok=True)
    if source != target:
        shutil.copy2(source, target)
    return target


def _ligand_sdf(project: Path) -> Path:
    data = _read_json(Path(project) / "project.json")
    path = Path(data.get("ligand_source", ""))
    if not path.is_file():
        path = Path(project) / "input" / "compounds.sdf"
    if not path.is_file():
        raise FileNotFoundError("Project compound SDF not found")
    return path


def prepared_ligands(project: Path) -> list[Path]:
    return sorted((root(project) / "ligands" / "prepared").glob("*.pdbqt"))


def prepare(project: Path) -> dict:
    """Independently prepare the complete SDF library; never run docking.

    An interrupted run may leave successful individual ligands. Keep those
    outputs and resume missing compounds without deleting or overwriting any.
    Only write a final manifest when every expected ligand is verified.
    """
    from .guided_workflow import compound_profile
    project = Path(project).resolve()
    out = root(project)
    profile = compound_profile(project)
    if not profile.get("available") or profile.get("errors"):
        raise RuntimeError(
            "SDF must contain valid, uniquely identified compounds: "
            + "; ".join(profile.get("errors", []))
        )
    sdf = _ligand_sdf(project)
    fingerprint = _sha(sdf)
    expected = [re.sub(r"[^A-Za-z0-9_.-]+", "_", n) + ".pdbqt"
                for n in profile["names"]]
    expected_names = set(expected)
    manifest_path = out / "ligands" / "manifest.json"
    manifest = _read_json(manifest_path)
    identity_path = out / "ligands" / "preparation_identity.json"
    prior_identity = _read_json(identity_path)
    old = prepared_ligands(project)
    present = {p.name for p in old}
    minimized_dir = out / "ligands" / "source"
    old_minimized = list(minimized_dir.glob("*.sdf"))
    saved_hash = prior_identity.get("source_sha256")
    if saved_hash is not None and saved_hash != fingerprint:
        raise RuntimeError(
            "The source SDF changed since the previous preparation attempt. "
            "Existing minimized and PDBQT files were preserved."
        )
    if (saved_hash is None and (old or old_minimized)
            and manifest.get("source_sha256") != fingerprint):
        raise RuntimeError(
            "Existing partial ligand outputs have no matching source-SDF provenance. "
            "Nothing was overwritten; review the files before resuming."
        )
    if old and not present.issubset(expected_names):
        raise RuntimeError(
            "Unexpected existing prepared ligands were found. Preserved all files; "
            "review the conflicting source SDF and manual_docking/ligands/prepared."
        )
    if old and manifest.get("source_sha256") not in (None, fingerprint):
        raise RuntimeError(
            "Ligand source changed; existing manual prepared ligands preserved. "
            "Review the source and outputs before retrying."
        )
    def valid_pdbqt(path: Path) -> bool:
        return (path.is_file() and path.stat().st_size > 0
                and any(line.startswith(("ATOM", "HETATM")) for line in
                        path.read_text(encoding="utf-8", errors="replace").splitlines()))
    if (expected_names == present
            and all(valid_pdbqt(p) for p in old)
            and manifest.get("source_sha256") == fingerprint
            and set(manifest.get("compounds", [])) == expected_names):
        return {"status": "REUSED", "prepared": len(expected),
                "ligands": sorted(expected_names)}
    if saved_hash is None:
        _save_json(identity_path, {
            "source_sha256": fingerprint,
            "expected_compounds": sorted(expected_names),
        })
    sources = split_and_minimize_sdf(
        sdf, out / "ligands" / "source", reuse_existing=True
    )
    actual_names = {p.stem + ".pdbqt" for p in sources}
    if actual_names != expected_names or len(sources) != len(expected):
        missing = expected_names - actual_names
        raise RuntimeError(
            f"Minimization incomplete: generated {len(sources)} of {len(expected)} SDFs. "
            "Missing: " + ", ".join(sorted(missing))
            + ". Existing files were preserved; review SDF/embedding failures."
        )
    prepared = []
    for source in sources:
        dest = out / "ligands" / "prepared" / (source.stem + ".pdbqt")
        if dest.is_file() and not valid_pdbqt(dest):
            raise RuntimeError(
                f"Existing ligand PDBQT is invalid: {dest}; preserved for inspection. "
                "Fix or archive it explicitly before retrying."
            )
        if not dest.is_file():
            prepare_ligand(source, dest, out / "logs" / ("prepare_" + source.stem + ".log"))
        if not valid_pdbqt(dest):
            raise RuntimeError(f"Ligand output is empty or invalid: {dest}")
        prepared.append(dest.name)
    if set(prepared) != expected_names:
        raise RuntimeError("Preparation did not produce all expected PDBQTs")
    _save_json(manifest_path, {
        "source_sha256": fingerprint,
        "compounds": sorted(prepared),
        "preparation_method": "OPENBABEL_LIGAND_SDF_PREPARATION",
    })
    return {"status": "PREPARED", "prepared": len(prepared), "ligands": sorted(prepared)}


def selected_receptor_for_run(project: Path) -> tuple[Path, str]:
    """Use hash-approved exploratory candidate, or leave standard receptor intact."""
    from .guided_workflow import load_workflow, approved
    from .exploratory_approval import MODE
    entry = load_workflow(project).get("approvals", {}).get("2", {})
    if entry.get("mode") == MODE:
        if not approved(project, 2):
            raise RuntimeError(
                "Exploratory receptor approval expired: source, model, parameters "
                "or smoke-test audit changed. Review and approve Step 2 again."
            )
        return root(project) / "inputs/receptor_openbabel_exploratory.pdbqt", MODE
    return root(project) / "receptor.pdbqt", "MANUAL_AUTODOCK4"


def _run_identity(project: Path) -> tuple[str, list[Path]]:
    project = Path(project).resolve()
    receptor, mode = selected_receptor_for_run(project)
    if not receptor.is_file():
        raise RuntimeError("Import or technically approve a prepared receptor PDBQT first")
    ligands = prepared_ligands(project)
    if not ligands:
        raise RuntimeError("Prepare the SDF ligands first")
    m = _read_json(root(project) / "ligands" / "manifest.json")
    if set(m.get("compounds", [])) != set(p.name for p in ligands):
        raise RuntimeError("Prepared-ligand manifest does not match files")
    if m.get("source_sha256") != _sha(_ligand_sdf(project)):
        raise RuntimeError("SDF has changed since ligand preparation")
    current = settings(project)
    if _input_metals(project) and not _pdbqt_metals(project, receptor) and not current["allow_missing_metal"]:
        raise RuntimeError(
            "The original structure contains catalytic nickel but this receptor PDBQT does not. "
            "Import a metal-parameterized receptor, or explicitly enable exploratory metal-free docking."
        )
    values = {k: current[k] for k in ("center_x","center_y","center_z",
                "npts_x","npts_y","npts_z","spacing", "ga_runs",
                "ga_num_evals", "ga_pop_size", "ga_num_generations")}
    custom = parameter_file_path(project)
    metal_types = {x for x in pdbqt_atom_types(receptor) if x.upper() == "NI"}
    if metal_types and custom is None:
        raise RuntimeError(
            "Nickel atom type found in the receptor but no custom AutoDock4 parameter library is installed. "
            "Upload a scientifically justified complete AD4 parameter file in Step 4."
        )
    if custom is not None:
        from .nickel_reconstruction import _atom_parameter_types
        declared = _atom_parameter_types(custom)
        required = set(pdbqt_atom_types(receptor))
        for ligand in ligands:
            required.update(pdbqt_atom_types(ligand))
        missing_types = required - declared
        if missing_types:
            raise RuntimeError(
                "Custom AutoDock4 parameter library missing receptor/ligand atom types: "
                + ", ".join(sorted(missing_types))
            )
    metal_review = _read_json(root(project) / "inputs/receptor_review.json")
    if mode != "MANUAL_AUTODOCK4":
        # Technical exploratory mode is independently audited by Step 2.
        # The Ni reconstruction audit belongs only to the scientific receptor.
        metal_review = {}
    if metal_review.get("metal_model"):
        audit_path = root(project) / "inputs/Ni_reconstruction_audit.json"
        audit = _read_json(audit_path)
        if (custom is None
                or audit.get("AD4_parameters_sha256") != _sha(custom)
                or audit.get("reconstructed_sha256") != _sha(receptor)
                or metal_review.get("metal_audit_sha256") != _sha(audit_path)):
            raise RuntimeError(
                "Nickel reconstruction or custom force-field SHA changed. "
                "Reconstruct and vet the receptor before generating a grid or docking."
            )
    values["custom_parameters_sha256"] = _sha(custom) if custom is not None else None
    values["receptor_sha256"] = _sha(receptor)
    values["receptor_mode"] = mode
    if mode != "MANUAL_AUTODOCK4":
        from .guided_workflow import load_workflow
        values["exploratory_step2_digest"] = load_workflow(project)["approvals"]["2"]["digest"]
    values["ligands"] = {x.name: _sha(x) for x in ligands}
    token = hashlib.sha256(json.dumps(values, sort_keys=True).encode()).hexdigest()[:16]
    return token, ligands


def active_run(project: Path) -> Path | None:
    try:
        token, _ = _run_identity(project)
    except (OSError, RuntimeError, ValueError):
        return None
    return root(project) / "runs" / token


def grid(project: Path) -> dict:
    project = Path(project).resolve()
    token, ligands = _run_identity(project)
    _, model_status = selected_receptor_for_run(project)
    if model_status != "MANUAL_AUTODOCK4":
        from .guided_workflow import approved
        if not approved(project, 2) or not approved(project, 3):
            raise RuntimeError("Current exploratory Step 2 and Step 3 approvals required.")
    run = root(project) / "runs" / token
    grid_path = run / "grid" / "production" / "grid_provenance.json"
    if grid_path.is_file():
        return {"run_id": token, "status": "REUSED", "grid": str(grid_path)}
    target = run / "receptor" / "receptor.pdbqt"
    target.parent.mkdir(parents=True, exist_ok=True)
    active_receptor, model_status = selected_receptor_for_run(project)
    shutil.copy2(active_receptor, target)
    t = settings(project)
    box = GridBox(t["center_x"], t["center_y"], t["center_z"],
                  t["npts_x"], t["npts_y"], t["npts_z"], t["spacing"])
    types = sorted({atom_type for ligand in ligands for atom_type in pdbqt_atom_types(ligand)})
    try:
        result = run_grid_stage(target, run / "grid" / "production", box,
                                types, resolve_executable("autogrid4", "AutoGrid4"),
                                parameter_file=parameter_file_path(project))
    except Exception as exc:
        logfile = run / "grid" / "production" / "autogrid.log"
        details = logfile.read_text(encoding="utf-8", errors="replace")[-1600:] if logfile.exists() else ""
        raise RuntimeError(f"AutoGrid failed: {exc}. Last AutoGrid output: {details}") from exc
    _save_json(run / "protocol.json", {
        "run_id": token, "settings": t,
        "receptor_sha256": _sha(active_receptor),
        "model_status": model_status,
        "scientific_approval": False,  # Manual/technical operation does not prove scientific validation
        "ligand_sha256": {x.name: _sha(x) for x in ligands},
        "metal_atoms_original": _input_metals(project),
        "metal_atoms_prepared": _pdbqt_metals(project, active_receptor),
        "nickel_reconstruction_audit_sha256": _sha(root(project) / "inputs/Ni_reconstruction_audit.json")
              if (root(project) / "inputs/Ni_reconstruction_audit.json").is_file()
              and _read_json(root(project) / "inputs/receptor_review.json").get("metal_model")
              else None,
        "metal_model_validation": "NOT_PERFORMED",
        "validation": "NOT_PERFORMED",
        "custom_parameter_sha256": _sha(parameter_file_path(project)) if parameter_file_path(project) else None,
        "method": "AUTODOCK4_TECHNICAL_EXPLORATORY" if model_status != "MANUAL_AUTODOCK4"
                  else "MANUAL_AUTODOCK4",
    })
    return {"run_id": token, "status": "CREATED", "grid": str(result["gpf"])}


def complete_dlg_energy(path: Path) -> float | None:
    """Do not reuse truncated logs containing an early energy but no completion marker."""
    from .core import parse_best_energy
    if not path.is_file():
        return None
    content = path.read_text(encoding="utf-8", errors="replace")
    if "Successful Completion" not in content:
        return None
    return parse_best_energy(path)


def dock(project: Path) -> dict:
    project = Path(project).resolve()
    token, ligands = _run_identity(project)
    _, model_status = selected_receptor_for_run(project)
    if model_status != "MANUAL_AUTODOCK4":
        from .guided_workflow import ready_to_dock
        ready_to_dock(project)  # Includes fresh hash-bound Step 4 grid approval.
    run = root(project) / "runs" / token
    prov = run / "grid" / "production" / "grid_provenance.json"
    if not prov.is_file():
        raise RuntimeError("Generate AutoGrid maps for the current receptor/grid settings first")
    if _read_json(prov).get("receptor_sha256") != _sha(run / "grid" / "production" / "receptor.pdbqt"):
        raise RuntimeError("Grid receptor hash mismatch; refusing docking")
    t = settings(project)
    params = DockingParameters(
        ga_runs=t["ga_runs"], ga_num_evals=t["ga_num_evals"],
        ga_pop_size=t["ga_pop_size"],
        ga_num_generations=t["ga_num_generations"],
        parallel_jobs=t["parallel_jobs"],
    )
    dpf_dir = run / "docking" / "dpf"
    dlg_dir = run / "docking" / "dlg"
    work_dir = run / "docking" / "work"
    pending = [p for p in ligands if complete_dlg_energy(dlg_dir / (p.stem + ".dlg")) is None]
    for ligand in pending:
        old_dlg=dlg_dir / (ligand.stem + ".dlg")
        if old_dlg.is_file() and old_dlg.stat().st_size:
            archive=run / "docking" / "incomplete_archives"
            archive.mkdir(parents=True, exist_ok=True)
            fingerprint=_sha(old_dlg)[:12]
            backup=archive/(ligand.stem+"."+fingerprint+".incomplete.dlg")
            if not backup.is_file():
                shutil.copy2(old_dlg, backup)
    errors = []
    def one(ligand):
        return run_compound_docking(ligand, "receptor", work_dir, dpf_dir,
                                     dlg_dir, pdbqt_atom_types(ligand),
                                     params, resolve_executable("autodock4", "AutoDock4"),
                                     parameter_file=parameter_file_path(project))
    with ThreadPoolExecutor(max_workers=min(len(pending) or 1, t["parallel_jobs"])) as pool:
        futures = {pool.submit(one, ligand): ligand for ligand in pending}
        for future in as_completed(futures):
            ligand = futures[future]
            try:
                future.result()
            except Exception as exc:
                errors.append({"ligand": ligand.stem, "error": str(exc)})
    import csv
    energies = {}
    for ligand in ligands:
        dlg_file = dlg_dir / (ligand.stem + ".dlg")
        energies[ligand.stem] = complete_dlg_energy(dlg_file)
    results_csv = run / "results.csv"
    with results_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(("compound_id", "lowest_docking_energy_kcal_mol", "status",
                         "model_status", "scientific_approval", "dlg_path"))
        for ligand in ligands:
            energy = energies[ligand.stem]
            dlg_file = dlg_dir / (ligand.stem + ".dlg")
            writer.writerow((ligand.stem, energy if energy is not None else "",
                             "COMPLETE" if energy is not None else "FAILED",
                             selected_receptor_for_run(project)[1], False,
                             str(dlg_file) if dlg_file.exists() else ""))
    _, model_status = selected_receptor_for_run(project)
    result = {"run_id": token, "model_status": model_status,
              "scientific_approval": False,  # Manual/technical operation does not prove scientific validation
              "total": len(ligands),
              "already_finished": len(ligands) - len(pending),
              "launched": len(pending),
              "completed": sum(energy is not None for energy in energies.values()),
              "results_csv": str(results_csv),
              "errors": errors}
    _save_json(run / "docking_result.json", result)
    if errors or result["completed"] != result["total"]:
        raise RuntimeError(
            f"AutoDock completed {result['completed']}/{result['total']} compounds; "
            f"{len(errors)} explicit job errors. See {run/'docking_result.json'}"
        )
    return result


def job_path(project: Path) -> Path:
    return root(project) / "job.json"


@contextmanager
def _job_file_lock(project: Path):
    """Serialize child status and parent launch writes across processes."""
    folder = root(project)
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / ".job.lock").open("a", encoding="utf-8") as fh:
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh.fileno(), fcntl.LOCK_UN)


def job_status(project: Path) -> dict:
    from .job_store import process_alive
    job = _read_json(job_path(project))
    if not job:
        return {}
    if job.get("status") == "RUNNING" and not process_alive(job.get("pid")):
        return {**job, "status": "EXITED", "error": "Manual worker exited without a recorded completion"}
    return job


def launch(project: Path, action: str) -> dict:
    if action not in {"prepare", "receptor", "receptor_metal_free", "receptor_ni_intermediate", "receptor_openbabel", "receptor_openbabel_pipeline", "receptor_openbabel_charged", "receptor_openbabel_hydrogen_audit", "receptor_openbabel_exploratory", "exploratory_smoke", "grid", "dock"}:
        raise ValueError("Unsupported manual action")
    with _job_file_lock(project):
        if job_status(project).get("status") == "RUNNING":
            raise RuntimeError("A manual SIDock action is already running")
        folder = root(project)
        logfile = folder / "manual_worker.log"
        with logfile.open("a", encoding="utf-8") as out:
            proc = subprocess.Popen(
                [sys.executable, "-m", "sidock_ubuntu.manual_docking",
                 str(Path(project).resolve()), action],
                cwd=str(Path(__file__).resolve().parent.parent),
                stdin=subprocess.DEVNULL, stdout=out, stderr=subprocess.STDOUT,
                start_new_session=True,
            )
        job = {"action": action, "status": "RUNNING", "pid": proc.pid,
               "started": datetime.now(timezone.utc).isoformat()}
        _save_json(job_path(project), job)
    return job


def execute(project: Path, action: str) -> dict:
    path = job_path(project)
    try:
        from .guided_workflow import prepare_protein
        from .nickel_reconstruction import prepare_intermediate
        from .openbabel_receptor import (generate_candidate, generate_protonated_charge_candidate,
                                         audit_and_correct_hydrogen_types, prepare_exploratory_receptor, prepare_openbabel_receptor)
        from .exploratory_docking import run_exploratory_smoke
        results = {"prepare": prepare, "receptor": prepare_protein,
                   "receptor_metal_free": lambda p: prepare_protein(p, exploratory_metal_free=True),
                   "receptor_ni_intermediate": prepare_intermediate,
                   "receptor_openbabel": generate_candidate,
                   "receptor_openbabel_pipeline": prepare_openbabel_receptor,
                   "receptor_openbabel_charged": generate_protonated_charge_candidate,
                   "receptor_openbabel_hydrogen_audit": audit_and_correct_hydrogen_types,
                   "receptor_openbabel_exploratory": prepare_exploratory_receptor,
                   "exploratory_smoke": run_exploratory_smoke,
                   "grid": grid, "dock": dock}[action](project)
        with _job_file_lock(project):
            job = _read_json(path)
            job.update({"status": "SUCCEEDED", "result": results,
                        "error": None,
                        "finished": datetime.now(timezone.utc).isoformat()})
            _save_json(path, job)
        return results
    except Exception as exc:
        with _job_file_lock(project):
            job = _read_json(path)
            job.update({"status": "FAILED", "error": str(exc),
                        "finished": datetime.now(timezone.utc).isoformat()})
            _save_json(path, job)
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description="Explicit manual AutoDock4 runner")
    parser.add_argument("project", type=Path)
    parser.add_argument("action", choices=("prepare", "receptor", "receptor_metal_free", "receptor_ni_intermediate", "receptor_openbabel", "receptor_openbabel_pipeline", "receptor_openbabel_charged", "receptor_openbabel_hydrogen_audit", "receptor_openbabel_exploratory", "exploratory_smoke", "grid", "dock", "status"))
    args = parser.parse_args()
    if args.action == "status":
        print(json.dumps({"settings": settings(args.project),
                          "job": job_status(args.project)}, indent=2))
    else:
        print(json.dumps(execute(args.project, args.action), indent=2))


if __name__ == "__main__":
    main()
