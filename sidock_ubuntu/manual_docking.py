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
import os
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


def _pdbqt_metals(project: Path) -> int:
    path = root(project) / "receptor.pdbqt"
    if not path.exists():
        return 0
    return sum(1 for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
               if line.startswith(("ATOM", "HETATM"))
               and (line[77:79].strip().upper() == "NI"
                    or line[17:20].strip() == "NI"))


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
    project = Path(project).resolve()
    out = root(project)
    sdf = _ligand_sdf(project)
    fingerprint = _sha(sdf)
    manifest_path = out / "ligands" / "manifest.json"
    manifest = _read_json(manifest_path)
    old = prepared_ligands(project)
    if old and manifest.get("source_sha256") != fingerprint:
        raise RuntimeError("Ligand source changed; existing manual prepared ligands preserved. Use a separate project instead of overwriting them.")
    sources = split_and_minimize_sdf(sdf, out / "ligands" / "source")
    if not sources:
        raise RuntimeError("No usable named compounds could be extracted from SDF")
    prepared = []
    for source in sources:
        dest = out / "ligands" / "prepared" / (source.stem + ".pdbqt")
        if not dest.is_file() or dest.stat().st_size == 0:
            prepare_ligand(source, dest, out / "logs" / ("prepare_" + source.stem + ".log"))
        prepared.append(dest.name)
    _save_json(manifest_path, {"source_sha256": fingerprint, "compounds": prepared})
    return {"prepared": len(prepared), "ligands": prepared}


def _run_identity(project: Path) -> tuple[str, list[Path]]:
    project = Path(project).resolve()
    receptor = root(project) / "receptor.pdbqt"
    if not receptor.is_file():
        raise RuntimeError("Import a prepared receptor PDBQT first")
    ligands = prepared_ligands(project)
    if not ligands:
        raise RuntimeError("Prepare the SDF ligands first")
    m = _read_json(root(project) / "ligands" / "manifest.json")
    if set(m.get("compounds", [])) != set(p.name for p in ligands):
        raise RuntimeError("Prepared-ligand manifest does not match files")
    if m.get("source_sha256") != _sha(_ligand_sdf(project)):
        raise RuntimeError("SDF has changed since ligand preparation")
    current = settings(project)
    if _input_metals(project) and not _pdbqt_metals(project) and not current["allow_missing_metal"]:
        raise RuntimeError(
            "The original structure contains catalytic nickel but this receptor PDBQT does not. "
            "Import a metal-parameterized receptor, or explicitly enable exploratory metal-free docking."
        )
    values = {k: current[k] for k in ("center_x","center_y","center_z",
                "npts_x","npts_y","npts_z","spacing")}
    values["receptor_sha256"] = _sha(receptor)
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
    run = root(project) / "runs" / token
    grid_path = run / "grid" / "production" / "grid_provenance.json"
    if grid_path.is_file():
        return {"run_id": token, "status": "REUSED", "grid": str(grid_path)}
    target = run / "receptor" / "receptor.pdbqt"
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(root(project) / "receptor.pdbqt", target)
    t = settings(project)
    box = GridBox(t["center_x"], t["center_y"], t["center_z"],
                  t["npts_x"], t["npts_y"], t["npts_z"], t["spacing"])
    types = sorted({atom_type for ligand in ligands for atom_type in pdbqt_atom_types(ligand)})
    result = run_grid_stage(target, run / "grid" / "production", box,
                            types, resolve_executable("autogrid4", "AutoGrid4"))
    _save_json(run / "protocol.json", {
        "run_id": token, "settings": t,
        "receptor_sha256": _sha(root(project) / "receptor.pdbqt"),
        "ligand_sha256": {x.name: _sha(x) for x in ligands},
        "metal_atoms_original": _input_metals(project),
        "metal_atoms_prepared": _pdbqt_metals(project),
        "validation": "NOT_PERFORMED",
        "method": "MANUAL_AUTODOCK4",
    })
    return {"run_id": token, "status": "CREATED", "grid": str(result["gpf"])}


def dock(project: Path) -> dict:
    project = Path(project).resolve()
    token, ligands = _run_identity(project)
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
    pending = [p for p in ligands if not (dlg_dir / (p.stem + ".dlg")).is_file()
               or (dlg_dir / (p.stem + ".dlg")).stat().st_size == 0]
    errors = []
    def one(ligand):
        return run_compound_docking(ligand, "receptor", work_dir, dpf_dir,
                                     dlg_dir, pdbqt_atom_types(ligand),
                                     params, resolve_executable("autodock4", "AutoDock4"))
    with ThreadPoolExecutor(max_workers=min(len(pending) or 1, t["parallel_jobs"])) as pool:
        futures = {pool.submit(one, ligand): ligand for ligand in pending}
        for future in as_completed(futures):
            ligand = futures[future]
            try:
                future.result()
            except Exception as exc:
                errors.append({"ligand": ligand.stem, "error": str(exc)})
    result = {"run_id": token, "total": len(ligands),
              "already_finished": len(ligands) - len(pending),
              "launched": len(pending),
              "completed": sum((dlg_dir / (ligand.stem + ".dlg")).is_file() and
                               (dlg_dir / (ligand.stem + ".dlg")).stat().st_size > 0 for ligand in ligands),
              "errors": errors}
    _save_json(run / "docking_result.json", result)
    if errors:
        raise RuntimeError(f"{len(errors)} AutoDock jobs failed; see {run/'docking_result.json'}")
    return result


def job_path(project: Path) -> Path:
    return root(project) / "job.json"


def job_status(project: Path) -> dict:
    job = _read_json(job_path(project))
    if not job:
        return {}
    if job.get("status") == "RUNNING":
        pid = job.get("pid")
        try:
            os.kill(int(pid), 0)
        except (ValueError, TypeError, ProcessLookupError):
            job["status"] = "EXITED"
        except PermissionError:
            pass
    return job


def launch(project: Path, action: str) -> dict:
    if action not in {"prepare", "grid", "dock"}:
        raise ValueError("Unsupported manual action")
    if job_status(project).get("status") == "RUNNING":
        raise RuntimeError("A manual SIDock action is already running")
    folder = root(project)
    folder.mkdir(parents=True, exist_ok=True)
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
        results = {"prepare": prepare, "grid": grid, "dock": dock}[action](project)
        job = _read_json(path)
        job.update({"status": "SUCCEEDED", "result": results,
                    "finished": datetime.now(timezone.utc).isoformat()})
        _save_json(path, job)
        return results
    except Exception as exc:
        job = _read_json(path)
        job.update({"status": "FAILED", "error": str(exc),
                    "finished": datetime.now(timezone.utc).isoformat()})
        _save_json(path, job)
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description="Explicit manual AutoDock4 runner")
    parser.add_argument("project", type=Path)
    parser.add_argument("action", choices=("prepare", "grid", "dock", "status"))
    args = parser.parse_args()
    if args.action == "status":
        print(json.dumps({"settings": settings(args.project),
                          "job": job_status(args.project)}, indent=2))
    else:
        print(json.dumps(execute(args.project, args.action), indent=2))


if __name__ == "__main__":
    main()
