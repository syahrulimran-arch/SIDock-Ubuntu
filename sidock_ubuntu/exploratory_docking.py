"""Hash-bound technical AutoDock4 smoke tests for receptor candidates.

This module is deliberately separate from scientific docking approval. It runs
one small technical job only after validating the exploratory receptor and
parameter-library provenance.
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path

from .core import pdbqt_atom_types, resolve_executable
from .model import DockingParameters, GridBox
from .stages import run_compound_docking, run_grid_stage


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def _save_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")
    temp.replace(path)


def _atom_type_error(text: str) -> list[str]:
    patterns = (r"unknown atom type[^\n]*", r"atom type[^\n]*not found", r"ERROR[^\n]*atom", r"FATAL[^\n]*")
    found = []
    for pattern in patterns:
        found.extend(re.findall(pattern, text, flags=re.IGNORECASE))
    return sorted(set(found))


def _docked_ligand_atoms(dlg: Path) -> list[dict]:
    atoms = []
    for line in dlg.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("DOCKED: ATOM"):
            continue
        atom_line = line[8:]
        if len(atom_line) < 54:
            continue
        try:
            atoms.append({"xyz": tuple(float(atom_line[i:i + 8]) for i in (30, 38, 46)),
                          "type": atom_line[77:79].strip()})
        except (ValueError, IndexError):
            continue
    return atoms


def run_exploratory_smoke(project: Path) -> dict:
    """Run one small real AutoGrid4/AutoDock4 job using existing artifacts."""
    from .manual_docking import root, settings
    from .openbabel_receptor import prepare_exploratory_receptor
    from .core import run_command

    project = Path(project).resolve()
    home = root(project)
    receptor_info = prepare_exploratory_receptor(project)
    receptor = Path(receptor_info["candidate"])
    parameter_file = home / "custom_parameters.dat"
    ligand_source = home / "inputs/crystal_ligand.pdb"
    smoke = home / "exploratory_smoke"
    ligand_dir = smoke / "ligand"
    ligand = ligand_dir / "crystal_ligand.pdbqt"
    ligand_log = smoke / "openbabel_ligand.log"
    report_path = smoke / "exploratory_smoke_audit.json"
    if not ligand_source.is_file():
        raise RuntimeError("No saved crystallographic ligand is available for the technical smoke test")
    if not parameter_file.is_file():
        raise RuntimeError("custom_parameters.dat is required for both AutoGrid4 and AutoDock4")
    ligand_dir.mkdir(parents=True, exist_ok=True)
    if not ligand.is_file():
        obabel = shutil.which("obabel")
        if not obabel:
            raise RuntimeError("Open Babel is required to prepare the smoke-test ligand")
        code = run_command([obabel, str(ligand_source.resolve()), "-O", str(ligand.resolve()),
                            "-h", "-xn", "-xp", "--partialcharge", "gasteiger"],
                           ligand_dir, ligand_log)
        if code != 0 or not ligand.is_file() or ligand.stat().st_size == 0:
            raise RuntimeError(f"Open Babel ligand conversion failed; see {ligand_log}")
    ligand_types = sorted(set(pdbqt_atom_types(ligand)))
    receptor_types = sorted(set(pdbqt_atom_types(receptor)))
    declared = set()
    for line in parameter_file.read_text(encoding="utf-8", errors="strict").splitlines():
        parts = line.split()
        if parts and parts[0] == "atom_par" and len(parts) >= 12:
            declared.add(parts[1])
    missing = sorted((set(receptor_types) | set(ligand_types)) - declared)
    if missing:
        raise RuntimeError("custom_parameters.dat lacks receptor/ligand types: " + ", ".join(missing))

    identity = {
        "receptor_sha256": _sha(receptor), "parameter_sha256": _sha(parameter_file),
        "ligand_sha256": _sha(ligand), "ligand_source_sha256": _sha(ligand_source),
        "receptor_types": receptor_types, "ligand_types": ligand_types,
        "mode": "TECHNICAL_EXPLORATORY_ONLY",
    }
    old = _json(report_path)
    if old and old.get("identity") == identity and old.get("status") in {"PASSED", "FAILED"}:
        old["technical_exploratory_approval"] = {
            "status": "APPROVED_FOR_TECHNICAL_SMOKE_ONLY",
            "receptor_sha256": identity["receptor_sha256"],
            "parameter_sha256": identity["parameter_sha256"],
            "scientific_approval": False,
        }
        _save_json(report_path, old)
        return old
    grid_dir = smoke / "grid" / "production"
    docking_dir = smoke / "docking"
    grid = GridBox(settings(project)["center_x"], settings(project)["center_y"], settings(project)["center_z"], 20, 20, 20, 0.375)
    params = DockingParameters(ga_runs=1, ga_num_evals=1000, ga_pop_size=20,
                               ga_num_generations=100, parallel_jobs=1, random_seed=20250301)
    grid_result = run_grid_stage(receptor, grid_dir, grid, ligand_types,
                                 resolve_executable("autogrid4", "AutoGrid4"),
                                 parameter_file=parameter_file)
    dock_result = run_compound_docking(ligand, "receptor", docking_dir / "work",
                                       docking_dir / "dpf", docking_dir / "dlg", ligand_types,
                                       params, resolve_executable("autodock4", "AutoDock4"),
                                       parameter_file=parameter_file)
    glg = Path(grid_result["glg"])
    dlg = Path(dock_result["dlg"])
    glg_text = glg.read_text(encoding="utf-8", errors="replace")
    dlg_text = dlg.read_text(encoding="utf-8", errors="replace")
    maps = [Path(x) for x in grid_result["maps"]]
    grid_parameter_copy = grid_dir / "custom_ad4_parameters.dat"
    dock_parameter_copy = dock_result["work_dir"] / "custom_ad4_parameters.dat"
    errors = _atom_type_error(glg_text + "\n" + dlg_text)
    docked = _docked_ligand_atoms(dlg)
    ni = []
    for line in receptor.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and line[77:79].strip() == "Ni":
            ni.append(tuple(float(line[i:i + 8]) for i in (30, 38, 46)))
    import math
    ni_distances = [round(min(math.dist(n, a["xyz"]) for a in docked), 3)
                    for n in ni] if docked else []
    report = {
        "status": "PASSED" if ("Successful Completion" in glg_text and "Successful Completion" in dlg_text
                                  and not errors and all(x.is_file() for x in maps)
                                  and _sha(grid_parameter_copy) == identity["parameter_sha256"]
                                  and _sha(dock_parameter_copy) == identity["parameter_sha256"]) else "FAILED",
        "identity": identity, "technical_exploratory_only": True,
        "scientific_approval": False, "parameter_file": str(parameter_file),
        "technical_exploratory_approval": {
            "status": "APPROVED_FOR_TECHNICAL_SMOKE_ONLY",
            "receptor_sha256": identity["receptor_sha256"],
            "parameter_sha256": identity["parameter_sha256"],
            "scientific_approval": False,
        },
        "grid": {"gpf": str(grid_result["gpf"]), "glg": str(glg), "maps": [str(x) for x in maps],
                 "maps_exist": all(x.is_file() for x in maps), "parameter_file_line": "parameter_file custom_ad4_parameters.dat"},
        "docking": {"dpf": str(dock_result["dpf"]), "dlg": str(dlg),
                    "successful_completion": "Successful Completion" in dlg_text,
                    "parameter_file_line": "parameter_file custom_ad4_parameters.dat"},
        "atom_type_errors": errors, "ligand_types": ligand_types,
        "ni_nearest_ligand_distances_angstrom": ni_distances,
        "geometry_review": "TECHNICAL_DISTANCE_SCREEN_ONLY",
        "logs": {"ligand": str(ligand_log), "autogrid": str(grid_dir / "autogrid.log"),
                 "autodock": str(dock_result["work_dir"] / (ligand.stem + ".log"))},
    }
    _save_json(report_path, report)
    return report
