from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import copy
import shutil
try:
    import yaml
except ImportError:  # keep config commands usable before optional environment installation
    yaml = None

VALIDATION_MODES = ("required", "if_available", "skip")

DEFAULT_CONFIG = {
    "reference": {"ligand": None, "remote_enrichment": False},
    "validation": {"mode": "required", "rmsd_threshold": 2.0, "ga_runs": 10},
    "docking": {"ga_runs": 10, "parallel_jobs": "auto", "population_size": 150,
                 "evaluations": 2500000, "generations": 27000,
                 "mutation_rate": 0.02, "crossover_rate": 0.8, "local_search_frequency": 0.26},
    "interaction": {"primary_engine": "plip"},
    "output": {"figures": True, "retain_raw_dlg": True},
}

def deep_merge(base: dict, override: dict) -> dict:
    result = copy.deepcopy(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = value
    return result

def load_config(path: Path | None = None) -> dict:
    if path is None:
        path = Path("sidock.yaml")
    if not path.exists():
        return copy.deepcopy(DEFAULT_CONFIG)
    if yaml is None:
        raise RuntimeError("PyYAML is required to read sidock.yaml; run ./install.sh")
    with path.open(encoding="utf-8") as fh:
        return deep_merge(DEFAULT_CONFIG, yaml.safe_load(fh) or {})

def write_default_config(path: Path = Path("sidock.yaml")) -> None:
    if yaml is None:
        raise RuntimeError("PyYAML is required to write sidock.yaml; run ./install.sh")
    path.write_text(yaml.safe_dump(DEFAULT_CONFIG, sort_keys=False), encoding="utf-8")

def validate_config(config: dict) -> list[str]:
    errors = []
    mode=str(config["validation"].get("mode", "required")).lower()
    if mode not in VALIDATION_MODES:
        errors.append("validation.mode must be one of: required, if_available, skip")
    if config["validation"]["rmsd_threshold"] <= 0: errors.append("validation.rmsd_threshold must be positive")
    if config["validation"]["ga_runs"] < 1: errors.append("validation.ga_runs must be >= 1")
    if int(config["docking"].get("ga_runs", 0)) < 1: errors.append("docking.ga_runs must be >= 1")
    parallel=config["docking"].get("parallel_jobs", "auto")
    if str(parallel).lower() != "auto":
        try:
            if int(parallel) < 1: errors.append("docking.parallel_jobs must be auto or >= 1")
        except (TypeError, ValueError):
            errors.append("docking.parallel_jobs must be auto or >= 1")
    for key in ("population_size", "evaluations", "generations"):
        if config["docking"][key] < 1: errors.append(f"docking.{key} must be >= 1")
    if config["interaction"]["primary_engine"].lower() != "plip": errors.append("interaction.primary_engine must be plip")
    return errors

def summary(config: dict) -> str:
    d = config["docking"]
    return ("SIDock Ubuntu v0.6-dev\n"
            f"Validation mode: {config['validation'].get('mode','required')}\n"
            f"Validation RMSD: <={config['validation']['rmsd_threshold']:.1f} Å\n"
            f"Validation GA runs: {config['validation']['ga_runs']}\n"
            f"Docking GA runs per compound: {d.get('ga_runs',10)}\n"
            f"Parallel compound jobs: {d.get('parallel_jobs','auto')}\n"
            f"Population: {d['population_size']}\n"
            f"Evaluations: {d['evaluations']:,}\n"
            f"Interaction engine: {config['interaction']['primary_engine'].upper()}")