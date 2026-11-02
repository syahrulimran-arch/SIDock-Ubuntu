from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import copy
import shutil
try:
    import yaml
except ImportError:  # keep config commands usable before optional environment installation
    yaml = None

DEFAULT_CONFIG = {
    "reference": {"ligand": None, "remote_enrichment": False},
    "validation": {"rmsd_threshold": 2.0, "ga_runs": 10},
    "docking": {"population_size": 150, "evaluations": 2500000, "generations": 27000,
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
    if config["validation"]["rmsd_threshold"] <= 0: errors.append("validation.rmsd_threshold must be positive")
    if config["validation"]["ga_runs"] < 1: errors.append("validation.ga_runs must be >= 1")
    for key in ("population_size", "evaluations", "generations"):
        if config["docking"][key] < 1: errors.append(f"docking.{key} must be >= 1")
    if config["interaction"]["primary_engine"].lower() != "plip": errors.append("interaction.primary_engine must be plip")
    return errors

def summary(config: dict) -> str:
    d = config["docking"]
    return ("SIDock Ubuntu v0.4.0\n"
            f"Validation RMSD: <={config['validation']['rmsd_threshold']:.1f} Å\n"
            f"GA runs: {config['validation']['ga_runs']}\n"
            f"Population: {d['population_size']}\n"
            f"Evaluations: {d['evaluations']:,}\n"
            f"Interaction engine: {config['interaction']['primary_engine'].upper()}")
