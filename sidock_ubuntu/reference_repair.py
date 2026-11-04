"""Repair an incorrectly selected crystallographic reference without docking.

This deliberately stops at REFERENCE. It does not prepare a receptor, grid,
native redocking, or compound docking. Upstream project inputs are immutable.
"""
from __future__ import annotations

import argparse
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

from .config import load_config
from .executors import production_context_from_project, build_production_executors
from .intake import eligible_reference_groups
from .orchestrator import StageRunner
from .project_lock import project_execution_lock
from .state import _atomic_write, invalidate_from, project_status


def repair_reference(project: Path, *, expected: str, apply: bool = False) -> dict:
    project = Path(project).resolve()
    data = json.loads((project / "project.json").read_text(encoding="utf-8"))
    receptor = Path(data["receptor_source"])
    if not receptor.is_file():
        raise FileNotFoundError(f"Receptor source missing: {receptor}")
    candidates = [g.key for g in eligible_reference_groups(receptor)]
    if expected not in candidates:
        raise ValueError(f"Reference {expected!r} is not an eligible nonpolymer crystallographic ligand; eligible={candidates}")
    state = project_status(project)
    if state.get("compounds") or any((project / "docking/dlg").glob("*.dlg")):
        raise RuntimeError("Reference correction refused: compound or DLG outputs exist; requires an explicit scientific reanalysis plan.")
    old_reference = project / "receptor/reference_ligand.pdb"
    old = ""
    if old_reference.is_file():
        for line in old_reference.read_text(encoding="utf-8",errors="replace").splitlines():
            if line.startswith(("ATOM", "HETATM")):
                old = f"{line[17:20].strip()}:{line[21:22].strip() or '-'}:{line[22:27].strip()}"
                break
    report = {"project": str(project), "current_reference": old or None, "proposed_reference": expected,
              "eligible_references": candidates, "compound_docking_rerun": False, "applied": False}
    if not apply:
        return report
    with project_execution_lock(project, {"interface": "CLI", "action": "REFERENCE_CORRECTION"}):
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        backup = project / "provenance" / ("reference_correction_backup_" + timestamp)
        backup.mkdir(parents=True, exist_ok=False)
        for relative in ("project_state.json", "receptor/reference_ligand.pdb",
                         "validation/reference_selection.json", "validation/reference_ligand_candidates.csv",
                         "validation/reference_stage.json"):
            src = project / relative
            if src.is_file():
                dest = backup / relative
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dest)

        current = project_status(project)
        _atomic_write(project / "project_state.json", invalidate_from(current, "REFERENCE",
                        "REFERENCE_LIGAND_CORRECTED_FROM_MODIFIED_POLYMER"))
        config = load_config()
        config["reference"] = {**config.get("reference", {}), "ligand": expected}
        ctx = production_context_from_project(project, config)
        result = StageRunner(project, config, build_production_executors(ctx)).run_global(("REFERENCE",))
        if not result or result[-1][1] != "PASS":
            raise RuntimeError(f"Reference repair failed: {result}; backup={backup}")
        ref = project / "receptor/reference_ligand.pdb"
        actual = next((f"{line[17:20].strip()}:{line[21:22].strip() or '-'}:{line[22:27].strip()}"
                       for line in ref.read_text().splitlines() if line.startswith(("ATOM", "HETATM"))), None)
        if actual != expected:
            raise RuntimeError(f"Post-repair reference mismatch: expected={expected}; found={actual}")
        report.update({"applied": True, "backup": str(backup), "current_reference": actual,
                       "receptor_preparation_performed": False, "metal_parameterization_performed": False})
        (project / "provenance/reference_correction_audit.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description="Correct saved reference ligand without running docking.")
    parser.add_argument("project", type=Path)
    parser.add_argument("--expected", required=True, help="Explicit component:chain:number, e.g. HAE:C:800")
    parser.add_argument("--apply", action="store_true", help="Perform correction; without this flag, inspect only")
    args = parser.parse_args()
    print(json.dumps(repair_reference(args.project, expected=args.expected, apply=args.apply), indent=2))


if __name__ == "__main__":
    main()
