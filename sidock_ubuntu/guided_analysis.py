"""Run-scoped post-processing of saved AutoDock4 artifacts; no redocking.

Each attempt has a separate revision directory. A report can only be reused
when all recorded inputs and outputs still match their SHA256 hashes.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shutil
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from .analysis import extract_docking_poses, analyze_autodock_clusters
from .assay import extract_sdf_assay
from .dependencies import resolve_plip_executable
from .guided_analysis_figures import generate_figures
from .manual_docking import root, active_run, complete_dlg_energy, _ligand_sdf, _save_json, _read_json
from .plip_runner import run_plip_analysis
from .scientific_report import write_scientific_report

SCHEMA = "guided-analysis-1"


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _csv(path, fields, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def resolve_run(project, run_id=None):
    project = Path(project).resolve()
    if run_id is not None:
        if not re.fullmatch(r"[0-9a-f]{16}", run_id):
            raise ValueError("Expected a saved 16-character docking run ID")
        run = root(project) / "runs" / run_id
    else:
        run = active_run(project)
    if run is None or not (run / "protocol.json").is_file():
        raise RuntimeError("No saved docking run for the selected settings. Generate the grid and complete docking first.")
    if not run.resolve().is_relative_to((root(project) / "runs").resolve()):
        raise ValueError("Analysis run must be inside this project")
    return run


def analysis_status(project, run):
    """Read saved status and detect edited/missing inputs/outputs; no analysis."""
    result = _read_json(run / "analysis_result.json")
    if not result:
        return {}
    base = Path(project).resolve()
    for group in ("inputs", "artifacts"):
        for name, expected in result.get(group, {}).items():
            path = (base / name).resolve()
            if not path.is_relative_to(base):
                return {**result, "status": "STALE", "stale_reason": "Artifact is outside the project"}
            actual = _sha(path) if path.is_file() else None
            if actual != expected:
                return {**result, "status": "STALE", "stale_reason": "Changed or missing artifact: " + name}
    return result


def _inputs(project, run):
    protocol = _read_json(run / "protocol.json")
    if protocol.get("run_id") != run.name:
        raise RuntimeError("Saved protocol does not match this docking run")
    receptor = run / "grid/production/receptor.pdbqt"
    if not receptor.is_file() or _sha(receptor) != protocol.get("receptor_sha256"):
        raise RuntimeError("Saved run receptor differs from its protocol; refusing mixed-run analysis")
    files = {"protocol": run / "protocol.json", "receptor": receptor}
    ligand_sources = {}
    names = protocol.get("ligand_sha256", {})
    if not names:
        raise RuntimeError("Saved protocol has no prepared-ligand hashes")
    for name, expected in sorted(names.items()):
        if Path(name).name != name or not re.fullmatch(r"[A-Za-z0-9_.-]+\.pdbqt", name):
            raise ValueError("Invalid ligand name in saved protocol")
        cid = Path(name).stem
        staged = run / "docking/work" / cid / name
        prepared = root(project) / "ligands/prepared" / name
        # Prefer the frozen per-compound copy used by AutoDock.
        source = next((p for p in (staged, prepared) if p.is_file() and _sha(p) == expected), None)
        ligand_sources[cid] = source
        if source is not None:
            files["ligand:" + cid] = source
        files["dlg:" + cid] = run / "docking/dlg" / (cid + ".dlg")
        files["dpf:" + cid] = run / "docking/dpf" / (cid + ".dpf")
        minimized = root(project) / "ligands/source" / (cid + ".sdf")
        if minimized.is_file():
            files["sdf:" + cid] = minimized
    try:
        sdf = _ligand_sdf(project)
    except FileNotFoundError:
        sdf = None
    manifest = root(project) / "ligands/manifest.json"
    linked_assay = None
    if sdf is not None and sdf.is_file() and manifest.is_file():
        data = _read_json(manifest)
        if data.get("source_sha256") == _sha(sdf) and set(data.get("compounds", [])) == set(names):
            current = root(project) / "ligands/prepared"
            if all((current / name).is_file() and _sha(current / name) == expected for name, expected in names.items()):
                linked_assay = sdf
                files["assay_source"] = sdf
                files["ligand_manifest"] = manifest
    # All downloadable artifacts and provenance are confined to this project.
    for path in files.values():
        if not path.resolve().is_relative_to(project):
            raise RuntimeError("Analysis input is outside this project: " + str(path))
    inputs = {p.relative_to(project).as_posix(): _sha(p) if p.is_file() else None for p in files.values()}
    return protocol, receptor, ligand_sources, linked_assay, inputs


def analyze(project, *, run_id=None):
    project = Path(project).resolve()
    run = resolve_run(project, run_id)
    protocol, receptor, ligand_sources, sdf, inputs = _inputs(project, run)
    if not any(complete_dlg_energy(run / "docking/dlg" / (cid + ".dlg")) is not None for cid in ligand_sources):
        raise RuntimeError("No successfully completed DLG files to analyse. Docking is not rerun by analysis.")
    plip = resolve_plip_executable()
    identity = {"schema": SCHEMA, "inputs": inputs,
                "plip_available": plip["available"], "plip_version": plip.get("version"),
                "pymol_available": bool(shutil.which("pymol"))}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    previous = analysis_status(project, run)
    if (previous.get("fingerprint") == fingerprint and not previous.get("retryable")
            and previous.get("status") in {"COMPLETE", "COMPLETE_WITH_WARNINGS"}):
        return {**previous, "reused": True, "docking_rerun": False}
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = run / "analysis/revisions" / (fingerprint[:12] + "-" + stamp)
    out.mkdir(parents=True)
    warnings, rows = [], []
    if not plip["available"]:
        warnings.append("PLIP unavailable: " + plip["message"] + " Install/repair PLIP and rerun analysis; docking does not need to be repeated.")
    if sdf is None:
        warnings.append("No hash-linked input SDF is available for this saved run; experimental IC50 comparisons are unavailable.")
    for cid, ligand_source in ligand_sources.items():
        row = {"compound_id": cid, "status": "FAILED", "analysis_status": "FAILED",
               "lowest_docking_energy": None, "selected_pose": None,
               "cluster_count": None, "largest_cluster_population": None,
               "PLIP_H_bonds": None, "PLIP_hydrophobic": None,
               "PLIP_pi_stack": None, "PLIP_salt_bridge": None,
               "plip_status": "UNAVAILABLE", "interactions": [], "errors": [], "warnings": []}
        dlg = run / "docking/dlg" / (cid + ".dlg")
        try:
            energy = complete_dlg_energy(dlg)
            if energy is None or not math.isfinite(energy):
                raise RuntimeError("DLG missing, incomplete or has no finite docking score")
            if ligand_source is None:
                raise RuntimeError("Prepared ligand copy does not match the saved protocol")
            pose = extract_docking_poses(dlg, out / "poses", out / "complexes", receptor)
            if not pose.selected_pose_path or pose.lowest_energy is None or not math.isfinite(pose.lowest_energy):
                raise RuntimeError("Completed DLG contains no usable energy-ranked pose")
            if abs(pose.lowest_energy - energy) > .011:
                raise RuntimeError("DLG summary score and extracted best-pose score disagree")
            # Confirm that the selected model preserves the docked atom identities.
            def identities(path):
                return Counter((line[6:11].strip(), line[12:16].strip(), line[77:79].strip())
                               for line in path.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM")))
            if identities(ligand_source) != identities(pose.selected_pose_path):
                raise RuntimeError("Selected pose atom identities differ from the prepared ligand")
            cluster = analyze_autodock_clusters(dlg, out / ("clusters_" + cid + ".csv"))
            row.update(status="COMPLETE", analysis_status="PARTIAL", lowest_docking_energy=pose.lowest_energy,
                       selected_pose=str(pose.selected_pose_path), selected_model=pose.selected_model,
                       complex_file=str(next(p["complex"] for p in pose.poses if p["model"] == pose.selected_model)),
                       cluster_count=cluster.cluster_count if cluster.clusters else None,
                       largest_cluster_population=cluster.largest_cluster_population)
            if not cluster.clusters:
                row["warnings"].append("No cluster-population table recorded in DLG; cluster statistics are unavailable")
            if plip["available"]:
                try:
                    work = out / ("plip_" + cid)
                    contacts = out / ("interactions_plip_" + cid + ".csv")
                    run_plip_analysis(Path(row["complex_file"]), work, contacts, compound_id=cid,
                                      pose_id=str(pose.selected_model), ligand_identity=("LIG", "Z", "999"))
                    with contacts.open(newline="") as handle:
                        interactions = list(csv.DictReader(handle))
                    count = Counter(r["interaction_type"] for r in interactions)
                    row.update(plip_status="COMPLETE", analysis_status="COMPLETE", interactions=interactions,
                               PLIP_H_bonds=count["hydrogen_bond"], PLIP_hydrophobic=count["hydrophobic"],
                               PLIP_pi_stack=count["pi-stack"], PLIP_salt_bridge=count["salt_bridge"])
                    _save_json(work / "provenance.json", {"ligand_identity": ["LIG", "Z", "999"],
                                                          "dlg_sha256": _sha(dlg), "receptor_sha256": _sha(receptor),
                                                          "selected_pose_sha256": _sha(pose.selected_pose_path),
                                                          "scientific_approval": False})
                except Exception as exc:
                    row["errors"].append("PLIP failed: " + str(exc))
            else:
                row["errors"].append("PLIP is unavailable; residue interactions are not assigned")
        except Exception as exc:
            row["errors"].append(str(exc))
        rows.append(row)
    ranked = sorted([r for r in rows if r["status"] == "COMPLETE"], key=lambda r: (r["lowest_docking_energy"], r["compound_id"]))
    _csv(out / "ranking.csv", ["rank", "ligand", "energy_kcal_mol", "pose_file", "complex_file"],
         [{"rank": i, "ligand": r["compound_id"], "energy_kcal_mol": r["lowest_docking_energy"],
           "pose_file": r["selected_pose"], "complex_file": r["complex_file"]} for i, r in enumerate(ranked, 1)])
    activity_rows, structures = [], {}
    if sdf is not None:
        try:
            assay = extract_sdf_assay(sdf, out / "assay.csv")
            from rdkit import Chem
            source_molecules = list(Chem.SDMolSupplier(str(sdf), removeHs=False))
            # Match the compound naming rule used during ligand preparation.
            for row in assay:
                mol = source_molecules[int(row["input_record_index"]) - 1]
                name = mol.GetProp("Compound_ID") if mol.HasProp("Compound_ID") else mol.GetProp("_Name") if mol.HasProp("_Name") else ""
                row["ligand"] = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
            _csv(out / "assay.csv", ["ligand", "ic50_uM", "pIC50", "input_record_index", "final_status", "failure_reason"], assay)
            values = {r["ligand"]: r for r in assay if r["pIC50"] is not None}
            activity_rows = [{"ligand": r["compound_id"], "energy_kcal_mol": r["lowest_docking_energy"],
                              "pIC50": values[r["compound_id"]]["pIC50"], "ic50_uM": values[r["compound_id"]]["ic50_uM"]}
                             for r in ranked if r["compound_id"] in values]
            for mol in source_molecules:
                if mol is not None:
                    name = mol.GetProp("Compound_ID") if mol.HasProp("Compound_ID") else mol.GetProp("_Name") if mol.HasProp("_Name") else ""
                    cid = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
                    if cid in structures:
                        raise RuntimeError("Ambiguous compound IDs in source SDF")
                    structures[cid] = mol
            warnings.append("Bare numeric IC50 values use the input SDF's micromolar schema. Qualified values and unsupported units are excluded from correlations.")
        except Exception as exc:
            warnings.append("Experimental activity extraction unavailable: " + str(exc))
    figures, retryable_figures = [], False
    try:
        figures, figure_warnings = generate_figures(rows, receptor, out / "figures", activity_rows=activity_rows, structures=structures)
        warnings.extend(figure_warnings)
        retryable_figures = any(not warning.startswith("PyMOL is unavailable:") for warning in figure_warnings)
    except Exception as exc:
        retryable_figures = True
        warnings.append("Figure generation incomplete: " + str(exc))
    model = protocol.get("model_status") or "MANUAL_AUTODOCK4_UNVALIDATED"
    status = "COMPLETE" if rows and all(r["analysis_status"] == "COMPLETE" for r in rows) else "PARTIAL"
    if not ranked:
        status = "FAILED"
    if status == "COMPLETE" and (warnings or any(r["warnings"] for r in rows)):
        status = "COMPLETE_WITH_WARNINGS"
    columns = ["compound_id", "status", "analysis_status", "lowest_docking_energy", "cluster_count",
               "largest_cluster_population", "PLIP_H_bonds", "PLIP_hydrophobic", "PLIP_pi_stack",
               "PLIP_salt_bridge", "selected_pose", "plip_status", "model_status", "scientific_approval", "errors"]
    _csv(out / "results.csv", columns, [{**r, "model_status": model, "scientific_approval": False,
                                        "errors": "; ".join(r["errors"])} for r in rows])
    _csv(out / "activity_comparison.csv", ["ligand", "energy_kcal_mol", "ic50_uM", "pIC50"], activity_rows)
    for row in rows:
        warnings.extend(row["compound_id"] + ": " + item for item in row["warnings"] + row["errors"])
    evidence = {"schema_version": SCHEMA, "project_status": status, "analysis_status": status,
                "run_id": run.name, "model_status": model, "scientific_approval": False,
                "docking_rerun": False, "target": protocol.get("target") or {"pdb": "saved receptor.pdbqt"},
                "validation": {"validation_policy": "guided_exploratory", "sampling_result": "NOT_PERFORMED",
                               "docking_validation_status": "UNVALIDATED", "native_redocking_RMSD_A": None,
                               "scientific_approval": False},
                "compounds": rows, "figures": figures, "warnings": warnings,
                "provenance": identity}
    _save_json(out / "evidence.json", evidence)
    write_scientific_report(run, analysis_dir=out)
    _save_json(out / "provenance.json", {**identity, "model_status": model, "scientific_approval": False})
    # Archive only this revision and the original files it describes.
    bundle = out.parent / (out.name + ".zip")
    with zipfile.ZipFile(bundle, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(out.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(out).as_posix())
        archive.write(run / "protocol.json", "sources/protocol.json")
        archive.write(receptor, "sources/receptor.pdbqt")
        if sdf is not None:
            archive.write(sdf, "sources/compounds.sdf")
        for cid in ligand_sources:
            if ligand_sources[cid] is not None:
                archive.write(ligand_sources[cid], "sources/prepared/" + cid + ".pdbqt")
            for suffix, subdir in ((".dlg", "dlg"), (".dpf", "dpf")):
                path = run / "docking" / subdir / (cid + suffix)
                if path.is_file():
                    archive.write(path, "sources/" + subdir + "/" + path.name)
    artifacts = {p.relative_to(project).as_posix(): _sha(p) for p in out.rglob("*") if p.is_file()}
    artifacts[bundle.relative_to(project).as_posix()] = _sha(bundle)
    for name, expected in inputs.items():
        path = project / name
        if (_sha(path) if path.is_file() else None) != expected:
            raise RuntimeError("Input changed while analysing; retained this revision for inspection: " + name)
    result = {"schema": SCHEMA, "status": status, "run_id": run.name, "fingerprint": fingerprint,
              "retryable": retryable_figures or status in {"PARTIAL", "FAILED"},
              "analysis_directory": out.relative_to(project).as_posix(), "bundle": bundle.relative_to(project).as_posix(),
              "total": len(rows), "analysed": len(ranked), "plip_complete": sum(r["plip_status"] == "COMPLETE" for r in rows),
              "model_status": model, "scientific_approval": False, "docking_rerun": False,
              "warnings": warnings, "inputs": inputs, "artifacts": artifacts}
    _save_json(run / "analysis_result.json", result)
    return result
