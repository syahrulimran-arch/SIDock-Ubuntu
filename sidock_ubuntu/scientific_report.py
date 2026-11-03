"""Deterministic, source-grounded interpretation of a completed SIDock study.

Reads finalized evidence and persisted AutoDock/PLIP records. Never reruns
docking or assigns unrecorded residue interactions or biological activity.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import re
from collections import Counter
from pathlib import Path


def _json(path):
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def _csv(path):
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))
    except (OSError, csv.Error):
        return []


def _num(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def _integer(value):
    number = _num(value)
    return int(number) if number is not None and number >= 0 and number.is_integer() else None


def _fmt(value, places=2, missing="NA"):
    number = _num(value)
    return missing if number is None else f"{number:.{places}f}"


def _settings(path):
    keys = {"ga_run", "ga_pop_size", "ga_num_evals", "ga_num_generations", "rmstol", "seed"}
    out = {}
    try:
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            fields = line.split()
            if len(fields) > 1 and fields[0] in keys:
                out[fields[0]] = " ".join(fields[1:])
    except OSError:
        pass
    return out


def _residue(row):
    name = (row.get("protein_residue_name") or "").strip()
    number = (row.get("protein_residue_number") or "").strip()
    chain = (row.get("protein_chain") or "").strip()
    return f"{name}{number}" + (f" chain {chain}" if chain else "") if name and number else None


def _activity(project):
    values = {}
    for row in _csv(project / "analysis/assay.csv"):
        cid = (row.get("ligand") or "").strip()
        value = _num(row.get("ic50_uM"))
        if cid and value is not None and value > 0:
            values[cid] = value
    return values


def build_summary(project):
    project = Path(project).resolve()
    evidence_path = project / "analysis/evidence.json"
    if not evidence_path.is_file():
        raise FileNotFoundError(f"Finalized SIDock evidence missing: {evidence_path}")
    evidence = _json(evidence_path)
    activity = _activity(project)
    compounds = []
    for item in evidence.get("compounds", []):
        cid = str(item.get("compound_id", "")).strip()
        if not cid:
            continue
        plip_path = project / "analysis" / f"interactions_plip_{cid}.csv"
        interactions = _csv(plip_path)
        grouped = {}
        for row in interactions:
            kind = (row.get("interaction_type") or "").strip()
            residue = _residue(row)
            if kind and residue:
                grouped.setdefault(kind, Counter())[residue] += 1
        clusters = _csv(project / "analysis" / f"clusters_{cid}.csv")
        if not clusters:
            clusters = [x for x in _csv(project / "analysis/clusters.csv") if x.get("ligand") == cid]
        populations = [_integer(x.get("cluster_size")) for x in clusters]
        populations = [x for x in populations if x is not None]
        largest = max(populations) if populations else _integer(item.get("largest_cluster_population"))
        counted = sum(populations) if populations else None
        counts = {"hydrogen_bond": item.get("PLIP_H_bonds"), "hydrophobic": item.get("PLIP_hydrophobic"),
                  "pi-stack": item.get("PLIP_pi_stack"), "salt_bridge": item.get("PLIP_salt_bridge")}
        # When the actual PLIP CSV is absent, do not label summary counts as confirmed contacts.
        counts = {key: _integer(value) if plip_path.is_file() else None for key, value in counts.items()}
        compounds.append({
            "compound_id": cid, "status": item.get("status"),
            "energy_kcal_mol": _num(item.get("lowest_docking_energy")),
            "cluster_count": _integer(item.get("cluster_count")),
            "largest_cluster_population": largest,
            "counted_cluster_poses": counted,
            "largest_cluster_percentage": 100 * largest / counted if counted and largest is not None else None,
            "plip_available": plip_path.is_file(),
            "plip_counts": counts,
            "residues_by_interaction": {key: [{"residue": residue, "count": n} for residue, n in counter.most_common()]
                                        for key, counter in sorted(grouped.items())},
            "ic50_uM": activity.get(cid),
            "dpf": _settings(project / "docking/dpf" / f"{cid}.dpf"),
            "dlg_available": (project / "docking/dlg" / f"{cid}.dlg").is_file(),
            "selected_pose": item.get("selected_pose"),
        })
    ranked = sorted([r for r in compounds if r["status"] == "COMPLETE" and r["energy_kcal_mol"] is not None],
                    key=lambda x: (x["energy_kcal_mol"], x["compound_id"]))
    stats = {"n": 0, "pearson_r": None, "spearman_rho": None}
    paired = [r for r in ranked if r["ic50_uM"] is not None]
    stats["n"] = len(paired)
    if len(paired) >= 3:
        try:
            from scipy.stats import pearsonr, spearmanr
            x = [r["energy_kcal_mol"] for r in paired]
            y = [-math.log10(r["ic50_uM"] * 1e-6) for r in paired]
            if len(set(x)) > 1 and len(set(y)) > 1:
                stats["pearson_r"] = float(pearsonr(x, y).statistic)
                stats["spearman_rho"] = float(spearmanr(x, y).statistic)
        except (ImportError, ValueError, TypeError):
            pass
    return {
        "schema_version": "0.6-scientific-report-1",
        "evidence_sha256": hashlib.sha256(evidence_path.read_bytes()).hexdigest(),
        "project_status": evidence.get("project_status"),
        "target": evidence.get("target") or {},
        "validation": evidence.get("validation") or {},
        "activity_correlation": stats,
        "ranked_ids": [r["compound_id"] for r in ranked],
        "compounds": compounds,
    }


def draft_markdown(summary):
    compounds = summary["compounds"]
    by_id = {r["compound_id"]: r for r in compounds}
    ranked = [by_id[cid] for cid in summary["ranked_ids"]]
    target = summary["target"]
    validation = summary["validation"]
    params = ranked[0]["dpf"] if ranked else {}
    lines = [
        "# SIDock scientific analysis and draft Results and Discussion", "",
        "AUTOMATED EVIDENCE-BOUND DRAFT. Verify target identity, structures, pose geometry, interaction assignments and methods before publication.", "",
        "## Target, site and validation", "",
        "- Receptor source label: " + str(target.get("pdb") or "not recorded") + " (a local filename is not a verified enzyme name or PDB accession).",
        "- Crystallographic reference instance: " + str(target.get("reference_ligand") or "not recorded") + ".",
        "- Native redocking: " + str(validation.get("sampling_result") or "not available") +
        "; minimum RMSD " + _fmt(validation.get("native_redocking_RMSD_A"), 3) + " Angstrom; threshold " +
        _fmt(validation.get("rmsd_threshold")) + " Angstrom.",
        "- Validation policy: " + str(validation.get("validation_policy") or "not recorded") + ".",
        "- Project scientific status: " + str(summary["project_status"] or "not recorded") + ".", "",
        "## Methods, from persisted inputs", "",
        "AutoDock4 predicted docking poses and energies; the lowest reported predicted docking energy per compound is tabulated. "
        "The genetic-algorithm settings below were read from a production DPF, and residue contacts from PLIP CSVs. "
        "The original input, DPF, DLG and PLIP files remain the scientific primary record.", "",
    ]
    for label, key in [("GA runs per compound", "ga_run"), ("GA population", "ga_pop_size"),
                       ("Maximum evaluations", "ga_num_evals"), ("Maximum generations", "ga_num_generations"),
                       ("Clustering RMSD tolerance", "rmstol"), ("Seeds", "seed")]:
        if key in params:
            lines.append("- " + label + ": " + params[key] + ".")
    lines.extend(["", "## Docking ranking and PLIP overview", "",
                  "| Compound | Energy (kcal/mol) | Clusters | Largest cluster | H-bonds | Hydrophobic | IC50 (uM) |",
                  "|---|---:|---:|---:|---:|---:|---:|"])
    for row in ranked:
        def n(value):
            return "NA" if value is None else str(value)
        lines.append("| " + " | ".join([
            row["compound_id"], _fmt(row["energy_kcal_mol"]), n(row["cluster_count"]),
            n(row["largest_cluster_population"]), n(row["plip_counts"]["hydrogen_bond"]),
            n(row["plip_counts"]["hydrophobic"]), _fmt(row["ic50_uM"], 3)
        ]) + " |")
    lines.extend(["", "## Draft Results and Discussion", ""])
    if ranked:
        top = ranked[0]
        sentence = (f"{len(ranked)} of {len(compounds)} compounds had completed docking and numeric scores. "
                    f"{top['compound_id']} showed the most favorable AutoDock4 score "
                    f"({_fmt(top['energy_kcal_mol'])} kcal/mol).")
        if len(ranked) > 1:
            second = ranked[1]
            sentence += (f" The next-ranked compound was {second['compound_id']} "
                         f"({_fmt(second['energy_kcal_mol'])} kcal/mol).")
        lines.append(sentence + " These are scoring-function estimates, not experimental binding free energies or confirmed inhibition.")
    else:
        lines.append("No complete compounds with numerical scores were recorded. No ranking is inferred.")
    lines.extend(["", "The cluster count describes the number of distinct pose groups in the docking search; "
                  "it does not on its own establish binding stability, convergence or potency. "
                  "Cluster populations and alternative top poses require inspection.", "",
                  "## Residue-level PLIP contact evidence", ""])
    for row in ranked:
        lines.extend(["### " + row["compound_id"], ""])
        if not row["plip_available"]:
            lines.append("PLIP CSV absent; interactions cannot be assigned from the summary alone.")
        elif not row["residues_by_interaction"]:
            lines.append("No protein residue identities were present in the recorded PLIP CSV.")
        else:
            for kind, entries in row["residues_by_interaction"].items():
                lines.append("- " + kind + ": " + ", ".join(x["residue"] + f" ({x['count']})" for x in entries) + ".")
        if row["counted_cluster_poses"]:
            lines.append("- Largest cluster population: " + str(row["largest_cluster_population"]) +
                         " of " + str(row["counted_cluster_poses"]) + " counted poses (" +
                         _fmt(row["largest_cluster_percentage"], 1) + "%).")
        lines.append("")
    stats = summary["activity_correlation"]
    lines.extend(["## Experimental activity comparison", ""])
    if stats["pearson_r"] is not None:
        lines.append(f"{stats['n']} compounds had numeric experimental IC50 (uM) values in the recorded assay table. "
                     f"Exploratory docking energy versus pIC50 Pearson r = {_fmt(stats['pearson_r'], 3)}; "
                     f"Spearman rho = {_fmt(stats['spearman_rho'], 3)}. Correlation is not predictive validation.")
    elif stats["n"]:
        lines.append(f"Only {stats['n']} compounds had paired docking and recorded IC50 values; "
                     "a meaningful correlation was not established.")
    else:
        lines.append("No numeric experimental IC50 data were present in analysis/assay.csv. "
                     "No experimental activity, inhibition, SAR or potency claims are made.")
    lines.extend(["", "## Limitations and next validation steps", "",
                  "Verify PDB provenance and exact target identity, biological site, ligand protonation and stereochemistry, "
                  "receptor/ligand coordinate frame, native validation mode and mapping, selected poses, "
                  "PLIP atom-level contacts/distances, and robustness to repeated seeds or receptor conformations.",
                  "Docking scores must not be converted directly into IC50 values.", "",
                  "## Audit sources", "",
                  "- analysis/evidence.json",
                  "- docking/dpf and docking/dlg: AutoDock settings and immutable original docking outputs",
                  "- analysis/clusters_<compound>.csv: clusters",
                  "- analysis/interactions_plip_<compound>.csv: residue-level contacts",
                  "- analysis/assay.csv: experimental values, when present",
                  "", "SHA256(evidence.json): " + summary["evidence_sha256"], ""])
    return "\n".join(lines)


def _html(markdown):
    parts = ["<!doctype html><html><head><meta charset='utf-8'><title>SIDock scientific report</title>"
             "<style>body{font:15px system-ui;line-height:1.5;max-width:1100px;margin:35px auto;color:#243047}"
             "table{border-collapse:collapse;width:100%}td,th{border:1px solid #ccd;padding:7px;text-align:left}"
             "th{background:#eef}h2{margin-top:28px}p{margin:10px 0}</style></head><body>"]
    table = False
    for line in markdown.splitlines():
        s = line.strip()
        if s.startswith("|") and s.endswith("|"):
            columns = [c.strip() for c in s.strip("|").split("|")]
            if all(re.fullmatch(r":?-{3,}:?", c or "") for c in columns):
                continue
            if not table:
                parts.append("<table>")
                table = True
            tag = "th" if columns[0] == "Compound" else "td"
            parts.append("<tr>" + "".join("<" + tag + ">" + html.escape(c) + "</" + tag + ">" for c in columns) + "</tr>")
            continue
        if table:
            parts.append("</table>")
            table = False
        if not s:
            continue
        if s.startswith("### "):
            parts.append("<h3>" + html.escape(s[4:]) + "</h3>")
        elif s.startswith("## "):
            parts.append("<h2>" + html.escape(s[3:]) + "</h2>")
        elif s.startswith("# "):
            parts.append("<h1>" + html.escape(s[2:]) + "</h1>")
        else:
            parts.append("<p>" + html.escape(s) + "</p>")
    if table:
        parts.append("</table>")
    return "\n".join(parts + ["</body></html>"])


def write_scientific_report(project):
    project = Path(project).resolve()
    summary = build_summary(project)
    markdown = draft_markdown(summary)
    out = project / "analysis"
    out.mkdir(parents=True, exist_ok=True)
    data = {"scientific_report.md": markdown, "scientific_report.html": _html(markdown),
            "scientific_report.json": json.dumps(summary, indent=2, sort_keys=True)}
    for name, content in data.items():
        path = out / name
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(content + "\n", encoding="utf-8")
        temporary.replace(path)
    return {"project": str(project), "files": [str(out / name) for name in data],
            "ranked_compounds": len(summary["ranked_ids"]), "docking_rerun": False}


def main():
    parser = argparse.ArgumentParser(description="Build SIDock analysis and draft write-up from existing artifacts, without docking.")
    parser.add_argument("project", type=Path)
    args = parser.parse_args()
    print(json.dumps(write_scientific_report(args.project), indent=2))


if __name__ == "__main__":
    main()
