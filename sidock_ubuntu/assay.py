from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

IC50_KEYS = ("IC50", "IC50_UM", "IC50_uM", "IC50 (uM)", "IC50_MICROMOLAR", "IC50 µM",
             "IC50_nM", "IC50 (nM)", "IC50_mM", "IC50 (mM)")
ID_KEYS = ("Compound_ID", "COMPOUND_ID", "ID", "Name", "NAME")


def _float(value: str | None) -> float | None:
    if value is None:
        return None
    text = str(value).strip().replace("µ", "u").replace("μ", "u")
    # Qualified ranges (e.g. >10), uncertainties and unrecognized units are
    # not exact measurements and must not enter a numerical correlation.
    m = re.fullmatch(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*(uM|nM|mM|M)?", text)
    if not m:
        return None
    value = float(m.group(1)) * {None: 1.0, 'uM': 1.0, 'nM': 0.001, 'mM': 1000.0, 'M': 1e6}[m.group(2)]
    return value if math.isfinite(value) and value > 0 else None


def extract_sdf_assay(sdf: Path, destination: Path) -> list[dict[str, object]]:
    from rdkit import Chem
    rows: list[dict[str, object]] = []
    supplier = Chem.SDMolSupplier(str(sdf), removeHs=False)
    for idx, mol in enumerate(supplier, start=1):
        if mol is None:
            continue
        props = {k: mol.GetProp(k) for k in mol.GetPropNames()}
        name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
        ligand = next((props[k].strip() for k in ID_KEYS if k in props and props[k].strip()), name)
        key = next((k for k in IC50_KEYS if k in props), None)
        raw = props.get(key) if key else None
        if raw and re.fullmatch(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", raw.strip()):
            units = next((props[k].strip() for k in ("IC50_unit", "IC50_units", "IC50_UNIT", "IC50_UNITS") if k in props), None)
            if units is None and key:
                units = "nM" if "nM" in key else "mM" if "mM" in key else None
            if units:
                raw = raw.strip() + " " + units
        ic50 = _float(raw)
        rows.append({"ligand": ligand, "ic50_uM": ic50, "pIC50": (-math.log10(ic50 * 1e-6) if ic50 and ic50 > 0 else None), "input_record_index": idx, "final_status": "ACCEPTED" if ligand else "REJECTED", "failure_reason": "" if ligand else "missing compound ID"})
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["ligand", "ic50_uM", "pIC50", "input_record_index", "final_status", "failure_reason"])
        writer.writeheader(); writer.writerows(rows)
    return rows


def correlate_with_ranking(assay_rows: list[dict[str, object]], ranking_csv: Path, destination: Path) -> dict[str, object]:
    with ranking_csv.open(encoding="utf-8") as fh:
        ranking = list(csv.DictReader(fh))
    assay = {str(r["ligand"]): r for r in assay_rows}
    rows = []
    for r in ranking:
        a = assay.get(str(r["ligand"]))
        if not a or a.get("pIC50") is None:
            continue
        rows.append({"ligand": r["ligand"], "energy_kcal_mol": float(r["energy_kcal_mol"]), "ic50_uM": float(a["ic50_uM"]), "pIC50": float(a["pIC50"])})
    stats: dict[str, object] = {"n": len(rows), "rows": rows, "pearson_energy_vs_pIC50": None, "spearman_energy_vs_pIC50": None, "r2_linear": None}
    if len(rows) >= 3:
        import scipy.stats as st
        import numpy as np
        x = np.array([r["energy_kcal_mol"] for r in rows], dtype=float)
        y = np.array([r["pIC50"] for r in rows], dtype=float)
        pearson = st.pearsonr(x, y); spearman = st.spearmanr(x, y)
        coef = np.polyfit(x, y, 1); pred = np.polyval(coef, x)
        ss_res = float(np.sum((y - pred) ** 2)); ss_tot = float(np.sum((y - y.mean()) ** 2))
        stats.update({"pearson_energy_vs_pIC50": float(pearson.statistic), "pearson_p": float(pearson.pvalue), "spearman_energy_vs_pIC50": float(spearman.statistic), "spearman_p": float(spearman.pvalue), "r2_linear": (1 - ss_res / ss_tot) if ss_tot else None})
    destination.write_text(json.dumps(stats, indent=2), encoding="utf-8")
    return stats
