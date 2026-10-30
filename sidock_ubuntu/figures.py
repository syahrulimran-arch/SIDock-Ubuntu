from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path


def activity_scatter(stats: dict[str, object], destination: Path) -> Path | None:
    rows = stats.get("rows", [])
    if not isinstance(rows, list) or len(rows) < 2:
        return None
    import matplotlib.pyplot as plt
    x = [float(r["energy_kcal_mol"]) for r in rows]
    y = [float(r["pIC50"]) for r in rows]
    labels = [str(r["ligand"]) for r in rows]
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.scatter(x, y)
    for a, b, label in zip(x, y, labels):
        ax.annotate(label, (a, b), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("AutoDock4 binding energy (kcal/mol)")
    ax.set_ylabel("pIC50")
    ax.set_title("Docking score versus experimental activity")
    fig.tight_layout(); destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=300); plt.close(fig)
    return destination


def pymol_binding_figures(analysis_dir: Path, figures_dir: Path, limit: int = 5) -> list[Path]:
    figures_dir.mkdir(parents=True, exist_ok=True)
    ranking = analysis_dir / "ranking.csv"
    if not ranking.exists():
        return []
    with ranking.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))[:limit]
    outputs: list[Path] = []
    pymol = shutil.which("pymol")
    for row in rows:
        rel = row.get("complex_file", "")
        if not rel:
            continue
        complex_path = analysis_dir / rel
        out = figures_dir / f"{row['ligand']}_binding_pose.png"
        script = figures_dir / f"{row['ligand']}_binding_pose.pml"
        script.write_text(
            f"load {complex_path}\nhide everything\nshow cartoon, polymer\nshow sticks, organic\nselect pocket, byres (polymer within 4.0 of organic)\nshow sticks, pocket\nbg_color white\nset ray_opaque_background, off\norient organic\nzoom organic, 8\nray 1800, 1400\npng {out}, dpi=300\nquit\n",
            encoding="utf-8",
        )
        if pymol:
            subprocess.run([pymol, "-cq", str(script)], check=False)
            if out.exists():
                outputs.append(out)
    return outputs
