from __future__ import annotations

import csv
import re
from pathlib import Path
from .core import parse_best_energy
from .contacts import detect_contacts, write_contacts


def parse_dlg_clusters(dlg: Path) -> list[dict[str, object]]:
    rows = []
    text = dlg.read_text(encoding="utf-8", errors="replace").splitlines()
    in_hist = False
    for line in text:
        if "CLUSTERING HISTOGRAM" in line:
            in_hist = True
            continue
        if in_hist and "Number of multi-member" in line:
            break
        m = re.match(r"^\s*(\d+)\s*\|\s*([-+]?\d+\.\d+)\s*\|\s*(\d+)\s*\|\s*([-+]?\d+\.\d+)\s*\|\s*(\d+)\s*\|", line)
        if m:
            rank, low, run, mean, size = int(m.group(1)), float(m.group(2)), int(m.group(3)), float(m.group(4)), int(m.group(5))
            rows.append({"ligand": dlg.stem, "cluster_rank": rank, "cluster_size": size, "cluster_percentage": 100.0 * size / _dlg_run_count(text), "lowest_binding_energy": low, "mean_binding_energy": mean, "run_number": run, "pose_model": run})
    return rows


def _dlg_run_count(lines: list[str]) -> int:
    for line in lines:
        m = re.search(r"Number of conformations\s*=\s*(\d+)", line)
        if m: return int(m.group(1))
    return sum(1 for line in lines if line.startswith("DOCKED: MODEL"))


def _extract_models(dlg: Path) -> list[tuple[int, float | None, list[str]]]:
    models: list[tuple[int, float | None, list[str]]] = []
    current_number: int | None = None
    current_energy: float | None = None
    current_atoms: list[str] = []
    for line in dlg.read_text(encoding="utf-8", errors="replace").splitlines():
        model_match = re.match(r"DOCKED:\s+MODEL\s+(\d+)", line)
        if model_match:
            current_number = int(model_match.group(1))
            current_energy = None
            current_atoms = []
            continue
        if current_number is None:
            continue
        energy_match = re.search(r"Estimated Free Energy of Binding\s*=\s*([+-]?\d+(?:\.\d+)?)", line)
        if energy_match:
            current_energy = float(energy_match.group(1))
        if line.startswith("DOCKED: ATOM") or line.startswith("DOCKED: HETATM"):
            current_atoms.append(line[8:])
        if line.startswith("DOCKED: ENDMDL"):
            models.append((current_number, current_energy, current_atoms))
            current_number = None
    return models


def _write_pose_and_complexes(dlg: Path, poses_dir: Path, complexes_dir: Path, receptor: Path) -> list[dict[str, object]]:
    poses_dir.mkdir(parents=True, exist_ok=True)
    complexes_dir.mkdir(parents=True, exist_ok=True)
    receptor_text = receptor.read_text(encoding="utf-8", errors="replace").rstrip() + "\n"
    outputs = []
    for model_number, energy, atoms in _extract_models(dlg):
        if not atoms:
            continue
        stem = f"{dlg.stem}_model_{model_number}"
        pose = poses_dir / f"{stem}.pdbqt"
        pose.write_text(f"MODEL {model_number}\n" + "\n".join(atoms) + "\nENDMDL\n", encoding="utf-8")
        complex_path = complexes_dir / f"{stem}_complex.pdbqt"
        complex_path.write_text(receptor_text + pose.read_text(encoding="utf-8"), encoding="utf-8")
        outputs.append({"model": model_number, "energy": energy, "pose": pose, "complex": complex_path})
    return outputs


def analyze_docking_results(dlg_dir: Path, output_dir: Path, top_n: int = 20) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    poses_dir = output_dir / "poses"
    complexes_dir = output_dir / "complexes"
    receptor = dlg_dir.parent.parent / "receptor" / "receptor.pdbqt"
    rows=[]
    cluster_rows=[]
    for dlg in sorted(dlg_dir.glob('*.dlg')):
        e=parse_best_energy(dlg)
        if e is not None:
            generated = _write_pose_and_complexes(dlg, poses_dir, complexes_dir, receptor)
            cluster_rows.extend(parse_dlg_clusters(dlg))
            best = min((item for item in generated if item["energy"] is not None), key=lambda item: float(item["energy"]), default=None)
            row = {'ligand':dlg.stem,'energy_kcal_mol':f'{e:.3f}','dlg':str(dlg)}
            if best is not None:
                row.update({'pose_file': str(best['pose']), 'complex_file': str(best['complex'])})
            rows.append(row)
    rows.sort(key=lambda r: float(r['energy_kcal_mol']))
    for i,r in enumerate(rows,1): r['rank']=i
    ranking=output_dir/'ranking.csv'
    with ranking.open('w',newline='',encoding='utf-8') as fh:
        w=csv.DictWriter(fh,fieldnames=['rank','ligand','energy_kcal_mol','dlg','pose_file','complex_file'], extrasaction='ignore'); w.writeheader(); w.writerows(rows)
    contacts=output_dir/'contacts.csv'
    all_contacts=[]
    for row in rows:
        pose = Path(row.get('pose_file', ''))
        if pose.exists(): all_contacts.extend(detect_contacts(receptor, pose, row['ligand']))
    write_contacts(all_contacts, contacts)
    clusters=output_dir/'clusters.csv'
    with clusters.open('w', newline='', encoding='utf-8') as fh:
        fields=['ligand','cluster_rank','cluster_size','cluster_percentage','lowest_binding_energy','mean_binding_energy','run_number','pose_model']
        w=csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(cluster_rows)
    comparison=output_dir/'comparison.csv'
    with comparison.open('w',newline='',encoding='utf-8') as fh:
        w=csv.DictWriter(fh,fieldnames=['rank','ligand','energy_kcal_mol']); w.writeheader(); w.writerows([{k:r[k] for k in ('rank','ligand','energy_kcal_mol')} for r in rows[:top_n]])
    return {'ranking':ranking,'contacts':contacts,'clusters':clusters,'comparison':comparison,'poses':poses_dir,'complexes':complexes_dir,'successful':len(rows)}
