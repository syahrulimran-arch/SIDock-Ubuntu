from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from .core import parse_best_energy
from .contacts import detect_contacts, write_contacts

@dataclass
class PoseExtractionResult:
    dlg_path: Path
    poses: list[dict]
    selected_pose_path: Path | None
    selected_model: int | None
    lowest_energy: float | None
    energy_ranking: list[dict]
    mapping_status: str
    warnings: list[str]
    provenance: dict

@dataclass
class ClusterAnalysisResult:
    cluster_count: int
    clusters: list[dict]
    largest_cluster_rank: int | None
    largest_cluster_population: int | None
    lowest_energy_cluster_rank: int | None
    cluster_table_path: Path | None
    warnings: list[str]

def parse_dlg_models(dlg: Path) -> list[dict]:
    """Read AutoDock models without contacts, figures, or summary side effects."""
    return [{'model': n, 'energy': e, 'atoms': atoms} for n,e,atoms in _extract_models(dlg)]

def rank_poses_by_energy(poses: list[dict]) -> list[dict]:
    return sorted([p for p in poses if p.get('energy') is not None], key=lambda p: (float(p['energy']), int(p['model'])))

def extract_docking_poses(dlg: Path, poses_dir: Path, complexes_dir: Path, receptor: Path) -> PoseExtractionResult:
    models=parse_dlg_models(dlg); poses_dir.mkdir(parents=True,exist_ok=True); complexes_dir.mkdir(parents=True,exist_ok=True)
    receptor_text=receptor.read_text(encoding='utf-8',errors='replace').rstrip()+'\n'; outputs=[]
    for item in models:
        if not item['atoms']: continue
        model=item['model']; stem=f'{dlg.stem}_model_{model}'; pose=poses_dir/f'{stem}.pdbqt'; complex_path=complexes_dir/f'{stem}_complex.pdbqt'
        pose.write_text(f'MODEL {model}\n'+'\n'.join(item['atoms'])+'\nENDMDL\n',encoding='utf-8')
        complex_path.write_text(receptor_text+pose.read_text(encoding='utf-8'),encoding='utf-8')
        outputs.append({'model':model,'energy':item['energy'],'pose':pose,'complex':complex_path,'mapping_status':'traversal_order_preserved_for_corrected_mapper'})
    ranked=rank_poses_by_energy(outputs); best=ranked[0] if ranked else None
    return PoseExtractionResult(dlg,outputs,best['pose'] if best else None,best['model'] if best else None,float(best['energy']) if best else None,ranked,'PASS' if outputs else 'FAIL',[],{'model_count':len(outputs)})

def analyze_autodock_clusters(dlg: Path, destination: Path | None = None) -> ClusterAnalysisResult:
    rows=parse_dlg_clusters(dlg)
    largest=max(rows,key=lambda r:int(r['cluster_size']),default=None)
    lowest=min(rows,key=lambda r:float(r['lowest_binding_energy']),default=None)
    if destination:
        destination.parent.mkdir(parents=True,exist_ok=True)
        fields=['ligand','cluster_rank','cluster_size','cluster_percentage','lowest_binding_energy','mean_binding_energy','run_number','pose_model']
        with destination.open('w',newline='',encoding='utf-8') as fh:
            w=csv.DictWriter(fh,fieldnames=fields); w.writeheader(); w.writerows(rows)
    return ClusterAnalysisResult(len(rows),rows,int(largest['cluster_rank']) if largest else None,int(largest['cluster_size']) if largest else None,int(lowest['cluster_rank']) if lowest else None,destination,[])


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
    if rows:
        return rows
    # AutoDock4 emits a model-based summary instead of the pipe-delimited
    # histogram when only one conformation is available (common in smoke
    # tests and deliberately small technical runs). Preserve that valid
    # result rather than silently producing an empty cluster table.
    total = _dlg_run_count(text)
    current = None
    for line in text:
        model = re.match(r"^\s*MODEL\s+(\d+)", line)
        if model:
            if current and {"cluster_rank", "cluster_size", "lowest_binding_energy"}.issubset(current):
                current["cluster_percentage"] = 100.0 * int(current["cluster_size"]) / max(total, 1)
                current["mean_binding_energy"] = current["lowest_binding_energy"]
                current["pose_model"] = current["run_number"]
                rows.append(current)
            current = {"ligand": dlg.stem, "pose_model": int(model.group(1))}
            continue
        if current is None:
            continue
        run = re.search(r"USER\s+Run\s*=\s*(\d+)", line)
        rank = re.search(r"USER\s+Cluster Rank\s*=\s*(\d+)", line)
        size = re.search(r"USER\s+Number of conformations in this cluster\s*=\s*(\d+)", line)
        energy = re.search(r"USER\s+Estimated Free Energy of Binding\s*=\s*([+-]?\d+(?:\.\d+)?)", line)
        if run:
            current["run_number"] = int(run.group(1))
        if rank:
            current["cluster_rank"] = int(rank.group(1))
        if size:
            current["cluster_size"] = int(size.group(1))
        if energy and "lowest_binding_energy" not in current:
            current["lowest_binding_energy"] = float(energy.group(1))
    if current and {"cluster_rank", "cluster_size", "lowest_binding_energy"}.issubset(current):
        current["cluster_percentage"] = 100.0 * int(current["cluster_size"]) / max(total, 1)
        current["mean_binding_energy"] = current["lowest_binding_energy"]
        current["run_number"] = current.get("run_number", current.get("pose_model", 1))
        current["pose_model"] = current["run_number"]
        rows.append(current)
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


def analyze_docking_results(dlg_dir: Path, output_dir: Path, top_n: int = 20) -> dict[str, object]:
    """Deprecated compatibility orchestrator over canonical pose/cluster APIs."""
    output_dir.mkdir(parents=True, exist_ok=True)
    poses_dir = output_dir / "poses"
    complexes_dir = output_dir / "complexes"
    receptor = dlg_dir.parent.parent / "receptor" / "receptor.pdbqt"
    rows=[]
    cluster_rows=[]
    for dlg in sorted(dlg_dir.glob('*.dlg')):
        e=parse_best_energy(dlg)
        if e is not None:
            pose_result = extract_docking_poses(dlg, poses_dir, complexes_dir, receptor)
            cluster_result = analyze_autodock_clusters(dlg)
            generated = pose_result.poses
            cluster_rows.extend(cluster_result.clusters)
            best = next((x for x in generated if x.get('model') == pose_result.selected_model), None)
            row = {'ligand':dlg.stem,'energy_kcal_mol':f'{e:.3f}','dlg':str(dlg)}
            if best is not None:
                row.update({'pose_file': str(best['pose']), 'complex_file': str(best['complex'])})
            rows.append(row)
    rows.sort(key=lambda r: float(r['energy_kcal_mol']))
    for i,r in enumerate(rows,1): r['rank']=i
    ranking=output_dir/'ranking.csv'
    with ranking.open('w',newline='',encoding='utf-8') as fh:
        w=csv.DictWriter(fh,fieldnames=['rank','ligand','energy_kcal_mol','dlg','pose_file','complex_file'], extrasaction='ignore'); w.writeheader(); w.writerows(rows)
    # Explicitly secondary legacy diagnostic; never primary PLIP evidence.
    contacts=output_dir/'interaction_crosscheck'/'interactions_sidock_experimental.csv'
    contacts.parent.mkdir(parents=True, exist_ok=True)
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
