from __future__ import annotations

import csv
from pathlib import Path
from .core import parse_best_energy


def analyze_docking_results(dlg_dir: Path, output_dir: Path, top_n: int = 20) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows=[]
    for dlg in sorted(dlg_dir.glob('*.dlg')):
        e=parse_best_energy(dlg)
        if e is not None:
            rows.append({'ligand':dlg.stem,'energy_kcal_mol':f'{e:.3f}','dlg':str(dlg)})
    rows.sort(key=lambda r: float(r['energy_kcal_mol']))
    for i,r in enumerate(rows,1): r['rank']=i
    ranking=output_dir/'ranking.csv'
    with ranking.open('w',newline='',encoding='utf-8') as fh:
        w=csv.DictWriter(fh,fieldnames=['rank','ligand','energy_kcal_mol','dlg']); w.writeheader(); w.writerows(rows)
    contacts=output_dir/'contacts.csv'
    contacts.write_text('ligand,residue,interaction_type,distance_A\n',encoding='utf-8')
    clusters=output_dir/'clusters.csv'
    clusters.write_text('ligand,cluster_rank,members,best_energy_kcal_mol\n',encoding='utf-8')
    comparison=output_dir/'comparison.csv'
    with comparison.open('w',newline='',encoding='utf-8') as fh:
        w=csv.DictWriter(fh,fieldnames=['rank','ligand','energy_kcal_mol']); w.writeheader(); w.writerows([{k:r[k] for k in ('rank','ligand','energy_kcal_mol')} for r in rows[:top_n]])
    return {'ranking':ranking,'contacts':contacts,'clusters':clusters,'comparison':comparison,'successful':len(rows)}
