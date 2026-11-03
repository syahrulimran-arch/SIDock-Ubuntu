from __future__ import annotations
import csv, html, json
from pathlib import Path

def _display(value, fallback='not recorded'):
    return fallback if value in (None,'') else value

def write_results(project: Path) -> Path:
    evidence = project / "analysis" / "evidence.json"
    out = project / "results.csv"
    rows = []
    if evidence.exists():
        data = json.loads(evidence.read_text())
        for c in data.get("compounds", []):
            rows.append({"compound_id": c.get("compound_id", c.get("ligand", "")), "status": c.get("status", c.get("final_status", c.get("preparation_status", "COMPLETE"))),
                         "lowest_energy": c.get("lowest_docking_energy", ""), "cluster_count": c.get("cluster_count", ""),
                         "largest_cluster_population": c.get("largest_cluster_population", ""), "selected_pose": c.get("selected_pose", ""),
                         "hbond_count": c.get("PLIP_H_bonds", ""), "hydrophobic_count": c.get("PLIP_hydrophobic", ""),
                         "pi_stack_count": c.get("PLIP_pi_stack", ""), "salt_bridge_count": c.get("PLIP_salt_bridge", ""), "warnings": c.get("warnings", "")})
    fields = ["compound_id","status","lowest_energy","cluster_count","largest_cluster_population","selected_pose","hbond_count","hydrophobic_count","pi_stack_count","salt_bridge_count","warnings"]
    tmp=out.with_suffix('.csv.tmp')
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)
    tmp.replace(out)
    return out

def write_dashboard(project: Path) -> Path:
    evidence = project / "analysis" / "evidence.json"
    data = json.loads(evidence.read_text()) if evidence.exists() else {}
    target = data.get("target", {})
    validation = data.get("validation", {})
    rows = []
    for c in data.get("compounds", []):
        cluster=_display(c.get('cluster_count'),'unavailable')
        hbond=c.get('PLIP_H_bonds'); hydrophobic=c.get('PLIP_hydrophobic')
        interactions='unavailable' if hbond is None or hydrophobic is None else f'{hbond} H / {hydrophobic} hydrophobic'
        rows.append(f"<tr><td>{html.escape(str(c.get('compound_id','')))}</td><td>{html.escape(str(c.get('status',c.get('final_status',c.get('preparation_status','COMPLETE')))))}</td><td>{html.escape(str(_display(c.get('lowest_docking_energy'),'unavailable')))}</td><td>{html.escape(str(cluster))}</td><td>{html.escape(str(interactions))}</td></tr>")
    out = project / "report.html"
    content=("<!doctype html><meta charset='utf-8'><title>SIDock project dashboard</title><style>body{font:15px system-ui;max-width:1000px;margin:40px auto;color:#243047}table{border-collapse:collapse;width:100%}td,th{border:1px solid #ccd;padding:8px;text-align:left}code{background:#eef;padding:2px 5px}</style>" +
             f"<h1>SIDock Ubuntu project dashboard</h1><p><b>Target:</b> {html.escape(str(_display(target.get('pdb'))))}<br><b>Reference:</b> {html.escape(str(_display(target.get('reference_ligand'))))}<br><b>Validation:</b> {html.escape(str(_display(validation.get('sampling_result'))))}; RMSD {html.escape(str(_display(validation.get('native_redocking_RMSD_A'),'unavailable')))} Å</p><table><tr><th>Compound</th><th>Status</th><th>Docking energy</th><th>Clusters</th><th>PLIP interactions</th></tr>{''.join(rows)}</table>")
    if validation.get('sampling_result') == 'FAIL':
        rmsd = html.escape(str(_display(validation.get('native_redocking_RMSD_A'), 'unavailable')))
        threshold = html.escape(str(_display(validation.get('rmsd_threshold'), 'unavailable')))
        content += f"<p><b>Native redocking validation failed:</b> minimum RMSD {rmsd} Å exceeded the acceptance threshold {threshold} Å; compound docking was not performed.</p>"
    tmp=out.with_suffix('.html.tmp'); tmp.write_text(content,encoding='utf-8'); tmp.replace(out)
    return out
