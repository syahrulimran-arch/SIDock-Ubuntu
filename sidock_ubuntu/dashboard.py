from __future__ import annotations
import csv, html, json
from pathlib import Path

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
                         "pi_stack_count": c.get("PLIP_pi_stack", 0), "salt_bridge_count": c.get("PLIP_salt_bridge", 0), "warnings": c.get("warnings", "")})
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
        rows.append(f"<tr><td>{html.escape(str(c.get('compound_id','')))}</td><td>{html.escape(str(c.get('status',c.get('final_status',c.get('preparation_status','COMPLETE')))))}</td><td>{html.escape(str(c.get('lowest_docking_energy','')))}</td><td>{html.escape(str(c.get('cluster_count','')))}</td><td>{html.escape(str(c.get('PLIP_H_bonds','')))} H / {html.escape(str(c.get('PLIP_hydrophobic','')))} hydrophobic</td></tr>")
    out = project / "report.html"
    content=("<!doctype html><meta charset='utf-8'><title>SIDock project dashboard</title><style>body{font:15px system-ui;max-width:1000px;margin:40px auto;color:#243047}table{border-collapse:collapse;width:100%}td,th{border:1px solid #ccd;padding:8px;text-align:left}code{background:#eef;padding:2px 5px}</style>" +
             f"<h1>SIDock Ubuntu project dashboard</h1><p><b>Target:</b> {html.escape(str(target.get('pdb','not recorded')))}<br><b>Reference:</b> {html.escape(str(target.get('reference_ligand','not recorded')))}<br><b>Validation:</b> {html.escape(str(validation.get('sampling_result','not recorded')))}; RMSD {html.escape(str(validation.get('native_redocking_RMSD_A','not recorded')))} Å</p><table><tr><th>Compound</th><th>Status</th><th>Docking energy</th><th>Clusters</th><th>PLIP interactions</th></tr>{''.join(rows)}</table>")
    if validation.get('sampling_result') == 'FAIL':
        rmsd = html.escape(str(validation.get('native_redocking_RMSD_A', 'unavailable')))
        threshold = html.escape(str(validation.get('rmsd_threshold', 'unavailable')))
        content += f"<p><b>Native redocking validation failed:</b> minimum RMSD {rmsd} Å exceeded the acceptance threshold {threshold} Å; compound docking was not performed.</p>"
    tmp=out.with_suffix('.html.tmp'); tmp.write_text(content,encoding='utf-8'); tmp.replace(out)
    return out
