import csv
import hashlib
import json

from sidock_ubuntu.dashboard import write_results, write_dashboard
from sidock_ubuntu.scientific_report import build_summary, write_scientific_report
from sidock_ubuntu.stages import run_finalization_stage


def _fixture(tmp_path):
    p = tmp_path / "study"
    (p / "analysis").mkdir(parents=True)
    (p / "docking/dpf").mkdir(parents=True)
    (p / "docking/dlg").mkdir(parents=True)
    evidence = {
        "project_status": "COMPLETE",
        "target": {"pdb": "protein.pdb", "reference_ligand": "GLC:A:601"},
        "validation": {"sampling_result": "PASS", "native_redocking_RMSD_A": 0.340, "rmsd_threshold": 2.0},
        "compounds": [
            {"compound_id": "7i", "status": "COMPLETE", "lowest_docking_energy": -10.38,
             "cluster_count": 46, "largest_cluster_population": 12,
             "PLIP_H_bonds": 1, "PLIP_hydrophobic": 1},
            {"compound_id": "7j", "status": "COMPLETE", "lowest_docking_energy": -9.60,
             "cluster_count": 56, "largest_cluster_population": 9,
             "PLIP_H_bonds": 0, "PLIP_hydrophobic": 0},
        ],
    }
    (p / "analysis/evidence.json").write_text(json.dumps(evidence))
    (p / "docking/dpf/7i.dpf").write_text("ga_run 100\nga_pop_size 150\n")
    (p / "docking/dlg/7i.dlg").write_bytes(b"immutable docking log")
    with (p / "analysis/interactions_plip_7i.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "interaction_type", "protein_chain", "protein_residue_name", "protein_residue_number"])
        w.writeheader()
        w.writerow({"interaction_type": "hydrogen_bond", "protein_chain": "A",
                    "protein_residue_name": "ASP", "protein_residue_number": "215"})
        w.writerow({"interaction_type": "hydrophobic", "protein_chain": "A",
                    "protein_residue_name": "TRP", "protein_residue_number": "150"})
    with (p / "analysis/clusters_7i.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["ligand", "cluster_size"])
        w.writeheader()
        w.writerow({"ligand": "7i", "cluster_size": "12"})
        w.writerow({"ligand": "7i", "cluster_size": "8"})
    return p, evidence


def test_draft_uses_existing_evidence_and_does_not_modify_dlg(tmp_path):
    p, _ = _fixture(tmp_path)
    before = hashlib.sha256((p / "docking/dlg/7i.dlg").read_bytes()).hexdigest()
    output = write_scientific_report(p)
    assert output["docking_rerun"] is False
    assert output["ranked_compounds"] == 2
    assert hashlib.sha256((p / "docking/dlg/7i.dlg").read_bytes()).hexdigest() == before
    text = (p / "analysis/scientific_report.md").read_text()
    assert "7i showed the most favorable AutoDock4 score (-10.38 kcal/mol)" in text
    assert "ASP215 chain A" in text and "TRP150 chain A" in text
    assert "GA runs per compound: 100" in text
    assert "No numeric experimental IC50" in text
    assert (p / "analysis/scientific_report.html").exists()


def test_missing_plip_stays_unavailable_not_fabricated(tmp_path):
    p, _ = _fixture(tmp_path)
    data = build_summary(p)
    a = {x["compound_id"]: x for x in data["compounds"]}
    assert a["7j"]["plip_counts"]["hydrogen_bond"] is None
    assert a["7j"]["residues_by_interaction"] == {}
    assert a["7i"]["largest_cluster_percentage"] == 60.0
    assert data["ranked_ids"] == ["7i", "7j"]


def test_finalization_includes_scientific_report_and_dashboard_link(tmp_path):
    p, evidence = _fixture(tmp_path)
    run_finalization_stage(p, evidence, write_results_fn=write_results, write_dashboard_fn=write_dashboard)
    assert (p / "analysis/scientific_report.json").is_file()
    assert "analysis/scientific_report.html" in (p / "report.html").read_text()


def test_assay_uses_only_explicit_recorded_values(tmp_path):
    p, _ = _fixture(tmp_path)
    with (p / "analysis/assay.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["ligand", "ic50_uM"])
        w.writeheader()
        w.writerow({"ligand": "7i", "ic50_uM": "3.5"})
    result = build_summary(p)
    assert result["activity_correlation"]["n"] == 1
    assert result["activity_correlation"]["pearson_r"] is None
    assert next(x for x in result["compounds"] if x["compound_id"] == "7i")["ic50_uM"] == 3.5

def test_distance_rebuilt_from_saved_plip_xml_without_docking(tmp_path):
    p, _ = _fixture(tmp_path)
    raw = p / "analysis/plip_7i/plip"
    raw.mkdir(parents=True)
    (raw / "binding_report.xml").write_text(
        "<report><bindingsite><hydrogen_bonds><hydrogen_bond>"
        "<reschain>A</reschain><restype>ASP</restype><resnr>215</resnr>"
        "<dist_d-a>2.81</dist_d-a><don_angle>155</don_angle>"
        "</hydrogen_bond></hydrogen_bonds>"
        "<hydrophobic_interactions><hydrophobic_interaction>"
        "<reschain>A</reschain><restype>TRP</restype><resnr>150</resnr>"
        "<dist>3.67</dist>"
        "</hydrophobic_interaction></hydrophobic_interactions>"
        "</bindingsite></report>")
    original = (p / "analysis/interactions_plip_7i.csv").read_bytes()
    result = write_scientific_report(p)
    assert result["docking_rerun"] is False
    assert (p / "analysis/interactions_plip_7i.csv").read_bytes() == original
    report = (p / "analysis/scientific_report.md").read_text()
    assert "2.81" in report and "3.67" in report
    assert "155.0" in report
    details = (p / "analysis/interactions_detailed_7i.csv").read_text()
    assert "2.81" in details and "3.67" in details
