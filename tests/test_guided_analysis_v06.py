"""Synthetic DLG/PLIP fixtures verify orchestration, never scientific outcomes."""
import csv
import hashlib
import json
import shutil
import zipfile

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from sidock_ubuntu import guided_analysis as ga
from sidock_ubuntu.guided_workflow import create_guided_project
from sidock_ubuntu.manual_docking import root, execute
from sidock_ubuntu.plip_runner import normalize_plip_interactions


@pytest.fixture(autouse=True)
def optional_pymol_is_not_required_for_regressions(monkeypatch):
    original = shutil.which
    monkeypatch.setattr(shutil, "which", lambda name: None if name == "pymol" else original(name))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atom(serial, name, code, chain, resid, xyz, typ):
    x, y, z = xyz
    return (f"HETATM{serial:5d} {name:>4s} {code:>3s} {chain}{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00    {0.1:6.3f} {typ:>2s}\n")


def fixture_run(tmp_path, *, count=2):
    project = create_guided_project(tmp_path, "analysis-test")
    home = root(project)
    run = home / "runs" / ("a" * 16)
    for directory in ("docking/dlg", "docking/dpf", "grid/production"):
        (run / directory).mkdir(parents=True)
    receptor = run / "grid/production/receptor.pdbqt"
    receptor.write_text(atom(1, "N", "GLY", "A", 1, (0, 0, 0), "N") +
                        atom(2, "CA", "GLY", "A", 1, (1.4, 0, 0), "C") +
                        atom(3, "C", "GLY", "A", 1, (2, 1.3, 0), "C") +
                        atom(4, "O", "GLY", "A", 1, (1.7, 2.4, 0), "OA"))
    sdf = project / "input/compounds.sdf"
    ligand_hashes = {}
    with Chem.SDWriter(str(sdf)) as writer:
        for i in range(count):
            cid = "A-" + str(i + 1)
            mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
            AllChem.EmbedMolecule(mol, randomSeed=i + 123)
            mol.SetProp("Compound_ID", cid)
            mol.SetProp("_Name", cid)
            mol.SetProp("IC50", "250 nM" if i == 0 else ">10")
            writer.write(mol)
            ligand = (atom(1, "C1", "UNL", "X", 1, (1, 3.5, 0), "C") +
                      atom(2, "C2", "UNL", "X", 1, (2.4, 3.5, 0), "C") +
                      atom(3, "O1", "UNL", "X", 1, (3, 4.5, 0), "OA"))
            prepared = home / "ligands/prepared" / (cid + ".pdbqt")
            prepared.parent.mkdir(parents=True, exist_ok=True)
            prepared.write_text("ROOT\n" + ligand + "ENDROOT\nTORSDOF 0\n")
            staged = run / "docking/work" / cid / prepared.name
            staged.parent.mkdir(parents=True)
            shutil.copy2(prepared, staged)
            ligand_hashes[prepared.name] = sha(prepared)
            dlg = "Number of conformations = 2\n"
            for model, energy in ((1, -5.0), (2, -7.2 - i)):
                dlg += f"DOCKED: MODEL {model}\nDOCKED: USER Estimated Free Energy of Binding = {energy}\n"
                dlg += "".join("DOCKED: " + line + "\n" for line in ligand.splitlines())
                dlg += "DOCKED: ENDMDL\n"
            dlg += f"CLUSTERING HISTOGRAM\n  1 | {-7.2-i} | 2 | -6.0 | 2 |\nNumber of multi-member\nautodock4: Successful Completion\n"
            (run / "docking/dlg" / (cid + ".dlg")).write_text(dlg)
            (run / "docking/dpf" / (cid + ".dpf")).write_text("ga_run 2\nga_pop_size 150\n")
    (home / "ligands/manifest.json").write_text(json.dumps({
        "source_sha256": sha(sdf), "compounds": list(ligand_hashes)}))
    (run / "protocol.json").write_text(json.dumps({
        "run_id": run.name, "model_status": "EXPLORATORY_UNVALIDATED_NI_MODEL",
        "scientific_approval": False, "receptor_sha256": sha(receptor),
        "ligand_sha256": ligand_hashes, "target": {"pdb": "fixture.pdb", "reference_ligand": "LIG:A:100"}}))
    return project, run


def mock_plip(monkeypatch, *, zero=False):
    monkeypatch.setattr(ga, "resolve_plip_executable", lambda: {"available": True, "version": "fixture"})
    calls = []

    def run(complex_pdb, work, destination, **kwargs):
        calls.append(kwargs["compound_id"])
        report = work / "plip/report.xml"
        report.parent.mkdir(parents=True)
        contact = "" if zero else ("<hydrogen_bonds><hydrogen_bond><reschain>A</reschain>"
                                   "<restype>GLY</restype><resnr>1</resnr><dist_d-a>2.9</dist_d-a>"
                                   "</hydrogen_bond></hydrogen_bonds>")
        report.write_text("<report><bindingsite><identifiers><hetid>LIG</hetid><chain>Z</chain>"
                          "<position>999</position></identifiers>" + contact + "</bindingsite></report>")
        normalize_plip_interactions([report], destination, **kwargs)

    monkeypatch.setattr(ga, "run_plip_analysis", run)
    return calls


def test_saved_run_produces_reports_figures_archive_and_reuses_without_docking(tmp_path, monkeypatch):
    project, run = fixture_run(tmp_path)
    calls = mock_plip(monkeypatch)
    original = {p: sha(p) for p in run.rglob("*") if p.is_file()}
    result = ga.analyze(project, run_id=run.name)
    assert result["status"] == "COMPLETE_WITH_WARNINGS"
    assert result["analysed"] == result["plip_complete"] == 2
    out = project / result["analysis_directory"]
    evidence = json.loads((out / "evidence.json").read_text())
    assert evidence["scientific_approval"] is False
    assert evidence["validation"]["native_redocking_RMSD_A"] is None
    assert evidence["compounds"][0]["selected_model"] == 2
    report = (out / "scientific_report.md").read_text()
    assert "SCIENTIFIC APPROVAL: FALSE" in report
    assert "GLY1 chain A" in report and "2.90" in report
    assert "GA runs per compound: 2" in report
    assert "0.250" in report  # explicit nM -> uM, not the censored >10 value
    assert (out / "figures/docking_ranking.png").is_file()
    assert (out / "figures/A-1_pose.png").is_file()
    assert (out / "figures/A-1_structure.png").is_file()
    assert (out / "figures/A-1_interactions_2d.png").is_file()
    assert "A-1_interactions_2d.png" in report
    assert (out / "figures/A-1_binding_pose.py").is_file()
    assert "<img" in (out / "scientific_report.html").read_text()
    with zipfile.ZipFile(project / result["bundle"]) as archive:
        assert {"scientific_report.md", "evidence.json", "figures/docking_ranking.png",
                "figures/A-1_interactions_2d.png",
                "sources/dlg/A-1.dlg"}.issubset(archive.namelist())
    assert all(sha(p) == expected for p, expected in original.items())
    cached = ga.analyze(project, run_id=run.name)
    assert cached["reused"] is True
    assert calls == ["A-1", "A-2"]


def test_missing_plip_is_partial_and_counts_stay_unavailable(tmp_path, monkeypatch):
    project, run = fixture_run(tmp_path, count=1)
    monkeypatch.setattr(ga, "resolve_plip_executable", lambda: {
        "available": False, "message": "fixture PLIP unavailable", "version": None})
    result = ga.analyze(project, run_id=run.name)
    assert result["status"] == "PARTIAL"
    data = json.loads((project / result["analysis_directory"] / "evidence.json").read_text())
    assert data["compounds"][0]["PLIP_H_bonds"] is None
    assert data["compounds"][0]["lowest_docking_energy"] == -7.2
    assert not (project / result["analysis_directory"] / "interactions_plip_A-1.csv").exists()


def test_successful_zero_contact_plip_is_distinct_from_unavailable(tmp_path, monkeypatch):
    project, run = fixture_run(tmp_path, count=1)
    mock_plip(monkeypatch, zero=True)
    result = ga.analyze(project, run_id=run.name)
    data = json.loads((project / result["analysis_directory"] / "evidence.json").read_text())
    assert data["compounds"][0]["PLIP_H_bonds"] == 0
    assert data["compounds"][0]["plip_status"] == "COMPLETE"


def test_changed_dlg_invalidates_report_and_keeps_previous_revision(tmp_path, monkeypatch):
    project, run = fixture_run(tmp_path, count=1)
    mock_plip(monkeypatch)
    first = ga.analyze(project, run_id=run.name)
    dlg = run / "docking/dlg/A-1.dlg"
    dlg.write_text(dlg.read_text().replace("-7.2", "-9.1"))
    assert ga.analysis_status(project, run)["status"] == "STALE"
    second = ga.analyze(project, run_id=run.name)
    assert second["analysis_directory"] != first["analysis_directory"]
    assert (project / first["analysis_directory"] / "evidence.json").is_file()


def test_older_analysis_version_is_not_reused_without_new_figures(tmp_path, monkeypatch):
    project, run = fixture_run(tmp_path, count=1)
    mock_plip(monkeypatch)
    first = ga.analyze(project, run_id=run.name)
    pointer = run / "analysis_result.json"
    data = json.loads(pointer.read_text())
    data["schema"] = "guided-analysis-1"
    pointer.write_text(json.dumps(data))
    assert ga.analysis_status(project, run)["status"] == "STALE"
    second = ga.analyze(project, run_id=run.name)
    assert second["analysis_directory"] != first["analysis_directory"]
    assert (project / second["analysis_directory"] / "figures/A-1_interactions_2d.png").is_file()


def test_failed_figures_can_be_retried_without_changing_docking_inputs(tmp_path, monkeypatch):
    project, run = fixture_run(tmp_path, count=1)
    mock_plip(monkeypatch)
    monkeypatch.setattr(ga, "generate_figures", lambda *args, **kwargs: ([], ["fixture image unavailable"]))
    first = ga.analyze(project, run_id=run.name)
    assert first["retryable"] is True
    monkeypatch.setattr(ga, "generate_figures", lambda *args, **kwargs: ([], []))
    second = ga.analyze(project, run_id=run.name)
    assert second["analysis_directory"] != first["analysis_directory"]
    assert second["retryable"] is False


def test_incomplete_dlg_is_excluded_from_ranked_poses(tmp_path, monkeypatch):
    project, run = fixture_run(tmp_path)
    mock_plip(monkeypatch)
    dlg = run / "docking/dlg/A-2.dlg"
    dlg.write_text(dlg.read_text().replace("Successful Completion", "interrupted"))
    result = ga.analyze(project, run_id=run.name)
    assert result["status"] == "PARTIAL" and result["analysed"] == 1
    rows = json.loads((project / result["analysis_directory"] / "evidence.json").read_text())["compounds"]
    assert rows[1]["lowest_docking_energy"] is None
    assert rows[1]["status"] == "FAILED"


def test_receptor_and_pose_identity_mismatches_are_not_analysed(tmp_path, monkeypatch):
    project, run = fixture_run(tmp_path, count=1)
    mock_plip(monkeypatch)
    receptor = run / "grid/production/receptor.pdbqt"
    receptor.write_text(receptor.read_text() + "REMARK changed\n")
    with pytest.raises(RuntimeError, match="receptor differs"):
        ga.analyze(project, run_id=run.name)


def test_missing_original_sdf_does_not_force_redocking(tmp_path, monkeypatch):
    project, run = fixture_run(tmp_path, count=1)
    mock_plip(monkeypatch)
    (project / "input/compounds.sdf").unlink()
    result = ga.analyze(project, run_id=run.name)
    assert result["analysed"] == 1
    report = (project / result["analysis_directory"] / "scientific_report.md").read_text()
    assert "No numeric experimental IC50" in report


def test_wrong_pose_atom_identity_is_not_reported_as_a_valid_pose(tmp_path, monkeypatch):
    project, run = fixture_run(tmp_path, count=1)
    mock_plip(monkeypatch)
    dlg = run / "docking/dlg/A-1.dlg"
    dlg.write_text(dlg.read_text().replace("  C1", "  XX"))
    result = ga.analyze(project, run_id=run.name)
    assert result["status"] == "FAILED"
    evidence = json.loads((project / result["analysis_directory"] / "evidence.json").read_text())
    assert evidence["compounds"][0]["selected_pose"] is None
    assert evidence["compounds"][0]["lowest_docking_energy"] is None


def test_plip_filters_out_endogenous_ligand_and_rejects_missing_docked_site(tmp_path):
    report = tmp_path / "report.xml"
    contact = "<hydrogen_bonds><hydrogen_bond><resnr>123</resnr></hydrogen_bond></hydrogen_bonds>"
    report.write_text("<report><bindingsite><identifiers><hetid>OTHER</hetid><chain>A</chain>"
                      "<position>100</position></identifiers>" + contact + "</bindingsite>"
                      "<bindingsite><identifiers><hetid>LIG</hetid><chain>Z</chain>"
                      "<position>999</position></identifiers></bindingsite></report>")
    dest = normalize_plip_interactions([report], tmp_path / "contacts.csv", ligand_identity=("LIG", "Z", "999"))
    assert list(csv.DictReader(dest.open())) == []
    with pytest.raises(Exception, match="did not identify"):
        normalize_plip_interactions([report], tmp_path / "missing.csv", ligand_identity=("LIG", "Y", "999"))


def test_worker_automatically_runs_analysis_after_docking(tmp_path, monkeypatch):
    from sidock_ubuntu import manual_docking as manual
    project, run = fixture_run(tmp_path, count=1)
    calls = []
    monkeypatch.setattr(manual, "dock", lambda p: {"run_id": run.name, "completed": 1})
    monkeypatch.setattr(ga, "analyze", lambda p, **kw: calls.append(kw) or {"status": "PARTIAL"})
    result = execute(project, "dock")
    assert calls == [{"run_id": run.name}]
    assert result["completed"] == 1
    assert json.loads((root(project) / "job.json").read_text())["status"] == "SUCCEEDED_WITH_WARNINGS"


def test_docking_completion_is_preserved_when_analysis_crashes(tmp_path, monkeypatch):
    from sidock_ubuntu import manual_docking as manual
    project, run = fixture_run(tmp_path, count=1)
    monkeypatch.setattr(manual, "dock", lambda p: {"run_id": run.name, "completed": 1})
    def fail(*args, **kwargs):
        raise RuntimeError("fixture analysis failure")
    monkeypatch.setattr(ga, "analyze", fail)
    result = execute(project, "dock")
    assert result["completed"] == 1 and result["analysis"]["status"] == "FAILED"
    job = json.loads((root(project) / "job.json").read_text())
    assert "Docking completed; analysis failed" in job["error"]


def test_gui_runs_saved_analysis_without_docking_approval_and_links_downloads(tmp_path, monkeypatch):
    pytest.importorskip("httpx")
    from fastapi.testclient import TestClient
    from sidock_ubuntu.gui import create_app
    from sidock_ubuntu import manual_docking as manual
    project, run = fixture_run(tmp_path, count=1)
    mock_plip(monkeypatch)
    result = ga.analyze(project, run_id=run.name)
    calls = []
    monkeypatch.setattr(manual, "launch", lambda p, action, **kwargs: calls.append((action, kwargs)))
    with TestClient(create_app(tmp_path)) as client:
        response = client.post(f"/projects/{project.name}/guided/action", data={"action": "analyze", "run_id": run.name},
                               follow_redirects=False)
        assert response.status_code == 303
        html = client.get(response.headers["location"]).text
        assert "Download all analysis files (ZIP)" in html
        assert "2D protein–ligand interaction figures" in html
        assert "no redocking" in html
        download = client.get(f"/projects/{project.name}/artifacts/{result['bundle']}")
        assert download.status_code == 200
        assert download.content.startswith(b"PK")
        bad = client.post(f"/projects/{project.name}/guided/action", data={"action": "analyze", "run_id": "../escape"})
        assert bad.status_code == 409
    assert calls == [("analyze", {"run_id": run.name})]


@pytest.mark.parametrize("value, expected", [("250 nM", .25), ("0.5 mM", 500), (">10", None),
                                             ("2 +/- 1", None), ("3 ug/mL", None), ("-1", None)])
def test_activity_parser_does_not_invent_exact_measurements(value, expected):
    from sidock_ubuntu.assay import _float
    assert _float(value) == expected


def test_assay_respects_explicit_unit_metadata(tmp_path):
    from sidock_ubuntu.assay import extract_sdf_assay
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("Compound_ID", "A")
    mol.SetProp("IC50", "250")
    mol.SetProp("IC50_units", "nM")
    source = tmp_path / "assay.sdf"
    with Chem.SDWriter(str(source)) as writer:
        writer.write(mol)
    assert extract_sdf_assay(source, tmp_path / "assay.csv")[0]["ic50_uM"] == .25
