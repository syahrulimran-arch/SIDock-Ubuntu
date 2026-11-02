import json
from pathlib import Path

from sidock_ubuntu.config import DEFAULT_CONFIG, deep_merge, validate_config, summary
from sidock_ubuntu.dashboard import write_dashboard, write_results
from sidock_ubuntu.errors import ERROR_CODES, SIDockError, error_payload
from sidock_ubuntu.state import STAGES, project_status, update_state


def test_default_config_has_valid_gate():
    assert validate_config(DEFAULT_CONFIG) == []
    assert DEFAULT_CONFIG["validation"]["rmsd_threshold"] == 2.0

def test_config_merge_preserves_defaults():
    c = deep_merge(DEFAULT_CONFIG, {"validation": {"ga_runs": 20}})
    assert c["validation"]["ga_runs"] == 20
    assert c["docking"]["population_size"] == 150

def test_config_rejects_bad_threshold():
    c = deep_merge(DEFAULT_CONFIG, {"validation": {"rmsd_threshold": 0}})
    assert validate_config(c)

def test_config_rejects_bad_engine():
    c = deep_merge(DEFAULT_CONFIG, {"interaction": {"primary_engine": "foo"}})
    assert validate_config(c)

def test_config_summary_is_concise():
    s = summary(DEFAULT_CONFIG)
    assert "Population: 150" in s and "PLIP" in s

def test_error_codes_are_stable():
    assert ERROR_CODES["VALIDATION_RMSD_EXCEEDED"] == "E302"

def test_error_payload_is_machine_readable():
    payload = error_payload(SIDockError("INVALID_SDF", "bad input"))
    assert payload["code"] == "E002" and payload["error"] == "INVALID_SDF"

def test_state_has_all_stages(tmp_path):
    state = project_status(tmp_path)
    assert set(STAGES) == set(state["stages"])

def test_state_checkpoint_records_status(tmp_path):
    update_state(tmp_path, "INPUT", "PASS", warnings=["example"])
    state = project_status(tmp_path)
    assert state["stages"]["INPUT"]["status"] == "PASS"
    assert state["stages"]["INPUT"]["warnings"] == ["example"]

def test_state_checkpoint_records_hashes(tmp_path):
    update_state(tmp_path, "REFERENCE", "PASS", inputs={"pdb": "abc"})
    assert project_status(tmp_path)["stages"]["REFERENCE"]["input_hashes"]["pdb"] == "abc"

def test_results_schema(tmp_path):
    (tmp_path / "analysis").mkdir()
    (tmp_path / "analysis/evidence.json").write_text(json.dumps({"compounds": [{"compound_id": "x", "lowest_docking_energy": "-1"}]}))
    out = write_results(tmp_path)
    assert "compound_id" in out.read_text() and "x" in out.read_text()

def test_dashboard_contains_target(tmp_path):
    (tmp_path / "analysis").mkdir()
    (tmp_path / "analysis/evidence.json").write_text(json.dumps({"target": {"pdb": "TEST"}, "validation": {}, "compounds": []}))
    out = write_dashboard(tmp_path)
    assert "TEST" in out.read_text()

def test_dashboard_escapes_html(tmp_path):
    (tmp_path / "analysis").mkdir()
    (tmp_path / "analysis/evidence.json").write_text(json.dumps({"target": {"pdb": "<x>"}, "validation": {}, "compounds": []}))
    assert "&lt;x&gt;" in write_dashboard(tmp_path).read_text()

def test_resume_completed_state_is_detectable(tmp_path):
    for stage in STAGES:
        update_state(tmp_path, stage, "PASS")
    assert all(v["status"] == "PASS" for v in project_status(tmp_path)["stages"].values())

def test_failure_state_is_visible(tmp_path):
    update_state(tmp_path, "VALIDATION", "FAIL", errors=["E302"])
    assert project_status(tmp_path)["stages"]["VALIDATION"]["errors"] == ["E302"]

def test_results_preserve_warnings(tmp_path):
    (tmp_path / "analysis").mkdir()
    (tmp_path / "analysis/evidence.json").write_text(json.dumps({"compounds": [{"compound_id": "x", "warnings": "review"}]}))
    assert "review" in write_results(tmp_path).read_text()
