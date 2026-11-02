from pathlib import Path

from sidock_ubuntu import stages


def test_validation_stage_forwards_authoritative_ccd_to_redocking(monkeypatch, tmp_path):
    seen = {}

    def fake_prepare(native_pdb, output_dir, ccd_path=None):
        seen["stage_prepare"] = ccd_path
        return {"crystal": object()}

    def fake_redocking(**kwargs):
        seen["redocking_ccd"] = kwargs["ccd_path"]
        return {"status": "TEST"}

    ccd = tmp_path / "input" / "ccd.json"
    ccd.parent.mkdir()
    ccd.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(stages, "prepare_native_ligand", fake_prepare)
    monkeypatch.setattr(stages, "run_native_redocking", fake_redocking)

    result = stages.run_validation_stage(
        tmp_path / "native.pdb",
        tmp_path / "receptor.pdbqt",
        tmp_path / "grid",
        tmp_path / "work",
        tmp_path / "output",
        object(),
        object(),
        {},
        ccd_path=ccd,
    )

    assert result == {"status": "TEST"}
    assert seen == {"stage_prepare": ccd, "redocking_ccd": ccd}
