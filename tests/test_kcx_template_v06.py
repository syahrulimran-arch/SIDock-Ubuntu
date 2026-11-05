"""Offline KCX template regression tests for Meeko 0.8 compatibility."""

import json
from pathlib import Path

from meeko import ResidueChemTemplates

from sidock_ubuntu import core


TEMPLATE = Path(core.__file__).with_name("data") / "kcx_residue_template.json"


def test_kcx_template_has_ccd_side_chain_and_peptide_links():
    data = json.loads(TEMPLATE.read_text(encoding="utf-8"))
    record = data["residue_templates"]["KCX"]
    assert {"CX", "OQ1", "OQ2"}.issubset(record["atom_name"])
    assert record["link_labels"] == {"0": "C-term", "4": "N-term"}

    templates = ResidueChemTemplates.create_from_defaults()
    templates.add_json_file(str(TEMPLATE))
    kcx = templates.residue_templates["KCX"]
    assert kcx.atom_names[20:23] == ["CX", "OQ1", "OQ2"]
    assert kcx.link_labels == {0: "C-term", 4: "N-term"}
    assert kcx.mol.GetNumAtoms() == len(kcx.atom_names)


def test_prepare_receptor_adds_offline_kcx_template(tmp_path, monkeypatch):
    source = tmp_path / "protein.pdb"
    source.write_text(
        "HETATM    1  C   ACE A   0       0.000   0.000   0.000  1.00 20.00           C  \n"
        "HETATM    1  N   KCX A 220       0.000   0.000   0.000  1.00 20.00           N  \n",
        encoding="utf-8",
    )
    captured = {}

    monkeypatch.setattr(core, "resolve_meeko_receptor_preparer", lambda: {
        "available": True,
        "path": "/fake/mk_prepare_receptor.py",
    })

    def fake_run(cmd, cwd, log_path):
        captured["cmd"] = cmd
        log_path.write_text("fake", encoding="utf-8")
        return 1

    monkeypatch.setattr(core, "run_command", fake_run)
    try:
        core.prepare_receptor(source, tmp_path / "out.pdbqt", tmp_path / "meeko.log")
    except RuntimeError:
        pass
    assert "--add_templates" in captured["cmd"]
    assert str(TEMPLATE) in captured["cmd"]
    assert captured["cmd"][captured["cmd"].index("--blunt_ends") + 1] == "A:0=0"
