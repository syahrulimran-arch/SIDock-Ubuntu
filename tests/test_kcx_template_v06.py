"""Offline KCX template regression tests for Meeko 0.8 compatibility."""

import json
from pathlib import Path


from sidock_ubuntu import core


TEMPLATE = Path(core.__file__).with_name("data") / "kcx_residue_template.json"


def test_kcx_template_has_ccd_side_chain_and_peptide_links():
    data = json.loads(TEMPLATE.read_text(encoding="utf-8"))
    record = data["residue_templates"]["KCX"]
    assert {"CX", "OQ1", "OQ2"}.issubset(record["atom_name"])
    assert record["link_labels"] == {"0": "C-term", "4": "N-term"}
