from __future__ import annotations

from pathlib import Path


def export_ligand_complex(receptor_pdb: Path, ligand_template_pdb: Path, ligand_coordinates_pdb: Path, destination: Path, *, resname: str = "1N1", chain: str = "A", resseq: str = "501") -> Path:
    """Export a PLIP-compatible complex using HETATM identity from the crystal ligand."""
    template = [x for x in ligand_template_pdb.read_text().splitlines() if x.startswith(("ATOM", "HETATM")) and x[17:20].strip() == resname and x[21] == chain and x[22:27].strip() == resseq]
    coords = [x for x in ligand_coordinates_pdb.read_text().splitlines() if x.startswith(("ATOM", "HETATM"))]
    if len(template) != len(coords):
        raise ValueError("Ligand template and coordinate pose have different atom counts")
    ligand = [t[:30] + c[30:54] + t[54:] for t, c in zip(template, coords)]
    receptor = [x for x in receptor_pdb.read_text().splitlines() if not (x.startswith(("ATOM", "HETATM", "CONECT")) and x[17:20].strip() == resname and x[21] == chain and x[22:27].strip() == resseq)]
    conect = [x for x in ligand_template_pdb.read_text().splitlines() if x.startswith("CONECT")]
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(receptor + ligand + conect + ["END"]) + "\n")
    return destination
