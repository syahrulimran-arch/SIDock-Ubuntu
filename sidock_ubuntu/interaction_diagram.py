"""Ligand-centred 2D diagrams from PLIP evidence, without inferred contacts."""
from __future__ import annotations

import io
import math
import re
import xml.etree.ElementTree as ET
from collections import OrderedDict
from pathlib import Path

from rdkit import Chem

STYLES = {
    "hydrogen_bond": ("Hydrogen bond", "#238b45", "--"),
    "hydrophobic": ("Hydrophobic", "#ac7133", ":"),
    "pi-stack": ("π stacking", "#9851a5", "--"),
    "pi-cation": ("π–cation", "#c467b7", "-."),
    "salt_bridge": ("Salt bridge", "#cc4c55", "-."),
    "halogen_bond": ("Halogen bond", "#159b9d", "--"),
    "water_bridge": ("Water bridge", "#3281b5", ":"),
    "metal_contact": ("Metal contact (PLIP)", "#4975a4", "--"),
}


def _text(node, path):
    value = node.findtext(path)
    return value.strip() if value else ""


def _ligand_graph(row, structure):
    """Use PLIP's canonical atom order only if it matches the supplied graph.

    SMILES-to-PDB is PLIP's own recorded atom correspondence, not an assumed
    match between RDKit atom order and PLIP atom indices.
    """
    supplied = Chem.RemoveHs(Chem.Mol(structure)) if structure is not None else None
    for name in row.get("plip_reports", []):
        for site in ET.parse(name).getroot().findall(".//bindingsite"):
            if tuple(_text(site, "identifiers/" + key) for key in ("hetid", "chain", "position")) != ("LIG", "Z", "999"):
                continue
            smiles = _text(site, "identifiers/smiles")
            parser = Chem.SmilesParserParams()
            parser.removeHs = False
            perceived = Chem.MolFromSmiles(smiles, parser) if smiles else None
            if perceived is None:
                continue
            # Keep explicit atom order from the recorded SMILES when drawing.
            # Never map its indices onto a different supplied chemical graph.
            if supplied is not None and Chem.MolToSmiles(supplied) != Chem.MolToSmiles(Chem.RemoveHs(perceived)):
                return supplied, {}, "Supplied structure; PLIP atom mapping does not match its chemical graph"
            mapping = {}
            raw_mapping = _text(site, "mappings/smiles_to_pdb")
            try:
                for pair in raw_mapping.split(","):
                    idx, serial = map(int, pair.split(":"))
                    if idx < 1 or idx > perceived.GetNumAtoms() or serial in mapping or idx - 1 in mapping.values():
                        raise ValueError("Invalid PLIP atom mapping")
                    mapping[serial] = idx - 1
                heavy = {a.GetIdx() for a in perceived.GetAtoms() if a.GetAtomicNum() > 1}
                if not heavy.issubset(mapping.values()):
                    raise ValueError("Incomplete PLIP atom mapping")
            except ValueError:
                mapping = {}
            note = "PLIP SMILES-to-PDB atom mapping"
            if supplied is None:
                note += "; structure perceived by PLIP, input structure verification unavailable"
            return perceived, mapping, note
    if supplied is None:
        raise RuntimeError("No supplied structure or usable PLIP ligand SMILES for the 2D interaction diagram")
    return supplied, {}, "Supplied structure; no usable PLIP atom mapping"


def _indices(raw):
    value = str(raw or "").strip()
    if not re.fullmatch(r"\d+(?:\s*[,;]\s*\d+)*", value):
        return ()
    return tuple(sorted(set(map(int, re.split(r"\s*[,;]\s*", value)))))


def draw_interaction_diagrams(row, structure, destination, *, per_page=12):
    """Draw every recorded contact, paginating dense diagrams; no docking calls."""
    if row.get("plip_status") != "COMPLETE":
        raise RuntimeError("PLIP did not complete; residue contacts cannot be drawn")
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from PIL import Image
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    mol, mapping, method = _ligand_graph(row, structure)
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DCairo(1000, 800)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.padding = .16
    options.minFontSize = 16
    drawer.DrawMolecule(mol)
    atom_points = {serial: drawer.GetDrawCoords(idx) for serial, idx in mapping.items()}
    drawer.FinishDrawing()
    ligand_image = Image.open(io.BytesIO(drawer.GetDrawingText()))
    grouped = OrderedDict()
    for contact in row.get("interactions", []):
        residue = (str(contact.get("protein_residue_name") or ""),
                   str(contact.get("protein_residue_number") or ""),
                   str(contact.get("protein_chain") or ""))
        indices = _indices(contact.get("ligand_atom"))
        kind = str(contact.get("interaction_type") or "unclassified")
        grouped.setdefault((*residue, kind, indices), []).append(contact)
    groups = list(grouped.items())
    def angle_order(entry):
        indices = entry[0][-1]
        if indices and all(idx in atom_points for idx in indices):
            x = sum(atom_points[idx].x for idx in indices) / len(indices) - 500
            y = sum(atom_points[idx].y for idx in indices) / len(indices) - 400
            return (math.atan2(y, x) + math.pi / 2) % (2 * math.pi)
        return 2 * math.pi
    groups.sort(key=angle_order)
    page_count = max(1, math.ceil(len(groups) / per_page))
    outputs, source_rows = [], []
    for page in range(page_count):
        selected = groups[page * per_page:(page + 1) * per_page]
        fig = plt.figure(figsize=(11, 9))
        ax = fig.add_axes((0, 0, 1, 1))
        ax.set(xlim=(0, 2200), ylim=(1800, 0))
        ax.set_axis_off()
        ax.imshow(ligand_image, extent=(600, 1600, 1300, 500), zorder=2)
        ax.text(1100, 65, row["compound_id"] + " — 2D protein–ligand interactions" +
                (f" ({page + 1}/{page_count})" if page_count > 1 else ""),
                ha="center", va="center", fontsize=15, weight="bold")
        coarse = False
        used_kinds = []
        for i, (key, contacts) in enumerate(selected):
            resname, number, chain, kind, indices = key
            angle = -math.pi / 2 + 2 * math.pi * i / max(len(selected), 1)
            x, y = 1100 + 890 * math.cos(angle), 900 + 640 * math.sin(angle)
            atom_level = bool(indices) and all(idx in atom_points for idx in indices)
            if atom_level:
                anchor_x = 600 + sum(atom_points[idx].x for idx in indices) / len(indices)
                anchor_y = 500 + sum(atom_points[idx].y for idx in indices) / len(indices)
            else:
                # Unverified endpoints terminate at the ligand region boundary,
                # never an arbitrarily chosen atom in the drawing.
                anchor_x, anchor_y = 1100 + 470 * math.cos(angle), 900 + 350 * math.sin(angle)
                coarse = True
            label, color, style = STYLES.get(kind, (kind, "#777777", "--"))
            used_kinds.append(kind)
            ax.plot([anchor_x, x], [anchor_y, y], linestyle=style, color=color, linewidth=1.4, zorder=3)
            distances = []
            for contact in contacts:
                try:
                    value = float(contact.get("distance_A"))
                    if math.isfinite(value) and value >= 0:
                        distances.append(f"{value:.2f}")
                except (TypeError, ValueError):
                    pass
                source_rows.append({**contact, "page": page + 1,
                                    "endpoint": "PLIP_ATOM_OR_GROUP" if atom_level else "LIGAND_LEVEL",
                                    "mapping_method": method})
            distances = list(dict.fromkeys(distances))
            distance_label = " / ".join(distances[:3]) + " Å" if distances else "distance unreported"
            if len(distances) > 3:
                distance_label += " (+ more in CSV)"
            residue_label = resname + number + (" · " + chain if chain else "")
            text = (residue_label or "Residue identity unreported") + "\n" + label
            if len(contacts) > 1:
                text += " ×" + str(len(contacts))
            text += "\n" + distance_label
            ax.text(x, y, text, ha="center", va="center", fontsize=8.5, color="#243047", zorder=4,
                    bbox={"boxstyle": "round,pad=.55", "facecolor": "white", "edgecolor": color, "linewidth": 1.4})
        if not selected:
            ax.text(1100, 1450, "PLIP recorded no residue contacts for this pose.", ha="center", fontsize=11)
        handles = [Line2D([0], [0], color=STYLES.get(kind, (kind, "#777", "--"))[1],
                          linestyle=STYLES.get(kind, (kind, "#777", "--"))[2],
                          label=STYLES.get(kind, (kind, "#777", "--"))[0]) for kind in dict.fromkeys(used_kinds)]
        if handles:
            ax.legend(handles=handles, loc="center", bbox_to_anchor=(.5, .085), ncol=min(4, len(handles)),
                      fontsize=8, frameon=False)
        note = "Ligand-level lines have no verified atom endpoint." if coarse else "Atom/group endpoints follow the saved PLIP mapping."
        ax.text(1100, 1690, note + " Distances are PLIP values; lines show contacts, not ligand covalent bonds.", ha="center", fontsize=8)
        ax.text(1100, 1740, "Exploratory docking; native validation not performed", ha="center", fontsize=9)
        name = row["compound_id"] + "_interactions_2d" + (f"_p{page + 1}" if page else "") + ".png"
        try:
            fig.savefig(destination / name, dpi=300)
        finally:
            plt.close(fig)
        outputs.append({"path": "figures/" + name,
                        "caption": row["compound_id"] + " 2D PLIP interaction diagram. " + note +
                        " " + method + ". Exploratory docking; native validation not performed."})
    import json
    (destination / (row["compound_id"] + "_interactions_2d.json")).write_text(
        json.dumps({"compound_id": row["compound_id"], "mapping_method": method,
                    "scientific_approval": False, "contacts": source_rows}, indent=2), encoding="utf-8")
    return outputs
