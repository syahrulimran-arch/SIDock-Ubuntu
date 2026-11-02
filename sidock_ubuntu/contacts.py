from __future__ import annotations
import csv
import math
from pathlib import Path


def _atoms(path: Path, protein: bool) -> list[dict[str, object]]:
    out=[]
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")): continue
        typ=line[77:79].strip() or line.split()[-1]
        element={"A":"C","C":"C","OA":"O","NA":"N","HD":"H","N":"N","O":"O","S":"S","SA":"S","Cl":"CL","ZN":"ZN","MG":"MG","CA":"CA","FE":"FE","MN":"MN","CU":"CU"}.get(typ, typ[:1].upper())
        out.append({"atom":line[12:16].strip(),"element":element,"type":typ,"x":float(line[30:38]),"y":float(line[38:46]),"z":float(line[46:54]),"chain":line[21:22].strip(),"resname":line[17:20].strip(),"resnum":line[22:27].strip(),"charge":float(line[70:76]) if line[70:76].strip() else 0.0,"protein":protein})
    return out


def _distance(a,b): return math.sqrt(sum((a[k]-b[k])**2 for k in ("x","y","z")))


def _centroid(atoms):
    return {k: sum(a[k] for a in atoms) / len(atoms) for k in ("x", "y", "z")}


def write_contacts(rows: list[dict[str, object]], destination: Path) -> None:
    fields=["ligand","pose","interaction_type","protein_chain","residue_name","residue_number","protein_atom","ligand_atom","distance_A","angle_deg"]
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=fields); w.writeheader(); w.writerows(rows)


def detect_contacts(receptor: Path, ligand: Path, ligand_name: str, destination: Path | None = None) -> list[dict[str, object]]:
    protein=_atoms(receptor, True); lig=_atoms(ligand, False); rows=[]; seen=set()
    for p in protein:
        for l in lig:
            d=_distance(p,l); pt=p["type"]; lt=l["type"]
            kind=None
            if p["element"] in {"C","S"} and l["element"] in {"C","S"} and d <= 4.0: kind="hydrophobic"
            elif p["element"] in {"N","O","S"} and l["element"] in {"N","O","S"} and d <= 3.5 and ("H" in str(pt) or "O" in str(pt) or "N" in str(pt)) : kind="hydrogen_bond"
            elif p["charge"] * l["charge"] < 0 and d <= 4.0: kind="salt_bridge"
            if kind:
                key=(p["atom"],l["atom"],kind)
                if key not in seen:
                    seen.add(key); rows.append({"ligand":ligand_name,"pose":ligand.name,"interaction_type":kind,"protein_chain":p["chain"],"residue_name":p["resname"],"residue_number":p["resnum"],"protein_atom":p["atom"],"ligand_atom":l["atom"],"distance_A":f"{d:.3f}","angle_deg":""})
    ligand_ring = [a for a in lig if a["type"] == "A"]
    protein_rings = {}
    for a in protein:
        if a["type"] == "A": protein_rings.setdefault((a["chain"], a["resname"], a["resnum"]), []).append(a)
    if len(ligand_ring) >= 5:
        lc = _centroid(ligand_ring)
        for key, ring in protein_rings.items():
            if len(ring) < 5: continue
            d = _distance(lc, _centroid(ring))
            if 3.0 <= d <= 6.0:
                p = ring[0]; l = ligand_ring[0]
                rows.append({"ligand":ligand_name,"pose":ligand.name,"interaction_type":"pi-pi","protein_chain":key[0],"residue_name":key[1],"residue_number":key[2],"protein_atom":p["atom"],"ligand_atom":l["atom"],"distance_A":f"{d:.3f}","angle_deg":""})
        for p in protein:
            if p["element"] == "N" and p["charge"] > 0:
                d = _distance(lc, p)
                if d <= 6.0:
                    rows.append({"ligand":ligand_name,"pose":ligand.name,"interaction_type":"pi-cation","protein_chain":p["chain"],"residue_name":p["resname"],"residue_number":p["resnum"],"protein_atom":p["atom"],"ligand_atom":ligand_ring[0]["atom"],"distance_A":f"{d:.3f}","angle_deg":""})
    for p in protein:
        for l in lig:
            d = _distance(p, l)
            if l["element"] in {"F", "CL", "BR", "I"} and p["element"] in {"N", "O", "S"} and d <= 3.5:
                rows.append({"ligand":ligand_name,"pose":ligand.name,"interaction_type":"halogen_bond","protein_chain":p["chain"],"residue_name":p["resname"],"residue_number":p["resnum"],"protein_atom":p["atom"],"ligand_atom":l["atom"],"distance_A":f"{d:.3f}","angle_deg":""})
            if p["element"] in {"ZN", "MG", "CA", "FE", "MN", "CU"} and l["element"] in {"N", "O", "S"} and d <= 3.0:
                rows.append({"ligand":ligand_name,"pose":ligand.name,"interaction_type":"metal_coordination","protein_chain":p["chain"],"residue_name":p["resname"],"residue_number":p["resnum"],"protein_atom":p["atom"],"ligand_atom":l["atom"],"distance_A":f"{d:.3f}","angle_deg":""})
    if destination is not None:
        write_contacts(rows, destination)
    return rows
