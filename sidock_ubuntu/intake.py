from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import csv, math

WATER = {"HOH", "WAT", "H2O", "DOD"}
COMMON_IONS = {"NA", "K", "CL", "CA", "MG", "ZN", "MN", "FE", "CU", "CO", "NI", "CD", "HG", "BR", "IOD", "I", "LI", "CS", "SR", "BA", "AL", "AG", "AU", "PB"}
COMMON_BUFFER = {"SO4", "PO4", "NO3", "ACT", "ACE", "FMT", "GOL", "EDO", "PEG", "MPD", "TRS", "MES", "HEP", "BME", "DMS", "DTT", "NAG", "MAN"}
AMINO_HEAVY = {
 "ALA": {"N","CA","C","O","CB"}, "ARG": {"N","CA","C","O","CB","CG","CD","NE","CZ","NH1","NH2"},
 "ASN": {"N","CA","C","O","CB","CG","OD1","ND2"}, "ASP": {"N","CA","C","O","CB","CG","OD1","OD2"},
 "CYS": {"N","CA","C","O","CB","SG"}, "GLN": {"N","CA","C","O","CB","CG","CD","OE1","NE2"},
 "GLU": {"N","CA","C","O","CB","CG","CD","OE1","OE2"}, "GLY": {"N","CA","C","O"},
 "HIS": {"N","CA","C","O","CB","CG","ND1","CD2","CE1","NE2"}, "ILE": {"N","CA","C","O","CB","CG1","CG2","CD1"},
 "LEU": {"N","CA","C","O","CB","CG","CD1","CD2"}, "LYS": {"N","CA","C","O","CB","CG","CD","CE","NZ"},
 "MET": {"N","CA","C","O","CB","CG","SD","CE"}, "PHE": {"N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ"},
 "PRO": {"N","CA","C","O","CB","CG","CD"}, "SER": {"N","CA","C","O","CB","OG"},
 "THR": {"N","CA","C","O","CB","OG1","CG2"}, "TRP": {"N","CA","C","O","CB","CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"},
 "TYR": {"N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ","OH"}, "VAL": {"N","CA","C","O","CB","CG1","CG2"}}


@dataclass(frozen=True)
class HeteroGroup:
    resname: str
    chain: str
    resseq: str
    atoms: tuple[str, ...]

    @property
    def atom_count(self) -> int:
        return len(self.atoms)

    @property
    def key(self) -> str:
        return f"{self.resname}:{self.chain or '-'}:{self.resseq}"


def hetero_groups(pdb: Path) -> list[HeteroGroup]:
    grouped: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("HETATM"):
            continue
        key = (line[17:20].strip().upper(), line[21:22].strip(), line[22:27].strip())
        grouped[key].append(line)
    return [HeteroGroup(r, c, n, tuple(lines)) for (r, c, n), lines in grouped.items()]


def choose_reference_ligand(pdb: Path) -> HeteroGroup | None:
    candidates = []
    for group in hetero_groups(pdb):
        if group.resname in WATER or group.resname in COMMON_IONS or group.resname in COMMON_BUFFER:
            continue
        if group.atom_count < 6:
            continue
        candidates.append(group)
    return rank_reference_candidates(pdb, candidates)[0] if candidates else None

def _candidate_metrics(group: HeteroGroup, pdb_lines: list[str]) -> dict:
    vals=[]; contacts5=contacts8=0
    for line in group.atoms:
        try: vals.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except ValueError: pass
    for line in pdb_lines:
        if not line.startswith("ATOM"): continue
        try: xyz=(float(line[30:38]),float(line[38:46]),float(line[46:54]))
        except ValueError: continue
        if any(math.dist(xyz, v)<=5 for v in vals): contacts5 += 1
        if any(math.dist(xyz, v)<=8 for v in vals): contacts8 += 1
    occ=[]; bf=[]
    for line in group.atoms:
        try: occ.append(float(line[54:60])); bf.append(float(line[60:66]))
        except ValueError: pass
    atom_names=[line[12:16].strip() for line in group.atoms if line[76:78].strip().upper() not in {"H","D"}]
    return {"ligand_code":group.resname,"chain":group.chain,"residue_number":group.resseq,"atom_count":group.atom_count,"atom_names":atom_names,
            "heavy_atom_count":sum(1 for x in group.atoms if x[76:78].strip() not in {"H","D"}),
            "occupancy":sum(occ)/len(occ) if occ else None,"B_factor_summary":sum(bf)/len(bf) if bf else None,
            "ccd_id":group.resname,"mapping_fraction":None,"mapping_ambiguity":None,"missing_atoms":None,"extra_atoms":None,
            "alternate_location":any(line[16:17].strip() for line in group.atoms),"model_number":1,"RSCC":None,"RSR":None,
            "chemical_fingerprint":None,
            "CCD_available":"UNKNOWN","CCD_mapping_complete":"UNKNOWN","protein_contacts_5A":contacts5,
            "protein_contacts_8A":contacts8,"candidate_status":"CREDIBLE"}

def rank_reference_candidates(pdb: Path, candidates: list[HeteroGroup] | None = None) -> list[HeteroGroup]:
    candidates = candidates if candidates is not None else [g for g in hetero_groups(pdb) if g.resname not in WATER|COMMON_IONS|COMMON_BUFFER and g.atom_count >= 6]
    lines=pdb.read_text(encoding="utf-8",errors="replace").splitlines()
    metrics={g.key:_candidate_metrics(g,lines) for g in candidates}
    # Quality-first, deterministic tie-break; chemistry-equivalent copies are ranked by occupancy,
    # then complete atom count, then contact count, then stable identifier.
    return sorted(candidates, key=lambda g: (-metrics[g.key]["occupancy"] if metrics[g.key]["occupancy"] is not None else 0,
                                             -metrics[g.key]["heavy_atom_count"], -metrics[g.key]["protein_contacts_5A"], g.key))

def reference_candidates(pdb: Path) -> list[dict]:
    groups=[g for g in hetero_groups(pdb) if g.resname not in WATER|COMMON_IONS|COMMON_BUFFER and g.atom_count >= 6]
    ordered=rank_reference_candidates(pdb,groups); lines=pdb.read_text(encoding="utf-8",errors="replace").splitlines()
    return [{**_candidate_metrics(g,lines), "key":g.key, "rank":i+1, "recommended":i==0} for i,g in enumerate(ordered)]

def write_reference_candidates(pdb: Path, destination: Path) -> list[dict]:
    rows=reference_candidates(pdb); destination.parent.mkdir(parents=True,exist_ok=True)
    for row in rows:
        row["selection_reason"] = "recommended deterministic top candidate" if row["recommended"] else "equivalent candidate retained for review"
    fields=["ligand_code","chain","residue_number","atom_count","atom_names","heavy_atom_count","occupancy","B_factor_summary","ccd_id","mapping_fraction","mapping_ambiguity","missing_atoms","extra_atoms","alternate_location","model_number","RSCC","RSR","chemical_fingerprint","CCD_available","CCD_mapping_complete","protein_contacts_5A","protein_contacts_8A","candidate_status","key","rank","recommended","selection_reason"]
    with destination.open("w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=fields); w.writeheader(); w.writerows(rows)
    return rows

def receptor_completeness(pdb: Path, reference: HeteroGroup, *, peripheral_cutoff=8.0) -> list[dict]:
    residues={}
    lines=pdb.read_text(encoding="utf-8",errors="replace").splitlines()
    refxyz=[]
    for line in reference.atoms:
        try: refxyz.append((float(line[30:38]),float(line[38:46]),float(line[46:54])))
        except ValueError: pass
    for line in lines:
        if not line.startswith("ATOM"): continue
        key=(line[21:22].strip(),line[22:27].strip(),line[17:20].strip().upper())
        residues.setdefault(key,[]).append(line)
    rows=[]
    for (chain,num,name), atoms in residues.items():
        present={a[12:16].strip().upper() for a in atoms if a[76:78].strip().upper() not in {"H","D"}}
        expected=AMINO_HEAVY.get(name)
        if expected is None:
            classification="UNKNOWN_RESIDUE"; missing=set()
        else:
            missing=expected-present
            if not missing: classification="COMPLETE"
            else:
                coords=[]
                for a in atoms:
                    try: coords.append((float(a[30:38]),float(a[38:46]),float(a[46:54])))
                    except ValueError: pass
                distance=min((math.dist(x,y) for x in coords for y in refxyz),default=float('inf'))
                classification="INCOMPLETE_BINDING_SITE" if distance<=5 else ("INCOMPLETE_PERIPHERAL" if distance<=peripheral_cutoff else "INCOMPLETE_REMOTE")
        coords=[]
        for a in atoms:
            try: coords.append((float(a[30:38]),float(a[38:46]),float(a[46:54])))
            except ValueError: pass
        distance=min((math.dist(x,y) for x in coords for y in refxyz),default=float('inf'))
        rows.append({"chain":chain,"residue_name":name,"residue_number":num,"present_heavy_atoms":";".join(sorted(present)),"expected_heavy_atoms":";".join(sorted(expected or set())),"missing_heavy_atoms":";".join(sorted(missing)),"distance_to_reference_ligand":distance,"within_5A":distance<=5,"within_8A":distance<=8,"classification":classification})
    return rows

def write_receptor_completeness(rows: list[dict], destination: Path) -> None:
    destination.parent.mkdir(parents=True,exist_ok=True)
    fields=list(rows[0]) if rows else ["chain","residue_name","residue_number","present_heavy_atoms","expected_heavy_atoms","missing_heavy_atoms","distance_to_reference_ligand","within_5A","within_8A","classification"]
    with destination.open("w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=fields); w.writeheader(); w.writerows(rows)


def write_group(group: HeteroGroup, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(group.atoms) + "\nEND\n", encoding="utf-8")
    return destination


def protein_intake_summary(pdb: Path) -> dict[str, object]:
    from .structure_intake import classify_structure
    groups = hetero_groups(pdb)
    ref = choose_reference_ligand(pdb)
    semantic=classify_structure(pdb)
    return {
        "pdb": str(pdb.resolve()),
        "structure_classification": semantic.classification,
        "structure_reason": semantic.reason,
        "processed_indicators": list(semantic.processed_indicators),
        "reference_candidates": list(semantic.candidates),
        "ions": list(semantic.ions),
        "waters": list(semantic.waters),
        "hetero_groups": [{"key": g.key, "resname": g.resname, "chain": g.chain, "resseq": g.resseq, "atoms": g.atom_count} for g in groups],
        "selected_reference_ligand": ref.key if ref else None,
        "reference_ligand_atoms": ref.atom_count if ref else 0,
    }
