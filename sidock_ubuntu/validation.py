from __future__ import annotations

import csv
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem import Descriptors, rdMolDescriptors

from .core import pdbqt_atom_types, run_command, write_dpf
from .model import DockingParameters, GridBox


def crystal_ligand_to_sdf(pdb: Path, destination: Path) -> Chem.Mol:
    mol = Chem.MolFromPDBFile(str(pdb), removeHs=False, sanitize=False)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("Could not parse crystallographic ligand coordinates")
    try:
        from rdkit.Chem import rdDetermineBonds
        if mol.GetNumBonds() == 0:
            rdDetermineBonds.DetermineConnectivity(mol)
        Chem.SanitizeMol(mol)
    except Exception as exc:
        raise ValueError(f"Could not assign chemically valid ligand connectivity: {exc}") from exc
    mol.SetProp("_Name", "native_ligand_crystal")
    destination.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(destination)); writer.write(mol); writer.close()
    return mol


def prepare_native_ligand(native_pdb: Path, output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    crystal_pdb = output_dir / "native_ligand_crystal.pdb"; shutil.copy2(native_pdb, crystal_pdb)
    crystal_sdf = output_dir / "native_ligand_crystal.sdf"
    crystal = crystal_ligand_to_sdf(crystal_pdb, crystal_sdf)
    prepared_sdf = output_dir / "native_ligand_prepared.sdf"
    mol = Chem.AddHs(Chem.Mol(crystal))
    if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
        mol.RemoveAllConformers(); p = AllChem.ETKDGv3(); p.randomSeed = 0xC0FFEE; AllChem.EmbedMolecule(mol, p)
    if AllChem.MMFFHasAllMoleculeParams(mol): AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    w = Chem.SDWriter(str(prepared_sdf)); w.write(mol); w.close()
    native_pdbqt = output_dir / "native_ligand.pdbqt"
    exe = shutil.which("mk_prepare_ligand.py") or shutil.which("mk_prepare_ligand")
    if not exe: raise RuntimeError("Meeko ligand preparation command not found")
    prep_log = output_dir / "native_ligand_prepare.log"
    if run_command([exe, "-i", str(prepared_sdf.resolve()), "-o", str(native_pdbqt.resolve())], output_dir, prep_log) != 0:
        raise RuntimeError(f"Native ligand preparation failed; see {prep_log}")
    return {"crystal": crystal, "prepared": mol, "crystal_pdb": crystal_pdb, "crystal_sdf": crystal_sdf, "prepared_sdf": prepared_sdf, "pdbqt": native_pdbqt}


def _pose_mol_from_pdbqt(native: Chem.Mol, pose: Path) -> Chem.Mol:
    mol = Chem.Mol(native)
    conf = Chem.Conformer(mol.GetNumAtoms())
    atoms = []
    for line in pose.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atom_type=line[77:79].strip() or line.split()[-1]
            if atom_type in {"H", "HD", "HS"} or line[12:16].strip().startswith("H"):
                continue
            atoms.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if len(atoms) != mol.GetNumAtoms():
        raise ValueError("Redocked pose atom count does not match native ligand")
    for idx, xyz in enumerate(atoms): conf.SetAtomPosition(idx, xyz)
    mol.RemoveAllConformers(); mol.AddConformer(conf, assignId=True)
    return mol


def heavy_atom_rmsd(crystal: Chem.Mol, pose: Chem.Mol) -> float:
    crystal_h = Chem.RemoveHs(crystal)
    pose_h = Chem.RemoveHs(pose)
    mcs = rdFMCS.FindMCS([crystal_h, pose_h], atomCompare=rdFMCS.AtomCompare.CompareElements, bondCompare=rdFMCS.BondCompare.CompareOrder, ringMatchesRingOnly=True, completeRingsOnly=True, timeout=10)
    if mcs.canceled or not mcs.smartsString:
        raise ValueError("No chemically valid heavy-atom correspondence found")
    query = Chem.MolFromSmarts(mcs.smartsString)
    a = crystal_h.GetSubstructMatch(query); b = pose_h.GetSubstructMatch(query)
    if len(a) != crystal_h.GetNumAtoms() or len(b) != pose_h.GetNumAtoms():
        raise ValueError("Atom mapping does not cover all heavy atoms")
    return float(AllChem.GetBestRMS(crystal_h, pose_h, map=[list(zip(a, b))]))


def _atom_name(atom) -> str:
    info=atom.GetPDBResidueInfo()
    return info.GetName().strip() if info else f"atom_{atom.GetIdx()+1}"


def _write_mapping_audit(crystal: Chem.Mol, prepared: Chem.Mol, pose: Chem.Mol, destination: Path) -> int:
    c=Chem.RemoveHs(crystal); p=Chem.RemoveHs(pose)
    mcs=rdFMCS.FindMCS([c,p], atomCompare=rdFMCS.AtomCompare.CompareElements, bondCompare=rdFMCS.BondCompare.CompareOrder, ringMatchesRingOnly=True, completeRingsOnly=True, timeout=10)
    query=Chem.MolFromSmarts(mcs.smartsString); ca=c.GetSubstructMatch(query); pa=p.GetSubstructMatch(query)
    mapped=dict(zip(ca,pa)); fields=['crystal_atom_index','crystal_atom_name','element','prepared_atom_index','redocked_atom_index','mapping_method','included_in_rmsd','reason_if_excluded']
    with destination.open('w',newline='',encoding='utf-8') as fh:
        w=csv.DictWriter(fh,fieldnames=fields); w.writeheader()
        for idx, atom in enumerate(c.GetAtoms()):
            included=idx in mapped
            w.writerow({'crystal_atom_index':idx,'crystal_atom_name':_atom_name(atom),'element':atom.GetSymbol(),'prepared_atom_index':idx if idx < prepared.GetNumAtoms() else '','redocked_atom_index':mapped.get(idx,''),'mapping_method':'RDKit MCS graph/element/bond-order mapping' if included else 'unmapped','included_in_rmsd':included,'reason_if_excluded':'' if included else 'No unique chemically valid mapping'})
    return len(mapped)


def _write_preparation_audit(crystal: Chem.Mol, prepared: Chem.Mol, destination: Path) -> None:
    crystal_h=Chem.AddHs(Chem.Mol(crystal)); prepared_charge=Chem.GetFormalCharge(prepared); original_charge=Chem.GetFormalCharge(crystal)
    uncertainty=[]
    if original_charge == 0 and any(a.GetIsHeteroAtom() for a in Chem.RemoveHs(crystal).GetAtoms()):
        uncertainty.append('Formal charge/protonation state was not experimentally assigned; review before production docking.')
    if any(a.GetIsAromatic() for a in Chem.RemoveHs(crystal).GetAtoms()):
        uncertainty.append('No tautomer enumeration was performed; tautomer-dependent docking should be reviewed.')
    audit={'original_formal_charge':original_charge,'prepared_formal_charge':prepared_charge,'hydrogen_count':{'original':sum(a.GetSymbol()=='H' for a in crystal_h.GetAtoms()),'prepared':sum(a.GetSymbol()=='H' for a in prepared.GetAtoms())},'rotatable_bonds':int(rdMolDescriptors.CalcNumRotatableBonds(Chem.RemoveHs(prepared))),'aromatic_atom_count':sum(a.GetIsAromatic() for a in Chem.RemoveHs(prepared).GetAtoms()),'stereocenters':len(Chem.FindMolChiralCenters(Chem.RemoveHs(prepared),includeUnassigned=True)),'bond_order_changes':False,'protonation_changes':original_charge != prepared_charge,'disconnected_fragments_removed':False,'sanitization_warnings':[],'preparation_warnings':uncertainty,'formal_charge_uncertain':bool(uncertainty and original_charge == 0),'multiple_plausible_protonation_states':'review required' if uncertainty else False,'tautomer_warning':any('tautomer' in x for x in uncertainty)}
    destination.write_text(json.dumps(audit,indent=2),encoding='utf-8')


def _overlay(crystal: Chem.Mol, pose: Chem.Mol, rmsd: float, destination: Path) -> None:
    import matplotlib.pyplot as plt
    c = crystal.GetConformer(); p = pose.GetConformer()
    fig = plt.figure(figsize=(6, 5)); ax = fig.add_subplot(111, projection="3d")
    cx = [c.GetAtomPosition(i).x for i in range(crystal.GetNumAtoms())]; cy = [c.GetAtomPosition(i).y for i in range(crystal.GetNumAtoms())]; cz = [c.GetAtomPosition(i).z for i in range(crystal.GetNumAtoms())]
    px = [p.GetAtomPosition(i).x for i in range(pose.GetNumAtoms())]; py = [p.GetAtomPosition(i).y for i in range(pose.GetNumAtoms())]; pz = [p.GetAtomPosition(i).z for i in range(pose.GetNumAtoms())]
    ax.scatter(cx, cy, cz, c="royalblue", label="crystal", s=35); ax.scatter(px, py, pz, c="crimson", label="redocked", s=35)
    for i in range(len(cx)): ax.plot([cx[i], px[i]], [cy[i], py[i]], [cz[i], pz[i]], color="gray", linewidth=.7)
    ax.set_title(f"Native redocking overlay — heavy-atom RMSD {rmsd:.3f} Å"); ax.legend(); destination.parent.mkdir(parents=True, exist_ok=True); fig.tight_layout(); fig.savefig(destination, dpi=220); plt.close(fig)


def run_native_redocking(*, native_pdb: Path, receptor_pdbqt: Path, grid_dir: Path, work_dir: Path, output_dir: Path, grid: GridBox, params: DockingParameters, reference_info: dict[str, str], threshold: float = 2.0) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True); work_dir.mkdir(parents=True, exist_ok=True)
    native = prepare_native_ligand(native_pdb, output_dir)
    crystal = native["crystal"]; crystal_pdb = native["crystal_pdb"]; crystal_sdf = native["crystal_sdf"]; prepared_sdf = native["prepared_sdf"]; native_pdbqt = native["pdbqt"]
    local = work_dir / native_pdbqt.name; shutil.copy2(native_pdbqt, local)
    receptor = grid_dir / "receptor.pdbqt"; shutil.copy2(receptor, work_dir / receptor.name)
    for f in grid_dir.glob("receptor.*.map*"): shutil.copy2(f, work_dir / f.name)
    fld = grid_dir / "receptor.maps.fld"
    if fld.exists(): shutil.copy2(fld, work_dir / fld.name)
    dpf = work_dir / "native_redocking.dpf"; write_dpf("receptor", local, pdbqt_atom_types(local), params, dpf)
    dlg = output_dir / "native_redocking.dlg"
    autodock = shutil.which("autodock4")
    if not autodock or run_command([autodock, "-p", dpf.name, "-l", str(dlg.resolve())], work_dir, output_dir / "native_redocking.log") != 0:
        raise RuntimeError("Native redocking failed; see validation/native_redocking.log")
    from .analysis import _extract_models
    models = _extract_models(dlg)
    best = min((m for m in models if m[1] is not None), key=lambda m: float(m[1]), default=None)
    if best is None: raise RuntimeError("Native redocking produced no scored pose")
    pose_pdbqt = output_dir / "native_redocked_best.pdbqt"; pose_pdbqt.write_text("MODEL 1\n" + "\n".join(best[2]) + "\nENDMDL\n", encoding="utf-8")
    pose_mol = _pose_mol_from_pdbqt(crystal, pose_pdbqt)
    pose_pdb = output_dir / "native_redocked_best.pdb"
    Chem.MolToPDBFile(pose_mol, str(pose_pdb))
    rmsd = heavy_atom_rmsd(crystal, pose_mol)
    mapped_count=_write_mapping_audit(crystal, native['prepared'], pose_mol, output_dir/'atom_mapping.csv')
    _write_preparation_audit(crystal, native['prepared'], output_dir/'ligand_preparation_audit.json')
    overlay = output_dir / "redocking_overlay.png"; _overlay(crystal, pose_mol, rmsd, overlay)
    result = {"native_ligand": reference_info, "grid": grid.__dict__, "autodock_parameters": params.__dict__, "best_docking_energy_kcal_mol": float(best[1]), "rmsd_A": rmsd, "chemically_mapped_atom_count": mapped_count, "excluded_atoms": [], "validation_threshold_A": threshold, "status": "PASS" if rmsd <= threshold else "FAIL", "files": {"crystal": str(crystal_pdb), "crystal_sdf": str(crystal_sdf), "prepared_sdf": str(prepared_sdf), "pdbqt": str(native_pdbqt), "dlg": str(dlg), "best_pdbqt": str(pose_pdbqt), "best_pdb": str(pose_pdb), "overlay": str(overlay), "atom_mapping": str(output_dir/'atom_mapping.csv'), "preparation_audit": str(output_dir/'ligand_preparation_audit.json')}}
    (output_dir / "redocking_validation.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    with (output_dir / "redocking_validation.csv").open("w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerows([["status", "rmsd_A", "threshold_A", "best_energy_kcal_mol"], [result["status"], f"{rmsd:.3f}", threshold, result["best_docking_energy_kcal_mol"]]])
    return result
