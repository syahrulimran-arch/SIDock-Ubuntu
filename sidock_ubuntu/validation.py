from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem import Descriptors, rdMolDescriptors

from .core import pdbqt_atom_types, prepare_ligand, run_command, write_dpf
from .model import DockingParameters, GridBox
from .reference import load_local_ccd


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


def reconstruct_reference_from_ccd(native_pdb: Path, output_dir: Path, component_id: str | None = None, ccd_path: Path | None = None) -> dict[str, object]:
    """Build the reference ligand graph from the RCSB CCD and PDB coordinates."""
    try:
        import gemmi
    except ImportError as exc:
        raise RuntimeError("gemmi is required for authoritative CCD reference chemistry") from exc
    output_dir.mkdir(parents=True, exist_ok=True)
    pdb_lines = native_pdb.read_text(encoding="utf-8", errors="replace").splitlines()
    atoms_by_name = {line[12:16].strip(): line for line in pdb_lines if line.startswith(("ATOM", "HETATM"))}
    if not atoms_by_name:
        raise RuntimeError("REFERENCE_CHEMISTRY_UNRESOLVED: reference ligand has no atoms")
    component_id = component_id or next((line[17:20].strip() for line in pdb_lines if line.startswith("HETATM")), None)
    if not component_id:
        raise RuntimeError("REFERENCE_CHEMISTRY_UNRESOLVED: component ID was not found")
    local_json=Path(ccd_path) if ccd_path else output_dir.parent.parent/'input'/'ccd.json'
    if local_json.exists():
        source=load_local_ccd(local_json,component_id)
        if source.get('parse_status')!='PASS': raise RuntimeError(f"REFERENCE_CHEMISTRY_UNRESOLVED: local CCD {source.get('warnings')}")
        chemistry=source['chemistry']; atoms=chemistry.get('atoms',[]); rw=Chem.RWMol(); ids=[]
        for item in atoms:
            if str(item.get('element','')).upper() in {'H','D'}: continue
            atom=Chem.Atom(str(item.get('element','C'))); atom.SetFormalCharge(int(item.get('formal_charge',item.get('charge',0)) or 0)); atom.SetIsAromatic(bool(item.get('aromatic',False))); atom.SetProp('_CCDAtomId',str(item['id'])); rw.AddAtom(atom); ids.append(str(item['id']))
        index={x:i for i,x in enumerate(ids)}
        for bond in chemistry.get('bonds',[]):
            a,b=str(bond.get('a',bond.get('atom1'))),str(bond.get('b',bond.get('atom2'))); order=float(bond.get('order',1.0)); bt=Chem.BondType.AROMATIC if bond.get('aromatic') or order==1.5 else (Chem.BondType.DOUBLE if order==2 else Chem.BondType.TRIPLE if order==3 else Chem.BondType.SINGLE); rw.AddBond(index[a],index[b],bt)
        mol=rw.GetMol(); missing=[x for x in ids if x not in atoms_by_name]
        if missing: raise RuntimeError(f"REFERENCE_CHEMISTRY_UNRESOLVED: CCD atoms missing from crystallographic instance: {missing}")
        conf=Chem.Conformer(mol.GetNumAtoms())
        for i,atom_id in enumerate(ids):
            line=atoms_by_name[atom_id]; conf.SetAtomPosition(i,(float(line[30:38]),float(line[38:46]),float(line[46:54])))
        mol.AddConformer(conf,assignId=True); Chem.SanitizeMol(mol)
        return {'mol':mol,'component_id':component_id,'cif':local_json,'sdf':output_dir/f'{component_id}_ccd.sdf','metadata':{'component_id':component_id,'source':'local_json','source_url':str(local_json),'sha256':source.get('source_hash')},'mapping':{x:x for x in ids}}
    cif_path = output_dir / f"{component_id}_ccd.cif"
    url = f"https://files.rcsb.org/ligands/view/{component_id}.cif"
    if not cif_path.exists():
        try:
            urllib.request.urlretrieve(url, cif_path)
        except Exception as exc:
            raise RuntimeError(f"REFERENCE_CHEMISTRY_UNRESOLVED: could not retrieve {url}: {exc}") from exc
    block = gemmi.cif.read(str(cif_path)).sole_block()
    atom_cat = block.get_mmcif_category("_chem_comp_atom.")
    bond_cat = block.get_mmcif_category("_chem_comp_bond.")
    atom_ids = [atom_id for atom_id, element in zip(atom_cat["atom_id"], atom_cat["type_symbol"]) if element.upper() != "H"]
    if not atom_ids:
        raise RuntimeError("REFERENCE_CHEMISTRY_UNRESOLVED: CCD has no atom table")
    rw = Chem.RWMol()
    for atom_id, element, charge, aromatic in zip(atom_cat["atom_id"], atom_cat["type_symbol"], atom_cat["charge"], atom_cat["pdbx_aromatic_flag"]):
        if element.upper() == "H":
            continue
        element = element[:1].upper() + element[1:].lower()
        atom = Chem.Atom(element)
        if charge not in ("?", "."):
            atom.SetFormalCharge(int(charge))
        atom.SetIsAromatic(aromatic == "Y")
        atom.SetProp("_CCDAtomId", atom_id)
        rw.AddAtom(atom)
    order_map = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    id_to_idx = {atom_id: i for i, atom_id in enumerate(atom_ids)}
    for a, b, order, aromatic in zip(bond_cat["atom_id_1"], bond_cat["atom_id_2"], bond_cat["value_order"], bond_cat["pdbx_aromatic_flag"]):
        if a not in id_to_idx or b not in id_to_idx:
            continue
        if a not in id_to_idx or b not in id_to_idx:
            raise RuntimeError("REFERENCE_CHEMISTRY_UNRESOLVED: CCD bond references unknown atom")
        bond = rw.AddBond(id_to_idx[a], id_to_idx[b], order_map.get(order, Chem.BondType.SINGLE))
        if aromatic == "Y":
            rw.GetBondBetweenAtoms(id_to_idx[a], id_to_idx[b]).SetIsAromatic(True)
    mol = rw.GetMol()
    missing = [atom_id for atom_id in atom_ids if atom_id not in atoms_by_name]
    if missing:
        raise RuntimeError(f"REFERENCE_CHEMISTRY_UNRESOLVED: CCD atoms missing from crystallographic instance: {missing}")
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, atom_id in enumerate(atom_ids):
        line = atoms_by_name[atom_id]
        conf.SetAtomPosition(i, (float(line[30:38]), float(line[38:46]), float(line[46:54])))
    mol.AddConformer(conf, assignId=True)
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        raise RuntimeError(f"REFERENCE_CHEMISTRY_UNRESOLVED: CCD graph failed sanitization: {exc}") from exc
    sdf_path = output_dir / f"{component_id}_ccd.sdf"
    writer = Chem.SDWriter(str(sdf_path)); writer.write(mol); writer.close()
    digest = hashlib.sha256(cif_path.read_bytes()).hexdigest()
    metadata = {"component_id": component_id, "formula": block.find_value("_chem_comp.formula"), "formal_charge": block.find_value("_chem_comp.pdbx_formal_charge"), "source_url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat(), "sha256": digest, "atom_count": mol.GetNumAtoms(), "heavy_atom_count": Chem.RemoveHs(mol).GetNumAtoms()}
    (output_dir / "ccd_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return {"mol": mol, "component_id": component_id, "cif": cif_path, "sdf": sdf_path, "metadata": metadata, "mapping": {atom_id: atoms_by_name[atom_id] for atom_id in atom_ids}}


def prepare_native_ligand(native_pdb: Path, output_dir: Path, ccd_path: Path | None = None) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    crystal_pdb = output_dir / "native_ligand_crystal.pdb"; shutil.copy2(native_pdb, crystal_pdb)
    chemistry_dir = output_dir.parent.parent / "reference_chemistry"
    reconstructed = reconstruct_reference_from_ccd(crystal_pdb, chemistry_dir, ccd_path=ccd_path)
    crystal = reconstructed["mol"]
    mapping_path = output_dir / "ccd_coordinate_mapping.csv"
    with mapping_path.open("w", newline="", encoding="utf-8") as fh:
        fields = ["ccd_atom_id", "pdb_atom_name", "element", "ccd_bond_environment", "crystal_x", "crystal_y", "crystal_z", "mapping_method", "mapping_confidence"]
        writer = csv.DictWriter(fh, fieldnames=fields); writer.writeheader()
        conf = crystal.GetConformer()
        for i, atom in enumerate(crystal.GetAtoms()):
            atom_id = atom.GetProp("_CCDAtomId")
            pos = conf.GetAtomPosition(i)
            neighbors = ";".join(f"{n.GetProp('_CCDAtomId')}:{crystal.GetBondBetweenAtoms(i, n.GetIdx()).GetBondType()}" for n in atom.GetNeighbors())
            writer.writerow({"ccd_atom_id": atom_id, "pdb_atom_name": atom_id, "element": atom.GetSymbol(), "ccd_bond_environment": neighbors, "crystal_x": f"{pos.x:.3f}", "crystal_y": f"{pos.y:.3f}", "crystal_z": f"{pos.z:.3f}", "mapping_method": "authoritative CCD atom name + crystallographic atom name", "mapping_confidence": "high"})
    h = Chem.RemoveHs(crystal)
    chemistry_audit = {"component_id": reconstructed["component_id"], "formula": reconstructed["metadata"].get("formula"), "formal_charge": Chem.GetFormalCharge(crystal), "heavy_atom_count": h.GetNumAtoms(), "bond_count": h.GetNumBonds(), "single_bonds": sum(b.GetBondType() == Chem.BondType.SINGLE for b in h.GetBonds()), "double_bonds": sum(b.GetBondType() == Chem.BondType.DOUBLE for b in h.GetBonds()), "triple_bonds": sum(b.GetBondType() == Chem.BondType.TRIPLE for b in h.GetBonds()), "aromatic_bonds": sum(b.GetIsAromatic() for b in h.GetBonds()), "ring_count": int(rdMolDescriptors.CalcNumRings(h)), "rotatable_bonds": int(rdMolDescriptors.CalcNumRotatableBonds(h)), "stereochemistry": Chem.FindMolChiralCenters(h, includeUnassigned=True), "CCD_source": reconstructed["metadata"]["source_url"], "CCD_hash": reconstructed["metadata"]["sha256"], "connectivity_validation": "PASS"}
    (output_dir / "reference_chemistry_audit.json").write_text(json.dumps(chemistry_audit, indent=2), encoding="utf-8")
    crystal_sdf = output_dir / "native_ligand_crystal.sdf"
    writer = Chem.SDWriter(str(crystal_sdf)); writer.write(crystal); writer.close()
    prepared_sdf = output_dir / "native_ligand_prepared.sdf"
    # Native validation must preserve the experimentally observed heavy-atom
    # coordinates.  Study compounds are minimized elsewhere, but minimizing
    # the native ligand here changes the reference conformation before AD4
    # search and makes the validation conflate preparation with redocking.
    mol = Chem.AddHs(Chem.Mol(crystal), addCoords=True)
    if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
        raise ValueError("Native ligand has no 3D crystallographic coordinates")
    w = Chem.SDWriter(str(prepared_sdf)); w.write(mol); w.close()
    native_pdbqt = output_dir / "native_ligand.pdbqt"
    prep_log = output_dir / "native_ligand_prepare.log"
    prepare_ligand(prepared_sdf, native_pdbqt, prep_log)
    return {"crystal": crystal, "prepared": mol, "crystal_pdb": crystal_pdb, "crystal_sdf": crystal_sdf, "prepared_sdf": prepared_sdf, "pdbqt": native_pdbqt, "reference_chemistry": reconstructed}


def _pose_mol_from_pdbqt(native: Chem.Mol, pose: Path, reference_pdbqt: Path | None = None) -> Chem.Mol:
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
    if reference_pdbqt is None:
        for idx, xyz in enumerate(atoms): conf.SetAtomPosition(idx, xyz)
    else:
        ref_atoms = []
        for line in reference_pdbqt.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith(("ATOM", "HETATM")):
                atom_type = line[77:79].strip() or line.split()[-1]
                if atom_type in {"H", "HD", "HS"} or line[12:16].strip().startswith("H"):
                    continue
                ref_atoms.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        if len(ref_atoms) != len(atoms):
            raise ValueError("Reference PDBQT and pose have different heavy-atom counts")
        nconf = native.GetConformer(); unused = set(i for i in range(native.GetNumAtoms()) if native.GetAtomWithIdx(i).GetSymbol() != "H"); order = []
        for xyz in ref_atoms:
            candidates = sorted((math.dist(xyz, (nconf.GetAtomPosition(i).x, nconf.GetAtomPosition(i).y, nconf.GetAtomPosition(i).z)), i) for i in unused)
            if not candidates or candidates[0][0] > 0.25: raise ValueError("Could not coordinate-validate PDBQT-to-CCD atom mapping")
            idx = candidates[0][1]; order.append(idx); unused.remove(idx)
        for idx, xyz in zip(order, atoms): conf.SetAtomPosition(idx, xyz)
    mol.RemoveAllConformers(); mol.AddConformer(conf, assignId=True)
    return mol


def fixed_frame_symmetry_corrected_heavy_atom_rmsd(crystal: Chem.Mol, pose: Chem.Mol) -> float:
    """Direct Cartesian RMSD; symmetry changes atom correspondence only."""
    crystal_h = Chem.RemoveHs(crystal)
    pose_h = Chem.RemoveHs(pose)
    mcs = rdFMCS.FindMCS([crystal_h, pose_h], atomCompare=rdFMCS.AtomCompare.CompareElements, bondCompare=rdFMCS.BondCompare.CompareOrder, ringMatchesRingOnly=True, completeRingsOnly=True, timeout=10)
    if mcs.canceled or not mcs.smartsString:
        raise ValueError("No chemically valid heavy-atom correspondence found")
    query = Chem.MolFromSmarts(mcs.smartsString)
    crystal_matches = crystal_h.GetSubstructMatches(query, uniquify=False)
    pose_matches = pose_h.GetSubstructMatches(query, uniquify=False)
    if not crystal_matches or not pose_matches or len(crystal_matches[0]) != crystal_h.GetNumAtoms() or len(pose_matches[0]) != pose_h.GetNumAtoms():
        raise ValueError("Atom mapping does not cover all heavy atoms")
    cconf, pconf = crystal_h.GetConformer(), pose_h.GetConformer()
    best = float("inf")
    for a in crystal_matches:
        for b in pose_matches:
            value = math.sqrt(sum(sum((cconf.GetAtomPosition(i)[k] - pconf.GetAtomPosition(j)[k]) ** 2 for k in range(3)) for i, j in zip(a, b)) / len(a))
            best = min(best, value)
    return float(best)


def heavy_atom_rmsd(crystal: Chem.Mol, pose: Chem.Mol) -> float:
    """Backward-compatible name for the fixed-frame metric."""
    return fixed_frame_symmetry_corrected_heavy_atom_rmsd(crystal, pose)


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
    if original_charge == 0 and any(a.GetSymbol() not in {'C', 'H'} for a in Chem.RemoveHs(crystal).GetAtoms()):
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


def run_native_redocking(*, native_pdb: Path, receptor_pdbqt: Path, grid_dir: Path, work_dir: Path, output_dir: Path, grid: GridBox, params: DockingParameters, reference_info: dict[str, str], threshold: float = 2.0, ccd_path: Path | None = None) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True); work_dir.mkdir(parents=True, exist_ok=True)
    # Keep the authoritative project CCD provenance all the way through the
    # redocking helper.  This function is also called independently by the
    # production executor, so relying on the caller's earlier preparation
    # would allow a second, legacy reconstruction to replace the chemistry.
    native = prepare_native_ligand(native_pdb, output_dir, ccd_path=ccd_path)
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
    pose_mol = _pose_mol_from_pdbqt(crystal, pose_pdbqt, native_pdbqt)
    pose_pdb = output_dir / "native_redocked_best.pdb"
    Chem.MolToPDBFile(pose_mol, str(pose_pdb))
    rmsd = heavy_atom_rmsd(crystal, pose_mol)
    mapped_count=_write_mapping_audit(crystal, native['prepared'], pose_mol, output_dir/'atom_mapping.csv')
    _write_preparation_audit(crystal, native['prepared'], output_dir/'ligand_preparation_audit.json')
    overlay = output_dir / "redocking_overlay.png"; _overlay(crystal, pose_mol, rmsd, overlay)
    result = {"native_ligand": reference_info, "grid": grid.__dict__, "autodock_parameters": params.__dict__, "best_docking_energy_kcal_mol": float(best[1]), "rmsd_A": rmsd, "chemically_mapped_atom_count": mapped_count, "excluded_atoms": [], "validation_threshold_A": threshold, "status": "PASS" if rmsd <= threshold else "FAIL", "files": {"crystal": str(crystal_pdb), "crystal_sdf": str(crystal_sdf), "prepared_sdf": str(prepared_sdf), "pdbqt": str(native_pdbqt), "dlg": str(dlg), "best_pdbqt": str(pose_pdbqt), "best_pdb": str(pose_pdb), "overlay": str(overlay), "atom_mapping": str(output_dir/'atom_mapping.csv'), "preparation_audit": str(output_dir/'ligand_preparation_audit.json'), "reference_chemistry": {"cif": str(native['reference_chemistry']['cif']), "sdf": str(native['reference_chemistry']['sdf']), "metadata": str(native['reference_chemistry']['cif'].parent/'ccd_metadata.json')}}}
    (output_dir / "redocking_validation.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    with (output_dir / "redocking_validation.csv").open("w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerows([["status", "rmsd_A", "threshold_A", "best_energy_kcal_mol"], [result["status"], f"{rmsd:.3f}", threshold, result["best_docking_energy_kcal_mol"]]])
    return result
