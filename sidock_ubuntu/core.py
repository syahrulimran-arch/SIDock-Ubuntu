from __future__ import annotations

import math
import re
import shutil
import subprocess
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from .model import DockingParameters, GridBox


def resolve_executable(name: str, label: str) -> str:
    path = shutil.which(name)
    if not path:
        raise RuntimeError(f"{label} executable not found on PATH: {name}")
    return path


def run_command(cmd: list[str], cwd: Path, log_path: Path) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        p = subprocess.run(cmd, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, text=True)
    return p.returncode


def box_from_reference(reference_pdb: Path, padding: float = 8.0, spacing: float = 0.375) -> GridBox:
    xyz = []
    for line in reference_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                pass
    if not xyz:
        raise ValueError("Reference ligand contains no readable coordinates.")
    mins = [min(v[i] for v in xyz) for i in range(3)]
    maxs = [max(v[i] for v in xyz) for i in range(3)]
    centers = [(mins[i] + maxs[i]) / 2 for i in range(3)]
    npts = []
    for i in range(3):
        span = (maxs[i] - mins[i]) + 2 * padding
        n = max(20, min(126, math.ceil(span / spacing)))
        if n % 2:
            n += 1
        npts.append(min(n, 126))
    return GridBox(centers[0], centers[1], centers[2], npts[0], npts[1], npts[2], spacing)


def clean_receptor_pdb(source: Path, destination: Path, ref_resname: str, ref_chain: str, ref_resseq: str) -> Path:
    kept = []
    for line in source.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("HETATM") and line[17:20].strip().upper() in {"HOH", "WAT", "H2O", "DOD"}:
            continue
        same_ref = line.startswith("HETATM") and line[17:20].strip().upper() == ref_resname and line[21:22].strip() == ref_chain and line[22:27].strip() == ref_resseq
        if same_ref:
            continue
        kept.append(line)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(kept) + "\n", encoding="utf-8")
    return destination


def prepare_receptor(source: Path, destination: Path, log_path: Path) -> Path:
    exe = shutil.which("mk_prepare_receptor.py") or shutil.which("mk_prepare_receptor")
    if not exe:
        raise RuntimeError("Meeko receptor preparation command was not found (mk_prepare_receptor.py).")
    base = destination.with_suffix("")
    code = run_command([exe, "-i", str(source.resolve()), "-o", str(base.resolve()), "-p"], destination.parent, log_path)
    candidates = [base.with_suffix(".pdbqt"), Path(str(base) + ".pdbqt")]
    out = next((p for p in candidates if p.exists()), None)
    if code != 0 or out is None:
        raise RuntimeError(f"Receptor preparation failed; see {log_path}")
    return out


def split_and_minimize_sdf(source: Path, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    supplier = Chem.SDMolSupplier(str(source), removeHs=False)
    outputs = []
    for idx, mol in enumerate(supplier, start=1):
        if mol is None:
            continue
        name = mol.GetProp("Compound_ID") if mol.HasProp("Compound_ID") else (mol.GetProp("_Name") if mol.HasProp("_Name") else f"ligand_{idx}")
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip()) or f"ligand_{idx}"
        mol = Chem.AddHs(mol)
        if mol.GetNumConformers() == 0:
            params = AllChem.ETKDGv3(); params.randomSeed = 0xC0FFEE
            if AllChem.EmbedMolecule(mol, params) != 0:
                continue
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol):
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            else:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            pass
        out = out_dir / f"{safe}.sdf"
        w = Chem.SDWriter(str(out)); w.write(mol); w.close(); outputs.append(out)
    return outputs


def prepare_ligand(source: Path, destination: Path, log_path: Path) -> Path:
    exe = shutil.which("mk_prepare_ligand.py") or shutil.which("mk_prepare_ligand")
    if not exe:
        raise RuntimeError("Meeko ligand preparation command was not found (mk_prepare_ligand.py).")
    code = run_command([exe, "-i", str(source.resolve()), "-o", str(destination.resolve())], destination.parent, log_path)
    if code != 0 or not destination.exists():
        raise RuntimeError(f"Ligand preparation failed for {source.name}; see {log_path}")
    return destination


def pdbqt_atom_types(path: Path) -> list[str]:
    types = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            bits = line.split()
            if bits:
                types.add(bits[-1])
    return sorted(types)


def write_gpf(receptor: Path, ligand_types: list[str], grid: GridBox, destination: Path) -> None:
    rec_types = " ".join(pdbqt_atom_types(receptor))
    lig_types = " ".join(ligand_types)
    lines = [
        f"npts {grid.npts_x} {grid.npts_y} {grid.npts_z}",
        f"gridfld {receptor.stem}.maps.fld",
        f"spacing {grid.spacing}",
        f"receptor_types {rec_types}",
        f"ligand_types {lig_types}",
        f"receptor {receptor.name}",
        f"gridcenter {grid.center_x:.3f} {grid.center_y:.3f} {grid.center_z:.3f}",
        "smooth 0.5",
    ]
    for t in ligand_types:
        lines.append(f"map {receptor.stem}.{t}.map")
    lines += [f"elecmap {receptor.stem}.e.map", f"dsolvmap {receptor.stem}.d.map"]
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_dpf(receptor_stem: str, ligand: Path, ligand_types: list[str], params: DockingParameters, destination: Path) -> None:
    lines = [
        "autodock_parameter_version 4.2", "outlev 1", "intelec", "seed pid time",
        f"ligand_types {' '.join(ligand_types)}", f"fld {receptor_stem}.maps.fld",
    ]
    for t in ligand_types:
        lines.append(f"map {receptor_stem}.{t}.map")
    lines += [f"elecmap {receptor_stem}.e.map", f"desolvmap {receptor_stem}.d.map", f"move {ligand.name}",
              "about 0.0 0.0 0.0", "tran0 random", "quaternion0 random", "dihe0 random", "torsdof 0",
              "rmstol %.2f" % params.rmstol, "ga_pop_size %d" % params.ga_pop_size,
              "ga_num_evals %d" % params.ga_num_evals, "ga_num_generations %d" % params.ga_num_generations,
              "ga_run %d" % params.ga_runs, "analysis"]
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_best_energy(dlg: Path) -> float | None:
    text = dlg.read_text(encoding="utf-8", errors="replace")
    values = []
    for pat in (r"Estimated Free Energy of Binding\s*=\s*([+-]?\d+(?:\.\d+)?)", r"DOCKED: USER\s+Estimated Free Energy of Binding\s*=\s*([+-]?\d+(?:\.\d+)?)"):
        values.extend(float(v) for v in re.findall(pat, text))
    return min(values) if values else None
