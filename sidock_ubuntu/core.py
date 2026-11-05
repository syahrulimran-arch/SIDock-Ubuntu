from __future__ import annotations

import math
import re
import shutil
import subprocess
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from .model import DockingParameters, GridBox
from .dependencies import resolve_meeko_receptor_preparer, resolve_meeko_script

SUPPORTED_AD4_TYPES = {'C','A','N','NA','OA','OS','HD','HS','S','SA','Cl','Br','I','F','P','Mg','Ca','Zn','Fe'}

def validate_ad4_types(types: list[str]) -> list[str]:
    ordered=sorted(set(types))
    bad=[x for x in ordered if x not in SUPPORTED_AD4_TYPES]
    if bad: raise RuntimeError(f'E501 UNSUPPORTED_AUTODOCK_TYPES: {", ".join(bad)}')
    return ordered


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
        if line.startswith("CONECT"):
            continue
        if line.startswith("HETATM") and line[17:20].strip().upper() in {"HOH", "WAT", "H2O", "DOD"}:
            continue
        # Remove all crystallographic copies of the reference ligand from the
        # receptor; the selected copy is retained separately as the native
        # validation ligand.  Keeping sibling ligand copies would introduce
        # artificial clashes in a biological oligomer.
        same_ref = line.startswith("HETATM") and line[17:20].strip().upper() == ref_resname
        if same_ref:
            continue
        kept.append(line)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(kept) + "\n", encoding="utf-8")
    return destination


def prepare_receptor(source: Path, destination: Path, log_path: Path, *, allow_bad_res: bool = False) -> Path:
    resolved = resolve_meeko_receptor_preparer(); exe = resolved.get('path')
    if not resolved['available']:
        raise RuntimeError("Meeko receptor preparation command was not found (mk_prepare_receptor.py).")
    base = destination.with_suffix("")
    cmd = [exe, "--read_pdb", str(source.resolve()), "--default_altloc", "A", "-o", str(base.resolve()), "-p"]
    if allow_bad_res:
        cmd.append("--allow_bad_res")
    code = run_command(cmd, destination.parent, log_path)
    candidates = [base.with_suffix(".pdbqt"), Path(str(base) + ".pdbqt")]
    out = next((p for p in candidates if p.exists()), None)
    if code != 0 or out is None:
        from .receptor_diagnostics import analyze_meeko_log
        log = log_path.read_text(encoding="utf-8", errors="replace") if log_path.is_file() else ""
        detail = analyze_meeko_log(log)
        reasons = []
        if detail["unknown_components"]:
            reasons.append("unknown components: " + ", ".join(detail["unknown_components"]))
        if detail["template_mismatches"]:
            reasons.append("unmatched residues: " + ", ".join(detail["template_mismatches"]))
        if not reasons:
            reasons.append("Meeko exited without a valid PDBQT; inspect its output")
        raise RuntimeError(
            f"Receptor preparation failed (Meeko exit={code}; {'; '.join(reasons)}). "
            f"Open Step 2 for the detailed explanation or see {log_path}"
        )
    return out


def split_and_minimize_sdf(source: Path, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    supplier = Chem.SDMolSupplier(str(source), removeHs=False)
    outputs = []
    for idx, mol in enumerate(supplier, start=1):
        if mol is None:
            continue
        name = mol.GetProp("Compound_ID") if mol.HasProp("Compound_ID") else (mol.GetProp("_Name") if mol.HasProp("_Name") else "")
        if not name.strip():
            continue
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip()) or f"ligand_{idx}"
        mol = Chem.AddHs(mol)
        if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
            if mol.GetNumConformers():
                mol.RemoveAllConformers()
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
    exe = resolve_meeko_script('mk_prepare_ligand.py') or resolve_meeko_script('mk_prepare_ligand')
    if not exe:
        raise RuntimeError("Meeko ligand preparation command was not found (mk_prepare_ligand.py).")
    # Step 3 is allowed to run before receptor preparation. Its PDBQT
    # output folder must therefore be created explicitly on first use.
    destination.parent.mkdir(parents=True, exist_ok=True)
    code = run_command([exe, "-i", str(source.resolve()), "-o", str(destination.resolve())], destination.parent, log_path)
    if code != 0 or not destination.is_file() or destination.stat().st_size == 0:
        raise RuntimeError(
            f"Ligand preparation failed for {source.name} (Meeko exit={code}); "
            f"inspect the unmodified source and log: {log_path}"
        )
    if not any(line.startswith(("ATOM", "HETATM")) for line in
               destination.read_text(encoding="utf-8", errors="replace").splitlines()):
        raise RuntimeError(f"Invalid ligand PDBQT for {source.name}: no atoms; see {log_path}")
    return destination


def pdbqt_atom_types(path: Path) -> list[str]:
    types = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            bits = line.split()
            if bits:
                types.add(bits[-1])
    return sorted(types)


def pdbqt_center(path: Path) -> tuple[float, float, float]:
    xyz=[]
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not xyz: raise ValueError(f"No coordinates found in {path}")
    return tuple(sum(v[i] for v in xyz) / len(xyz) for i in range(3))


def pdbqt_torsdof(path: Path) -> int:
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("TORSDOF"):
            return int(line.split()[-1])
    return 0


def write_gpf(receptor: Path, ligand_types: list[str], grid: GridBox, destination: Path) -> None:
    ligand_types=validate_ad4_types(ligand_types)
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


def _dpf_prefix(receptor_stem: str, ligand: Path, ligand_types: list[str], params: DockingParameters) -> list[str]:
    center=pdbqt_center(ligand)
    lines = [
        "autodock_parameter_version 4.2", "outlev 1", "intelec", f"seed {params.random_seed} {params.random_seed + 1}",
        f"ligand_types {' '.join(ligand_types)}", f"fld {receptor_stem}.maps.fld",
    ]
    for t in ligand_types:
        lines.append(f"map {receptor_stem}.{t}.map")
    lines += [f"elecmap {receptor_stem}.e.map", f"desolvmap {receptor_stem}.d.map", f"move {ligand.name}",
              "about %.3f %.3f %.3f" % center, "torsdof %d" % pdbqt_torsdof(ligand),
              ]
    return lines


def write_dpf(receptor_stem: str, ligand: Path, ligand_types: list[str], params: DockingParameters, destination: Path) -> None:
    lines = _dpf_prefix(receptor_stem, ligand, ligand_types, params)
    lines += ["tran0 random", "quaternion0 random", "dihe0 random",
              "rmstol %.2f" % params.rmstol, "ga_pop_size %d" % params.ga_pop_size,
              "ga_num_evals %d" % params.ga_num_evals, "ga_num_generations %d" % params.ga_num_generations,
              "ga_elitism 1", "ga_mutation_rate 0.02", "ga_crossover_rate 0.8", "ga_window_size 10", "ga_cauchy_alpha 0.0", "ga_cauchy_beta 1.0", "set_ga",
              "sw_max_its 300", "sw_max_succ 4", "sw_max_fail 4", "sw_rho 1.0", "sw_lb_rho 0.01", "ls_search_freq 0.26", "set_psw1",
              "ga_run %d" % params.ga_runs, "analysis"]
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_native_start_dpf(receptor_stem: str, ligand: Path, ligand_types: list[str], params: DockingParameters, destination: Path, *, epdb: bool = False, local_runs: int | None = None) -> None:
    """Write a DPF whose absolute tran0 is the ligand's crystal center."""
    center = pdbqt_center(ligand)
    lines = _dpf_prefix(receptor_stem, ligand, ligand_types, params)
    lines += ["tran0 %.3f %.3f %.3f" % center, "quaternion0 0.0 0.0 0.0 1.0"]
    torsdof = pdbqt_torsdof(ligand)
    if torsdof:
        lines.append("dihe0 " + " ".join(["0.0"] * torsdof))
    lines += ["rmstol %.2f" % params.rmstol, "sw_max_its 300", "sw_max_succ 4", "sw_max_fail 4", "sw_rho 1.0", "sw_lb_rho 0.01", "ls_search_freq 1.0", "set_psw1"]
    if epdb:
        lines += ["epdb"]
    elif local_runs is not None:
        lines += [f"do_local_only {local_runs}", "analysis"]
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_best_energy(dlg: Path) -> float | None:
    text = dlg.read_text(encoding="utf-8", errors="replace")
    values = []
    for pat in (r"Estimated Free Energy of Binding\s*=\s*([+-]?\d+(?:\.\d+)?)", r"DOCKED: USER\s+Estimated Free Energy of Binding\s*=\s*([+-]?\d+(?:\.\d+)?)"):
        values.extend(float(v) for v in re.findall(pat, text))
    return min(values) if values else None
