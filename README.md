# SIDock Ubuntu

Headless Ubuntu automation derived from SIDock AD4 0.7.6. It is designed for Codex to orchestrate a reproducible docking workflow while chemistry and docking remain deterministic.

## Current v0.1 workflow

`protein.pdb + compounds.sdf` -> protein/ligand intake -> co-crystal ligand grid -> ligand 3D generation/minimization -> Meeko PDBQT -> AutoGrid4 -> AutoDock4 -> pose/contact analysis -> IC50/docking statistics -> figures -> evidence bundle for Codex scientific writing.

The workflow deliberately stops when a suitable bound reference ligand cannot be detected. Blind AI-selected pockets are not enabled in v0.1.

## Ubuntu setup

```bash
git clone https://github.com/YOUR_USER/SIDock-Ubuntu.git
cd SIDock-Ubuntu
./install.sh
source .venv/bin/activate
./check_install.sh
```

AutoGrid4 and AutoDock4 must be installed and available as `autogrid4` and `autodock4`. PyMOL is optional but recommended for automatic binding-pose PNGs.

## Run

```bash
sidock-ubuntu run protein.pdb compounds.sdf --project my_study --ga-runs 20 --jobs 4
```

SDF properties recognized for activity include `IC50`, `IC50_uM`, `IC50_UM`, `IC50 (uM)`, and `IC50 µM`.

## Important scientific safeguards

- Original input files are not silently overwritten.
- Ligand preparation uses RDKit 3D embedding and MMFF94 by default and records geometry audit metadata.
- A co-crystallized ligand defines the grid automatically.
- Docking results and interactions are parsed from generated files, not inferred by the AI.
- Codex receives an `evidence.json` bundle and a strict writing prompt.
- Correlation claims are calculated from matched experimental IC50 and docking outputs.

## Origin

The core AutoDock4, preparation, and analysis functions were adapted from SIDock AD4 0.7.6 for a headless Ubuntu workflow.
