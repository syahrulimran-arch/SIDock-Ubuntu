# SIDock Ubuntu

Headless Ubuntu automation derived from SIDock AD4 0.7.6. It is designed for Codex to orchestrate a reproducible docking workflow while chemistry and docking remain deterministic.

## Current v0.1 workflow

`protein.pdb + compounds.sdf` -> protein/ligand intake -> co-crystal ligand grid -> RDKit 3D generation/minimization -> Meeko PDBQT -> AutoGrid4 -> AutoDock4 -> docking-energy analysis -> IC50/docking statistics -> evidence bundle for Codex scientific writing.

The workflow deliberately stops when a suitable bound reference ligand cannot be detected. Blind AI-selected pockets are not enabled in v0.1.

**Current status:** this is the first Ubuntu/headless scaffold and must be integration-tested on the target Ubuntu workstation before scientific use. Full SIDock 0.7.6 pose/contact parsing and publication binding-pose figures are the next migration step.

## Ubuntu setup

```bash
git clone https://github.com/syahrulimran-arch/SIDock-Ubuntu.git
cd SIDock-Ubuntu
chmod +x install.sh check_install.sh
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

## Scientific safeguards

- Original input files are not silently overwritten.
- Ligand preparation uses RDKit ETKDG and MMFF94 where parameters are available, with UFF fallback.
- A co-crystallized ligand defines the grid automatically.
- Docking values are parsed from generated AutoDock files, not invented by the AI.
- Codex receives an `evidence.json` bundle and strict writing instructions.
- Correlation claims are calculated from matched experimental IC50 and docking outputs.
- Failed preparation/docking stages stop or are logged rather than being silently converted into results.

## Origin

The workflow is derived from SIDock AD4 0.7.6 and is being migrated into a GUI-free Ubuntu backend for Codex-controlled automation.
