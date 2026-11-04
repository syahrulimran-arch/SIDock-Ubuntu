# SIDock reviewed five-step docking workflow (v0.6-dev)

This is the default route for new projects. No docking begins until each
scientific decision has been explicitly reviewed. The legacy fully automated
route remains available for old projects but is not required by this workflow.

## Launch

On Ubuntu, from the cloned repository with its virtual environment active:

    git pull --ff-only origin v0.6-dev
    python -m py_compile sidock_ubuntu/gui.py sidock_ubuntu/guided_ui.py sidock_ubuntu/guided_workflow.py sidock_ubuntu/manual_docking.py sidock_ubuntu/stages.py
    pytest -q tests/test_guided_workflow_v06.py tests/test_manual_docking_v06.py tests/test_custom_ad4_parameters_v06.py
    sidock-gui

Open http://localhost:8000/projects/new, or for an existing project:
http://localhost:8000/projects/PROJECT_NAME/guided/1.

## 1. Create project

Enter only the project name. No protein or SDF is required at creation.
New guided projects have their own stored review state.

## 2. Review protein and binding site

Upload a PDB (unless already present). Inspect candidate bound ligands, ions,
modified protein residues, waters and solvents.

- Choose an explicit **crystal ligand** to define the site. SIDock extracts
  those ligand coordinates, computes a bounding-box center, adds configurable
  padding, and proposes even AutoGrid point counts. This is a *geometric grid
  proposal*, not AI pocket prediction or proof of a binding site.
- Or choose **no crystal ligand** and enter X, Y and Z center coordinates
  yourself. Grid point counts and spacing can be adjusted in Step 4.
- Select which hetero-residues and solvents to remove. Water removal is a
  separate explicit checkbox. Modified protein residues are protected by
  default; ions are retained by default. The cleaned PDB and extracted
  reference are downloadable for review.
- Prepare PDBQT with Meeko using an **independent** action, or import a
  pre-prepared scientifically reviewed PDBQT. Failed receptor preparation
  stops only this step, with a log. It does not launch AutoGrid/AutoDock.
- Approve Step 2 only after reviewing the actual receptor PDBQT and
  retained components.

For metal-containing proteins, Meeko or AutoDock's default force field may
not cover the metal centre. The GUI will not invent nickel parameters or
silently delete critical atoms. Deliberately choosing a metal-free model
requires an explicit exploratory acknowledgement.

## 3. Define, minimize and approve compounds

Upload a single SDF library; compound IDs must be unique and present.
Click **Minimize and prepare compounds** to create minimized SDFs/PDBQTs.
Alternatively, copy existing project-prepared PDBQTs into the independent
manual work area without reminimization, where filenames match exactly.
Inspect per-compound files and counts before approving Step 3. This
workflow uses the existing RDKit MMFF/UFF minimization pipeline and Meeko
PDBQT preparation; success is not a guarantee that every conformer
represents the biological bound state.

## 4. Docking setup and grid validation

Manually set/edit:

- Grid center X/Y/Z; even grid-point counts along X/Y/Z (20–126);
  spacing (normally 0.375 Å).
- Genetic algorithm runs per compound; evaluations, population size,
  generations; number of simultaneous compound jobs.
- Optional **complete** custom AutoDock4 parameter library. This file is
  copied into both AutoGrid's GPF and each AutoDock DPF; the protocol
  fingerprints the exact content. A metal atom type needs documented
  parameters, not a guessed placeholder.

Click **Generate AutoGrid4 maps**. Review the written GPF and log, then
explicitly approve Step 4. No compounds are docked in this step.

## 5. Docking run

Step 5 refuses to run if the approval hashes for Steps 2–4 no longer
match the currently selected protein, ligand preparation or grid setup.

Completed compounds are *not* rerun: SIDock accepts an existing DLG only
when it contains a recognized energy and the AutoDock 'Successful
Completion' marker. Incomplete DLG files are backed up before retrying.
AutoDock jobs can run in parallel. Outputs include each DLG and
a manually scoped results.csv.

## Data isolation and limitations

Guided artifacts are stored under
\`PROJECT/manual_docking/\` including separate \`runs/RUN_ID/\` paths.
Changing a grid or receptor makes a different run ID; prior DLG files remain
untouched. The legacy \`project_state.json\` and old docking results are not
rewritten by guided execution.

A crystal-ligand-derived box is not itself native redocking validation.
No redocking RMSD should be claimed from this manual workflow without a
separate, documented redocking analysis. Metal-aware atom types, charges,
coordination and custom parameter file selection require scientific review.
This workflow provides explicit preparation choices; it does not
automatically validate a metalloprotein force field.

The default GUI listens on localhost only. Use SSH port forwarding for
remote access unless network authentication has been configured.
