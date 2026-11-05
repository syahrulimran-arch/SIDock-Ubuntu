# Guided GUI: hash-bound technical exploratory docking

**Purpose:** Reuse an *existing* Open Babel receptor PDBQT and an actual passed AutoGrid4/AutoDock4 smoke-test audit. This is not scientific or native-ligand validation of a metal-dependent binding model. No source PDB, reference ligand, previous receptor, SDF, completed DLG or prior run is overwritten.

## Existing-project steps

1. Open the SAME existing SIDock project that contains:
   - \`manual_docking/inputs/receptor_openbabel_exploratory.pdbqt\`
   - \`manual_docking/inputs/receptor_openbabel_exploratory_audit.json\`
   - earlier Open Babel structural, charged and hydrogen-typed candidates and audits
   - \`manual_docking/custom_parameters.dat\`
   - \`manual_docking/exploratory_smoke/exploratory_smoke_audit.json\` plus its real GLG/DLG/maps
2. Go to **Step 2**. The app checks the original protein/cleaned-PDB hashes, audit chain, unchanged receptor and AD4 parameter hashes, and the actual completed smoke logs/maps.
3. Tick the acknowledgement under **Separate exploratory receptor approval** and click **Approve for exploratory docking and continue**. The accepted state is \`EXPLORATORY_UNVALIDATED_NI_MODEL\`, \`scientific_approval=false\`. It does **not** copy the exploratory receptor onto \`manual_docking/receptor.pdbqt\`.
4. **Step 3** should show the existing prepared ligand set. If it says \`Verified preparations: Yes\`, click **Approve and continue**. No re-minimization is required. Files belong to individual projects: prepared ligands from another SIDock project are not silently imported.
5. In **Step 4**, review center, box coverage, spacing, GA settings and parameter library. The grid must cover the intended pocket and all tested ligand sizes; the previous 20×20×20 smoke-test box is NOT a production grid. Click **Generate AutoGrid4 maps** only when ready, inspect the GPF and then approve Step 4.
6. **Step 5** remains locked until Steps 2–4 are approved. Starting/resuming docking is a separate explicit click; no full-library docking begins by approving earlier stages.

## Integrity/safety

- Approvals expire when an upstream source, corrected receptor, parameter file, ligand manifest, smoke audit, GLG/DLG/map evidence, or relevant settings change.
- The active receptor for exploratory runs is taken from the hash-approved Open Babel candidate, **not** the standard scientific receptor file; it is staged into an independently hashed run directory.
- Grid, docking protocol, results and CSV include exploratory/unvalidated model flags. Energies are technical exploratory scores, not validated physical binding free energies.
- Generalized across structures by parsed artifacts and hashes. No 4UBP or protein-specific residue identifiers are used to select this workflow.
- The Open Babel conversion itself and any Ni-specific atom parameters are separate chemical choices; this workflow does not validate coordination chemistry or partial charges.
- The native scientific approval route remains independent and unchanged.

## Local verification

\`\`\`bash
cd ~/SIDock-Ubuntu
source .venv/bin/activate
python -m compileall -q sidock_ubuntu
pytest -q tests/test_exploratory_guided_approval_v06.py tests/test_openbabel_receptor_v06.py tests/test_guided_workflow_v06.py tests/test_manual_docking_v06.py
\`\`\`

Inspect the actual local project \`sidock_projects/taha-ali-almahroom5\`. Only claim end-to-end GUI PASS after direct local execution. A GitHub commit or mocked regression tests are not a live AutoGrid/AutoDock verification.
