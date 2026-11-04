# SIDock: Codex live Ubuntu testing and repair

Run this protocol with **Codex CLI on the actual Ubuntu host** in the cloned SIDock-Ubuntu repository. This is a local environment and scientific pipeline test, not merely code review or GitHub CI.

## Mission

Run real installed Python, Meeko, AutoGrid4, and AutoDock4 to diagnose and correct failures in the reviewed five-stage SIDock workflow. Repair source code as necessary, retest the failed stage, and write a factual local report. Never claim that mocked unit tests constitute scientific validation.

## Mandatory protections

- Start by inspecting existing project taha-ali-almahroom5 (and project taha-ali-almahroom2 only if relevant). Check which project has the current manually reviewed protein and selected ligand. Do not create new projects to evade previous failures.
- Record git status, current branch and HEAD. Preserve dirty worktree, original source files, prior projects, generated maps and completed DLG files. Never git reset --hard, git clean, remove original data, or overwrite historical results.
- Work on copies of PDB, SDF and PDBQT inside a fresh mktemp -d temporary directory for live chemistry tests. Only change original project inputs after separate explicit approval.
- Never start the full 16-compound / 150-GA-run docking experiment as a test. One-ligand smoke tests only, and only when prerequisites are present.
- Never delete catalytic Ni or KCX to force a successful result. Never automatically guess nickel AutoDock parameters, nickel partial charges, or RMSD; never call an unvalidated fixed-centre metal model scientifically validated.
- Do not assume the conversation-side /mnt/data path is on the Ubuntu host. If AD4.1_bound_Ni.dat is missing locally, report it missing. Do not create a fabricated replacement.
- No sudo, unexpected external network access, unsafe permissions, pushing to GitHub or other changes outside the repository without human approval.

## Phase A: Host and project inspection

1. Capture date, hostname, OS, current working directory, git status, branch and SHA, Python/.venv versions, imports for RDKit and Meeko, and availability/path/version for mk_prepare_receptor.py, mk_prepare_ligand.py, autogrid4 and autodock4.
2. Inspect user project configuration and existing Step 2 artifacts, especially input/protein.pdb, manual_docking/inputs/protein_choice.json, protein_clean.pdb, manual_docking/job.json, logs/receptor_meeko.log, logs/receptor_meeko_Ni_intermediate.log, and input Ni reconstruction files if present.
3. Verify what the PDB actually contains: HAE:C:800 versus ACE:A:0, modified KCX C220, and two Ni C798/C799; do not infer this from folder names.
4. Record file hashes and existing DLG count before running any test.

## Phase B: Real local test suite

Run Python compilation of sidock_ubuntu/gui.py, guided_ui.py, guided_workflow.py, manual_docking.py, nickel_reconstruction.py, core.py, stages.py, and receptor_diagnostics.py.

Run pytest on:
- tests/test_guided_workflow_v06.py
- tests/test_manual_docking_v06.py
- tests/test_nickel_reconstruction_v06.py
- tests/test_custom_ad4_parameters_v06.py
- tests/test_receptor_diagnostics_v06.py

For any failure: inspect the complete traceback, reproduce, fix narrowly, add regression tests and rerun. These are unit/integration-code tests, not proof of actual Meeko/AutoDock scientific correctness.

## Phase C: Live Meeko

Use a temporary copy of the reviewed receptor. Run actual installed Meeko through the Ni-free intermediate path. Verify that Ni is removed only from the temporary input and that original PDB and cleaned PDB remain untouched. Verify KCX carbamyl atoms CX/OQ1/OQ2 survive preparation if it succeeds. Capture actual exit code and stdout/stderr. Do NOT install a metal-free intermediate as a native production receptor.

If Meeko reports A:0, B:111, B:119, C:326, C:403 or other residue mismatches, diagnose exact residue/bond/terminus/altloc issues from PDB and logs. Do not silently remove residues, use --allow_bad_res as a shortcut, or invent missing atoms.

## Phase D: Nickel restoration and custom AutoDock parameters

Only if Meeko succeeds AND a complete locally available user-approved parameter file and individually justified Ni charge model exist:

- Inspect the uploaded parameter library atom_par Ni/NI records for full atom-type coverage. Audit their source; do not treat the previous provisional community Ni values as validated for urease.
- Use the actual nickel_reconstruction.py route on temp/copied inputs; verify exact crystallographic Ni positions, partial charges, atom types, Ni-Ni distance, PDBQT fixed columns and no duplicate metal atoms.
- Inspect that identical parameter files are referenced from GPF and DPF and that the file follows to each run directory.
- Confirm approvals become stale if receptor, parameters or nickel audit change.
- A PDBQT with Ni is still an unvalidated metal-coordination model.

## Phase E: Limited AutoGrid/AutoDock smoke and validation gate

If the prerequisites are ready, perform one-ligand test with a tiny temporary AutoGrid/AutoDock setup, and inspect actual GPF, GLG, maps, DPF, DLG, exit codes and successful-completion marker. A smoke test is not a scientifically validated binding model. Native HAE redocking (appropriate heavy-atom RMSD and Ni-O geometry), original crystallographic ligand correctness and parameter provenance are separate required scientific validation steps.

If any prerequisite is missing, stop and explicitly mark the downstream test BLOCKED; never start the full compound library to compensate.

## Deliverables

Write docs/CODEX_LIVE_TEST_REPORT.md with date, host, git SHA, exact commands/exit codes, relevant file hashes and log paths. Make a PASS/FAIL/BLOCKED table for unit tests, live Meeko, Ni reconstruction, AutoGrid smoke, AutoDock smoke and native HAE validation. List code fixes made and tests rerun, scientific limitations still open, whether any project files were modified, and the single next action. Do not push, merge, or run all 16 compounds without approval.
