# SIDock: Codex live Ubuntu testing and repair

Run this protocol with **Codex CLI on the actual Ubuntu host** in the cloned SIDock-Ubuntu repository. This is a local environment and scientific pipeline test, not merely code review or GitHub CI.

## Mission

Implement and live-test a **general-purpose receptor preparation system** for the five-stage SIDock workflow. Read `docs/GENERAL_RECEPTOR_PREPARATION_SPEC.md` before coding. Use actual installed Python, Meeko, AutoGrid4, and AutoDock4 to diagnose failures in the user's selected project. The current 4UBP project is a regression fixture, **not a receptor type the application should assume**. Repair source code as necessary, retest the failed stage, and write a factual local report. Never claim that mocked unit tests constitute scientific validation.

## Mandatory protections

- Start by inspecting existing project taha-ali-almahroom5 (and project taha-ali-almahroom2 only if relevant). Check which project has the current manually reviewed protein and selected ligand. Do not create new projects to evade previous failures.
- Record git status, current branch and HEAD. Preserve dirty worktree, original source files, prior projects, generated maps and completed DLG files. Never git reset --hard, git clean, remove original data, or overwrite historical results.
- No PDB-ID-specific branches, no hardcoded chain/residue numbers, no assumption of exactly two Ni sites or any other fixed element count. Detect and classify structural features from each uploaded structure. Keep Ni, KCX and ACE support as optional chemistry-specific strategies rather than universal defaults.
- Work on copies of PDB, SDF and PDBQT inside a fresh mktemp -d temporary directory for live chemistry tests. Only change original project inputs after separate explicit approval.
- Never start the full 16-compound / 150-GA-run docking experiment as a test. One-ligand smoke tests only, and only when prerequisites are present.
- Never delete catalytic Ni or KCX to force a successful result. Never automatically guess nickel AutoDock parameters, nickel partial charges, or RMSD; never call an unvalidated fixed-centre metal model scientifically validated.
- Do not assume the conversation-side /mnt/data path is on Ubuntu. The user reports AD4.1_bound_Ni.dat is on their Desktop. First check "$HOME/Desktop/AD4.1_bound_Ni.dat" with test -f; if absent, use find "$HOME" -maxdepth 5 -type f -iname 'AD4.1_bound*Ni*.dat' -print to locate the actual host file (including case variations or redirected Desktop folders). Record the exact resolved path and SHA256. Never claim the file is missing without checking the Desktop first. Never create a fabricated replacement, and do not overwrite the original file.
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

## PRIORITY: Real Step 2 PDBQT conversion without Meeko (Open Babel)

The user explicitly requires a receptor PDBQT and does not want Meeko to block obtaining the file. SIDock now has a `receptor_openbabel` Step 2 action, implemented in `sidock_ubuntu/openbabel_receptor.py`. It generates `manual_docking/inputs/receptor_openbabel_candidate.pdbqt`, not an automatically approved final docking receptor.

1. Inspect the current requested project's actual `manual_docking/inputs/protein_clean.pdb` and `protein_choice.json` hashes (do not assume it is the Step 3 project). Confirm `obabel` is available with `command -v obabel`, and obtain permission before installing missing dependencies.
2. Test the new worker route directly on the user's current project: `python -m sidock_ubuntu.manual_docking <ABSOLUTE_PROJECT_PATH> receptor_openbabel` after reviewing whether its job manager already has an active process. Do not invoke Meeko, AutoGrid or docking.
3. Confirm Open Babel's rigid PDBQT options `-xr -xc -xn -xp -xh`. Check actual program exit status and log `manual_docking/logs/receptor_openbabel.log`, and confirm that a nonempty `.pdbqt` file exists.
4. Verify the audit's original-heavy-atom preservation, source/target SHA256, Ni/other metal retention and coordinate identity. If no atoms were lost, report conversion success even if chemical validation is still pending. If the audit reports loss, inspect precisely which residue/atom/coordinates were changed by Open Babel and fix a general serialization problem where feasible. Never silently delete these atoms.
5. Confirm the Step 2 GUI visibly shows the export button, link to the PDBQT, and link to the audit. Do **not** represent Gasteiger partial charges or raw rigid PDBQT conversion as a validated metalloprotein receptor; only docking approval requires subsequent chemistry review.
6. Run `pytest -q tests/test_openbabel_receptor_v06.py` and all affected guided workflow tests. Record exact PASS/FAIL/BLOCKED evidence in `docs/CODEX_LIVE_TEST_REPORT.md`. Preserve all original project data and do not rerun compound minimization or docking.

## Phase C: Live Meeko

Use a temporary copy of the reviewed receptor. Run actual installed Meeko through the Ni-free intermediate path. Verify that Ni is removed only from the temporary input and that original PDB and cleaned PDB remain untouched. Verify KCX carbamyl atoms CX/OQ1/OQ2 survive preparation if it succeeds. Capture actual exit code and stdout/stderr. Do NOT install a metal-free intermediate as a native production receptor.

Codex already reported a successful isolated offline KCX template test. Preserve that implementation and its regression tests. The subsequent full receptor test still identifies A:0, B:111, B:119, C:326, C:403 as unmatched; those are the immediate next blockers. Inspect their exact residue names, atom sets, adjacency, polymer connectivity and terminal/alternate-location annotations against the original PDB and Meeko's template matcher. Determine whether A:0 is an acetyl cap bonded to the protein rather than a free ligand, and inspect the other four residues independently. Fix template/representation problems narrowly on temporary copies, rerun only the Meeko intermediate stage and add targeted regression tests. Do not silently remove residues, use --allow_bad_res as a shortcut, or invent missing atoms.

## Separate Step 3 live verification: all ligands, even with receptor blocked

SIDock's SDF minimization and Meeko ligand-PDBQT conversion **must work before Step 2 is approved**, without pretending the metal-containing receptor is valid. Specifically inspect the existing project `sidock_projects/taha-ali-almahroom4` (the user's Step 3 screenshot: 16 compounds `AA-1` through `AA-16`). This is one regression project, not a hardcoded protein path or ligand-library special case.

1. Preserve its original `input/compounds.sdf`, all existing `ligands/`, `manual_docking/`, finished grid maps/DLG and project approvals. Record source fingerprint and current output counts/hashes before any preparation. Do not recreate the project.
2. Verify Step 3 now displays **Minimize and prepare compounds** above the compound table even without Step 2 approval. The POST handler must allow the `prepare` and `reuse_compounds` actions independent of protein approval.
3. With user-approved normal ligand preparation, test **real** installed RDKit/MMFF or UFF and Meeko on this exact 16-compound SDF, without launching receptor preparation, AutoGrid or AutoDock. Inspect individual Meeko logs and diagnose any actual compound errors. If an earlier manual ligand job is running, do not start a second.
4. Confirm exactly 16 distinct minimized SDF and 16 nonempty ligand PDBQT files with atoms and expected IDs, and a source-SDF-hash-matched ligand manifest. If interrupted, preserve completed outputs and test resume; never overwrite unexpected outputs.
5. Refresh Step 3 and confirm `Molecules: 16 | Prepared: 16 | Verified preparations: Yes`, or report the exact remaining blocker and file/log. Step 3 **approval** must remain pending until Step 2 receptor approval. No docking is authorized.
6. Test a second general-purpose, non-4UBP project to make sure this separation is truly independent of the protein, selected grid site and metal status. Run targeted automated tests and record local git commit/status in the report.
7. Record the difference between `taha-ali-almahroom2` (previous Step 2 tests) and `taha-ali-almahroom4` (Step 3 screenshot); do not report tests on one as successful testing of the other.

## Phase D: Nickel restoration and custom AutoDock parameters

Only if Meeko succeeds AND a complete locally available user-approved parameter file and individually justified Ni charge model exist:

- Read the user-provided Desktop parameter library at its **confirmed Ubuntu path** (usually "$HOME/Desktop/AD4.1_bound_Ni.dat"), verify its SHA256, both atom_par Ni/NI records, column counts and full receptor/ligand type coverage. The existence of the file removes the *file-not-found* blocker but does **not** certify its provisional community metal parameters as valid for urease. Audit the source and request a justified site-specific Ni charge model separately; do not fill missing charges by guesswork.
- Use the actual nickel_reconstruction.py route on temp/copied inputs; verify exact crystallographic Ni positions, partial charges, atom types, Ni-Ni distance, PDBQT fixed columns and no duplicate metal atoms.
- Inspect that identical parameter files are referenced from GPF and DPF and that the file follows to each run directory.
- Confirm approvals become stale if receptor, parameters or nickel audit change.
- A PDBQT with Ni is still an unvalidated metal-coordination model.

## Phase E: Limited AutoGrid/AutoDock smoke and validation gate

If the prerequisites are ready, perform one-ligand test with a tiny temporary AutoGrid/AutoDock setup, and inspect actual GPF, GLG, maps, DPF, DLG, exit codes and successful-completion marker. A smoke test is not a scientifically validated binding model. Native HAE redocking (appropriate heavy-atom RMSD and Ni-O geometry), original crystallographic ligand correctness and parameter provenance are separate required scientific validation steps.

If any prerequisite is missing, stop and explicitly mark the downstream test BLOCKED; never start the full compound library to compensate.

## Mandatory generalization acceptance tests

Before declaring SIDock generally ready, test at least:
- A standard non-metalloprotein lacking a crystal ligand and modified residues.
- A protein with a different ligand and chain/residue numbering from the current test case.
- A protein with nonstandard modified/covalent residues requiring an offline chemical component template.
- An ion/cofactor-containing protein whose ion is **not** Ni, to confirm the Ni-specific path is not launched or misapplied.
- A protein with zero-occupancy/missing side-chain atoms, to confirm the problem is reported and repair is reviewed, not silently forced.
- The current 4UBP test as an additional regression case only.

Document whether each scenario was unit tested, live tested or blocked. Do not claim general metal parameterization, automatic chemical template generation or reference-based structural repair is complete until it has been built and tested independently.

## Deliverables

Write docs/CODEX_LIVE_TEST_REPORT.md with date, host, git SHA, exact commands/exit codes, relevant file hashes and log paths. Make a PASS/FAIL/BLOCKED table for unit tests, live Meeko, Ni reconstruction, AutoGrid smoke, AutoDock smoke and native HAE validation. List code fixes made and tests rerun, scientific limitations still open, whether any project files were modified, and the single next action. Do not push, merge, or run all 16 compounds without approval.
