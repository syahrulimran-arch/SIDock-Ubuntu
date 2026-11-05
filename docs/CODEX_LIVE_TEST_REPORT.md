# SIDock Ubuntu live test report

Date: 2026-11-04

## Scope and baseline

The protocol in `docs/CODEX_LIVE_UBUNTU_TESTING.md` was executed on the actual
Ubuntu host. Project `taha-ali-almahroom5` does not exist in this checkout;
project `taha-ali-almahroom2` is the available reviewed project and was used for
inspection. No full compound-library or 150-run docking experiment was started.

Host: `imran-S2600CWR`

OS/kernel: Ubuntu host, Linux `5.4.0-216-generic`, `x86_64`

Working directory: `/home/imran/SIDock-Ubuntu`

Git branch/SHA: `v0.6-dev` / `b106953dfa008cbe10c371002768b22b0f41d80a`

Baseline worktree was clean. It remains clean after testing. Existing project
DLG count before testing: `0`.

Installed versions and tools:

- Python 3.13.12 (`python3` and `.venv/bin/python`)
- RDKit 2026.09.1
- Meeko 0.8.0
- AutoGrid4 4.2.6 (`/home/imran/.local/bin/autogrid4`)
- AutoDock4 4.2.6 (`/home/imran/.local/bin/autodock4`)
- `.venv/bin/mk_prepare_receptor.py` and `.venv/bin/mk_prepare_ligand.py`

The reviewed original PDB contains `HAE:C:800`, `ACE:A:0`, modified `KCX:C:220`,
and nickel atoms at `NI:C:798` and `NI:C:799`. The reviewed cleaned PDB hash
matches the saved review metadata.

Relevant baseline hashes:

| File | SHA-256 |
|---|---|
| `sidock_projects/taha-ali-almahroom2/input/protein.pdb` | `e395fca736871c98fbc456290740d6e606af72ebfe2a417cc166f7f31b3089b6` |
| `sidock_projects/taha-ali-almahroom2/manual_docking/inputs/protein_clean.pdb` | `7c394e4147c46628a87126ecec631d4d96f47c633e3b38ec4519c3d081aed9b2` |
| `sidock_projects/taha-ali-almahroom2/manual_docking/inputs/protein_choice.json` | `a4cd16524317f7cadb50527cd7f8315c08e7bf804507e3f3b8da532f5ce6e2e8` |
| `sidock_projects/taha-ali-almahroom2/manual_docking/logs/receptor_meeko.log` | `40c35462a7d22e7a77a7c5982cef1e609bd14f830790c2f6033d0517e487891f` |

## Commands and results

| Test | Exact command/result | Status |
|---|---|---|
| Python compilation | `.venv/bin/python -m py_compile sidock_ubuntu/gui.py sidock_ubuntu/guided_ui.py sidock_ubuntu/guided_workflow.py sidock_ubuntu/manual_docking.py sidock_ubuntu/nickel_reconstruction.py sidock_ubuntu/core.py sidock_ubuntu/stages.py sidock_ubuntu/receptor_diagnostics.py` — exit 0 | PASS |
| Required pytest suite | `.venv/bin/python -m pytest -q tests/test_guided_workflow_v06.py tests/test_manual_docking_v06.py tests/test_nickel_reconstruction_v06.py tests/test_custom_ad4_parameters_v06.py tests/test_receptor_diagnostics_v06.py` — `31 passed`, exit 0 | PASS |
| Live Meeko, Ni-free temporary copy | `.venv/bin/mk_prepare_receptor.py --read_pdb /tmp/sidock-live3-ikweyt/no_ni.pdb --default_altloc A -o /tmp/sidock-live3-ikweyt/receptor2 -p` — exit 1 | FAIL/BLOCKED |
| Ni reconstruction | Not run: no complete locally available, user-approved Ni AutoDock4 parameter library and explicit charge model | BLOCKED |
| AutoGrid smoke | Not run: no valid receptor PDBQT; Ni model prerequisites are absent | BLOCKED |
| AutoDock smoke | Not run: no valid receptor/maps; AutoGrid prerequisite is absent | BLOCKED |
| Native HAE validation | Not run: no validated receptor/maps/AutoDock result | BLOCKED |

The temporary live-test files were created under `/tmp/sidock-live3-ikweyt/`.
The temporary `no_ni.pdb` contained zero Ni atom records and retained KCX. The
copied original cleaned PDB hash was unchanged:
`7c394e4147c46628a87126ecec631d4d96f47c633e3b38ec4519c3d081aed9b2`.

## Actual failures and diagnosis

The prior project failure is recorded in
`sidock_projects/taha-ali-almahroom2/manual_docking/job.json` and
`sidock_projects/taha-ali-almahroom2/manual_docking/logs/receptor_meeko.log`:
standard Meeko preparation exits 1 because nickel is unsupported by the
standard templates. The log also identifies KCX as an unknown residue.

The permitted Ni-free live test removes only Ni from a temporary copy. Meeko
then reaches the retained KCX residue and attempts to generate its chemical
template. The complete live stderr reports:

`Could not download CIF file for KCX ... <urlopen error [Errno -2] Name or service not known>`

This is an environment/data prerequisite failure. KCX is a modified lysine in
the catalytic environment; it was not removed, guessed, or bypassed. The
original PDB and reviewed cleaned PDB were not modified. `--allow_bad_res` was
not used.

At the time of the earlier run, no local Ni library had been found. In this
continuation, `/home/imran/Desktop/AD4.1_bound_Ni.dat` was found with SHA-256
`38348aca5745036747c8cf9f0c39cd1caf03a8dbd97b5d8a9d3a48cbf297f58f`.
It has 36 well-formed `atom_par` records and includes both `Ni` and `NI`, but
the file itself labels those records “Provisional non H-bonding Nickel” and
states that atom-type parameters do not set nickel partial charges. Its notes
cite a ResearchGate discussion rather than a urease-specific calibration. No
site-specific `C:798`/`C:799` charge model or independent validation provenance
was found, so it is not scientifically sufficient for reconstruction.

## Offline KCX continuation

Meeko 0.8.0 has no KCX entry in its installed
`meeko/data/residue_chem_templates.json`, no KCX file in the local cache, and
no offline CCD cache entry. The installed CLI supports `--add_templates`.

The RCSB CCD definition used for the offline template is the [RCSB KCX ligand
summary](https://rcsb.org/ligand/KCX), with the [RCSB file download
service](https://www2.rcsb.org/docs/programmatic-access/file-download-services)
documenting the CCD definition endpoint. It identifies KCX as “LYSINE
NZ-CARBOXYLIC ACID”, a neutral `L-PEPTIDE LINKING` component with formula
C7H14N2O4. The template records the observed polymer atom names `CX`, `OQ1`,
and `OQ2`, and uses `C-term` on the backbone carbonyl C and `N-term` on
backbone N. The standalone CCD terminal hydroxyl is represented consistently
with the polymer form in the reviewed PDB; no Ni coordination bond is added.

Code changes:

- Added `sidock_ubuntu/data/kcx_residue_template.json`, an offline CCD-derived
  Meeko residue template.
- Updated `sidock_ubuntu/core.py` so `prepare_receptor` supplies that template
  only when KCX is actually present. It does not delete residues, use
  `--allow_bad_res`, or change coordinates.
- Added `tests/test_kcx_template_v06.py` covering the carbamyl atom names,
  peptide link labels, Meeko template loading, and command construction.

The isolated actual Meeko command succeeded:

`.venv/bin/mk_prepare_receptor.py --read_pdb /tmp/sidock-kcx-isolated-z3Vynt/kcx_only.pdb --default_altloc A -o /tmp/sidock-kcx-isolated-z3Vynt/receptor -p --add_templates sidock_ubuntu/data/kcx_residue_template.json`

Exit code was 0. The output contained all 12 input KCX heavy atoms plus three
Meeko-added hydrogens, including `CX`, `OQ1`, and `OQ2`; its PDBQT atom records
contained only `KCX`, so no other residue could have been silently deleted in
that focused test.

The full reviewed temporary Ni-free receptor was then tested through the
updated SIDock route using `.venv/bin/python` and
`sidock_ubuntu.core.prepare_receptor` on
`/tmp/sidock-kcx-live-A7fgCL/protein_for_meeko_without_Ni.pdb`. The route
correctly passed the KCX template: the Meeko log contains no `C:220` mismatch
and no KCX CIF/network fetch. It still exits 1 because the reviewed 4UBP input
has independent existing mismatches at `A:0`, `B:111`, `B:119`, `C:326`, and
`C:403`. No PDBQT was written. The original copied-clean PDB remained hash
`7c394e4147c46628a87126ecec631d4d96f47c633e3b38ec4519c3d081aed9b2` before and
after; the temporary Ni-free copy had 6,056 atoms versus 6,058 in the source
and zero Ni records.

Regression and rerun results:

- `.venv/bin/python -m py_compile sidock_ubuntu/core.py`: exit 0.
- `.venv/bin/python -m pytest -q tests/test_kcx_template_v06.py tests/test_guided_workflow_v06.py tests/test_manual_docking_v06.py tests/test_nickel_reconstruction_v06.py tests/test_custom_ad4_parameters_v06.py tests/test_receptor_diagnostics_v06.py`: `33 passed`, exit 0.

## Code changes and remaining gate

The genuine compatibility issue was the absent offline KCX template and
Meeko’s attempted network retrieval. That issue is fixed and KCX preparation
is proven in isolation and matched in the full receptor attempt. The ACE
terminal-cap compatibility issue is also fixed. Full receptor preparation
remains blocked by the four unrelated, pre-existing zero-occupancy residue
mismatches listed below; they require separate structural review and were not
silently deleted.

## Final gate

| Gate | Result |
|---|---|
| Unit/code tests | PASS — 33 tests |
| KCX offline template / isolated Meeko | PASS |
| Full reviewed Ni-free Meeko preparation | BLOCKED by B:111, B:119, C:326, C:403; ACE A:0 compatibility fixed |
| Ni reconstruction | BLOCKED |
| AutoGrid smoke | BLOCKED |
| AutoDock smoke | BLOCKED |
| Native HAE validation | BLOCKED |

Scientific limitations remain: a Ni-free receptor is not a validated urease
metal-site model; geometric contacts would remain putative; no binding score,
IC50, interaction claim, or validation statistic was generated.

Existing project files were not modified. The source tree now contains the
offline KCX template, compatibility fix, regression test, and this updated
report. Only temporary test files under `/tmp` were created; prior project logs
and results were preserved.

Single next action: obtain a scientifically approved decision for the four
zero-occupancy residues; after full receptor preparation succeeds, re-audit
the Desktop Ni library and obtain documented site-specific Ni charges before
any reconstruction or docking.

## Continuation: remaining full-receptor mismatches

The original PDB contains `REMARK 480 ZERO OCCUPANCY ATOM` entries identifying
the relevant side-chain atoms as potentially unreliable:

- `GLU B:111`: `CG CD OE1 OE2`
- `GLU B:119`: `CD OE1 OE2`
- `LYS C:326`: `CE NZ`
- `LEU C:403`: `CG CD1 CD2`

The raw Meeko/RDKit residue graphs confirm the failure is connectivity, not an
absent standard template. For example, the GLU side-chain carboxyl fragment is
disconnected from the beta/gamma chain in the generated graph, while the
standard GLU template requires a connected `CB-CG-CD-(OE1/OE2)` topology. The
PDB records have no alternate-location identifier to select. Their coordinates
are retained in the input, but the deposition explicitly warns that their
location/properties may not be reliable.

The reviewed ACE cap was a separate compatibility issue. Its PDB has an
explicit `LINK C ACE A:0 - N MET A:1` at 2.20 Å, outside Meeko’s automatic bond
cutoff. The temporary preparation route now passes `--blunt_ends A:0=0` for
this reviewed terminal cap, preserving its atoms and coordinates. ACE and KCX
then no longer appear in the mismatch log.

Actual temporary tests:

1. Without compatibility options,
   `.venv/bin/mk_prepare_receptor.py ... --add_templates ...` exited 1 with all
   five old mismatches.
2. With `--forgive_extra_bonds`, it still exited 1 with all four zero-occupancy
   side-chain mismatches; this option does not repair missing intra-residue
   connectivity.
3. Through the updated SIDock route,
   `.venv/bin/python` calling `sidock_ubuntu.core.prepare_receptor` on
   `/tmp/sidock-mismatch-live-qa4XGl/no_ni.pdb` exited 1 with only
   `B:111`, `B:119`, `C:326`, and `C:403`. The log contains no `A:0` or
   `C:220` mismatch and no PDBQT was written.

No safe chemically correct repair is available for the four zero-occupancy
residues without either modeling/untrustworthy coordinate replacement,
inventing bonds, or deleting atoms/residues. Those actions would change or
silently reinterpret the crystal structure and were not performed. Full
receptor Meeko therefore remains BLOCKED rather than being forced through QC.

The Desktop parameter audit was read-only. No copy was installed into a
project, and no nickel charges were invented.

The compatibility change and regression test were rerun with:

`.venv/bin/python -m py_compile sidock_ubuntu/gui.py sidock_ubuntu/guided_ui.py sidock_ubuntu/guided_workflow.py sidock_ubuntu/manual_docking.py sidock_ubuntu/nickel_reconstruction.py sidock_ubuntu/core.py sidock_ubuntu/stages.py sidock_ubuntu/receptor_diagnostics.py`

Exit 0, followed by the required tests plus `tests/test_kcx_template_v06.py`:
`33 passed`, exit 0. `git diff --check` also exited 0.

## Continuation: four zero-occupancy residues and modeling gate

The unchanged original PDB was analyzed for residue identity, occupancy,
geometry, and proximity to the crystallographic HAE and Ni atoms. Distances
below are minimum heavy-atom distances; the `zero-occ` column restricts the
minimum to atoms explicitly deposited with occupancy 0.00.

| Residue | Zero-occupancy atoms | Internal geometry issue | HAE min / zero-occ (Å) | Ni min / zero-occ (Å) | Binding-site assessment |
|---|---|---|---:|---:|---|
| GLU B:111 | CG, CD, OE1, OE2 | CG-CD = 2.117 Å; RDKit graph disconnects CD/OE carboxyl fragment | 74.089 / 77.729 | 72.922 / 76.535 | Distant structural residue |
| GLU B:119 | CD, OE1, OE2 | CG-CD = 2.057 Å; RDKit graph disconnects CD/OE carboxyl fragment | 75.031 / 77.818 | 75.031 / 77.818 | Distant structural residue |
| LYS C:326 | CE, NZ | CD-CE = 1.982 Å; terminal side chain is zero occupancy | 16.677 / 18.222 | 19.771 / 21.707 | Outside direct pocket; distant structural residue |
| LEU C:403 | CG, CD1, CD2 | CB-CG = 2.093 Å; RDKit graph disconnects CG/CD branch | 24.709 / 26.779 | 24.709 / 26.922 | Distant structural residue |

The HAE atom coordinates used were the five deposited atoms in `HAE:C:800`;
the Ni coordinates used were `NI:C:798` at (30.981, 72.613, 72.887) and
`NI:C:799` at (27.566, 73.489, 72.684). The nearest HAE/Ni distances are far
beyond direct contact distances, including for the unresolved atoms themselves.
This supports treating all four as structural/preparation issues rather than
active-site residues. It does not make the incomplete side chains chemically
valid for a production receptor.

## Experimental-structure comparison and modeling decision

The RCSB search identified compatible experimental urease structures including:

- `2UBP`, native *Sporosarcina pasteurii* urease, X-ray, 2.00 Å, with KCX and
  two Ni ions.
- `1IE7`, phosphate-inhibited *S. pasteurii* urease, X-ray, 1.85 Å, with KCX
  and two Ni ions.
- `1UBP`, *S. pasteurii* urease, X-ray, 1.65 Å.
- `1FWE`, a related urease/HAE/Ni structure with KCX in the deposited model.

These entries are documented by the [RCSB 2UBP structure
page](https://www.rcsb.org/structure/2UBP), [RCSB 1IE7 structure
page](https://www.rcsb.org/structure/1IE7), [RCSB 1UBP structure
page](https://www.rcsb.org/structure/1UBP), and the [wwPDB 1FWE validation
report](https://files.rcsb.org/validation/view/1fwe_full_validation.pdf).
However, direct coordinate retrieval from this Ubuntu host failed with DNS
exit 6 for `https://files.rcsb.org/download/1IE7.pdb`, `1FWE.pdb`, `1E9Y.pdb`,
and `2UBP.pdb`. No compatible coordinate file was present in the local home
directory/cache.

Because the missing side-chain atoms are explicitly zero occupancy and no
locally verifiable homologous coordinates or protein rotamer/energy package is
available, a restrained local optimization would be a newly modeled structure,
not an experimentally supported repair. A model could not be defensibly
established in this session. No atoms were moved, added, deleted, or assigned
new bonds; no repaired PDB was created.

The original experimental files remain unchanged:

- `input/protein.pdb` SHA-256:
  `e395fca736871c98fbc456290740d6e606af72ebfe2a417cc166f7f31b3089b6`
- `manual_docking/inputs/protein_clean.pdb` SHA-256:
  `7c394e4147c46628a87126ecec631d4d96f47c633e3b38ec4519c3d081aed9b2`

Full-receptor Meeko therefore remains blocked at B111, B119, C326 and C403.
The already verified KCX template and ACE handling were not rerun or changed.
No docking or compound-library work was started.

Ni readiness remains a separate blocker. The Desktop
`AD4.1_bound_Ni.dat` file is syntactically complete and has Ni/NI coverage,
but its own comments call the values provisional and state that charges are not
provided. No site-specific charges for `C:798` and `C:799`, urease-specific
calibration, or validated coordination model was supplied. Nickel
reconstruction remains BLOCKED.

## Continuation: Desktop-reference repair attempt

The requested coordinate files were present locally and were not downloaded:

| File | SHA-256 | Use |
|---|---|---|
| `/home/imran/Desktop/2UBP.cif` | `daa6c30a9c6ddd203d6b8f90179087ad3ec30bdb6b77e0dade864df47441237d` | compatible comparison |
| `/home/imran/Desktop/1UBP.cif` | `8f79dd211bbfc4c1fae76e09b363bbffc51497fa05ed06993b2fc2180e7d0b81` | compatible comparison |
| `/home/imran/Desktop/1IE7.cif` | `1b233ca21af2edc8f0dba794cb2f0fe4c26d04c4ed851f06c3380a05570aa6cc` | source for B:111/B:119 side chains |
| `/home/imran/Desktop/1FWE.cif` | `c47c5253213377034890f2ae1b2ce148bc8adc6b6eab397cc09ac75db13a7451` | rejected as non-corresponding |

Command used for the comparison: `.venv/bin/python /tmp/align_refs.py`, using
Gemmi, global sequence alignment and CA Kabsch fitting. The 4UBP-to-reference
chain correspondences were exact for B and C in 2UBP, 1UBP and 1IE7 (B: 122
residues; C: 570 residues). CA RMSDs were:

| Reference | chain B | chain C | assessment |
|---|---:|---:|---|
| 2UBP | 0.087 Å | 0.115 Å | same sequence/numbering; target side chains still contain zero-occupancy atoms |
| 1UBP | 0.104 Å | 0.219 Å | same sequence/numbering; target side chains still contain zero-occupancy atoms |
| 1IE7 | 0.121 Å | 0.154 Å | same sequence/numbering; B:111 and B:119 are complete and occupied |
| 1FWE | not corresponding | not corresponding | different chain lengths/sequence register; rejected |

The locally supplied 1IE7 coordinates were therefore used only for the two
defensible replacements. A fresh temporary model was generated by
`/tmp/repair_b111_b119.py`; the original 4UBP PDB was not edited. The global
chain-B fit was 0.1209 Å CA RMSD. The repaired model retained every original
atom and replaced these seven zero-occupancy atoms with transformed 1IE7
coordinates, anchored on the retained 4UBP CB atom:

- GLU B:111 — `CG`, `CD`, `OE1`, `OE2`
- GLU B:119 — `CD`, `OE1`, `OE2`

This anchoring preserves the reliable experimental CB coordinates while
retaining the experimentally observed 1IE7 side-chain rotamer. The resulting
bond lengths were B:111 CB-CG 1.508 Å, CG-CD 1.521 Å, CD-OE1 1.252 Å and
CD-OE2 1.258 Å; B:119 CB-CG 1.530 Å, CG-CD 1.502 Å, CD-OE1 1.256 Å and
CD-OE2 1.266 Å. Peptide C-N distances at both repaired residues were 1.332–
1.346 Å. The nearest nonlocal heavy-atom distances were 2.697 Å for B:111 and
3.045 Å for B:119; no sub-2.0 Å nonlocal heavy-atom clash was found.

The supplied `/home/imran/Desktop/4ubpclean-template.pdb` was inspected but
not used as experimental evidence. Its heavy-atom coordinates at all four
targets are identical to 4UBP; it changes the zero occupancies to 1.00 and
adds hydrogens without provenance showing that the unresolved heavy atoms were
experimentally rebuilt. It cannot establish a defensible repair for C:326 or
C:403.

Actual temporary Meeko command path:

`.venv/bin/python /tmp/run_partial_repair.py`

That script used `clean_receptor_pdb` to remove only HAE:C:800 for receptor
preparation, removed the two crystallographic Ni records for this intermediate
test, and invoked the repository's real `prepare_receptor` route. The result
was an expected failure only for `C:326` and `C:403`; B:111 and B:119 no longer
appear in the Meeko log. The log also contains no A:0 or KCX/C:220 mismatch.
No PDBQT was generated because the two unresolved C-chain residues still fail
standard LYS/LEU template matching. No `--allow_bad_res`, deletion option or
extra-bond bypass was used.

The temporary repaired PDB contains 6,809 coordinate records before HAE/Ni
receptor cleaning and still contains 12 KCX and 3 ACE records. The cleaned
Meeko input contains 6,056 records after removing only HAE and the two Ni
records; it contains no unexplained protein-residue deletion. The original
4UBP hash remains
`e395fca736871c98fbc456290740d6e606af72ebfe2a417cc166f7f31b3089b6`.

No defensible C:326 or C:403 model was established. 2UBP, 1UBP and 1IE7 all
retain zero-occupancy or incomplete coordinates for those target side chains,
and 1FWE is not a valid chain/sequence counterpart. A rotamer/energy rebuild
would introduce newly modeled heavy atoms rather than recover experimental
coordinates; it was not silently performed. Full receptor Meeko and PDBQT
generation therefore remain blocked at C:326 and C:403.

The affected regression set was rerun after the temporary repair work:
`.venv/bin/python -m pytest -q tests/test_kcx_template_v06.py
tests/test_guided_workflow_v06.py tests/test_manual_docking_v06.py
tests/test_nickel_reconstruction_v06.py tests/test_custom_ad4_parameters_v06.py
tests/test_receptor_diagnostics_v06.py` — `33 passed` in 1.46 seconds. No
repository code defect was identified in this repair step, and no existing
project or docking result was modified.

Ni reconstruction remains independently blocked. `/home/imran/Desktop/AD4.1_bound_Ni.dat`
has complete syntax and Ni/NI atom-type coverage, but its comments identify
the values as provisional and state that atom types do not define partial
charges. It provides no validated urease-site charges for C:798/C:799, no
coordination-state calibration, and no scientifically justified charge model.
The file is not sufficient to authorize Ni reconstruction or docking.

## Step 2 preparation/approval workflow repair

The live `taha-ali-almahroom2` state reproduced the GUI message. Its reviewed
choice and cleaned PDB exist, but
`manual_docking/receptor.pdbqt` and `inputs/receptor_review.json` do not. The
last job is `FAILED`; its actual Meeko log reports unknown NI residues and
states that NI appears to be a linking fragment requiring a parameterized
template. The GUI had nevertheless rendered an enabled Step 2 approval form,
so the user reached the generic approval exception instead of seeing the
workflow state.

The general workflow was corrected as follows:

- Step 2 approval now uses a read-only `receptor_approval_readiness` check and
  renders a disabled approval button whenever the final receptor PDBQT is
  absent, stale, unaudited, or has an unresolved metal-count mismatch. The
  server-side approval checks remain authoritative.
- Meeko diagnostics now expose residue/component-level failures plus generic
  corrective options: inspect occupancy, connectivity, alternate locations and
  terminal chemistry; use experimentally supported or separately validated
  repair coordinates; rerun Meeko and verify atom/residue retention; and use
  reviewed offline templates or separately audited metal treatment for unknown
  components. No deletion or `--allow_bad_res` path was added.
- Successful `prepare_protein` execution still runs independently of approval,
  then installs the generated PDBQT and now records preparation method, source
  PDB hash, output PDBQT hash, log path/hash, metal-free status and parsed
  diagnostics in `inputs/receptor_review.json`. Failed preparation cannot
  register a final receptor.
- Existing explicit nickel safeguards remain: a receptor with unsupported Ni
  is not approved as a native metal-site receptor until metal treatment is
  complete. No Ni charges or parameters were invented.

Live read-only regression command:

`.venv/bin/python - <<'PY' ... receptor_approval_readiness(Path('sidock_projects/taha-ali-almahroom2')) ... step2(...) ... PY`

It reported `ready: False`, no final receptor PDBQT, two input Ni ions and zero
prepared Ni ions. The rendered page included “Why Meeko failed”, the NI
unsupported-component diagnosis, corrective options, and a disabled
“Approve and continue” button. No project file changed during this check.

Affected regression command:

`.venv/bin/python -m pytest -q tests/test_guided_workflow_v06.py
tests/test_receptor_diagnostics_v06.py tests/test_kcx_template_v06.py
tests/test_nickel_reconstruction_v06.py tests/test_manual_docking_v06.py
tests/test_custom_ad4_parameters_v06.py`

Result: `35 passed` in 1.53 seconds. Targeted modules also passed Python
compilation and `git diff --check`.

The current 4UBP receptor is not ready. The remaining scientific prerequisites
are defensible repair of C:326 and C:403, followed by successful full-receptor
Meeko preparation, and an independently justified Ni parameter/charge model.
