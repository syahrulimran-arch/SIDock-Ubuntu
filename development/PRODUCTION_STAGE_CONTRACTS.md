# Production stage contracts

This document records the currently confirmed callable seams. Extraction is
only complete when each callable delegates to an operation already used by the
validated workflow.

| Stage | Callable/backend | Inputs | Outputs | Failure |
|---|---|---|---|---|
| Reference | `pipeline.prepare_project` + `intake` + `reference` | PDB, optional local mmCIF/CCD, reference config | chemically grouped/ranked ligand instance, candidate CSV, selection JSON, grid box | E101 |
| Receptor | `core.clean_receptor_pdb`, `core.prepare_receptor` | PDB/reference | receptor PDB/PDBQT | E201/E202 |
| Native validation | `validation.prepare_native_ligand`, `validation.run_native_redocking` | receptor, maps, reference | validation artifacts | E301/E302 |
| Grid | `core.box_from_reference`, `core.write_gpf`, `core.run_command` | receptor/reference | GPF, GLG, maps | E501 |
| Preparation | `core.split_and_minimize_sdf`, `core.prepare_ligand` | SDF | prepared SDF/PDBQT | E401 |
| Docking | `core.write_dpf`, `core.run_command` | ligand, maps, DPF | DPF/DLG | E502 |
| Pose extraction | `analysis.extract_docking_poses` via `stages.run_pose_extraction_stage` | DLG, receptor | poses, complexes, ranking | E602 |
| Clustering | `analysis.analyze_autodock_clusters` via `stages.run_cluster_stage` | raw DLG | cluster table/ranks | E601 |
| Figures | `figures.pymol_binding_figures` | analysis outputs | figures | E999 |
| Finalization | `dashboard.write_results`, `dashboard.write_dashboard`, pipeline evidence block | analysis/provenance | results/dashboard/evidence | E999 |

## PLIP status

## Reference selection contract

Reference selection is completed before receptor preparation and native
redocking. Candidate identity is based on a canonical CCD-style chemical graph
(element, charge, bond topology/order, aromaticity, and reliable stereo data),
with explicit observed-to-CCD mapping descriptors. Equivalent instances are
ranked by mapping completeness, alternate-location ambiguity, occupancy,
RSCC/RSR when deposited, and local binding-site completeness, followed by a
stable crystallographic identifier. Docking, RMSD, PLIP, and compound results
are not inputs to selection.

The canonical artifacts are `validation/reference_ligand_candidates.csv` and
`validation/reference_selection.json`. Automatic low-confidence selection
raises `E101`; an explicit component/chain/residue override is still checked
for existence and usable mapping. PDB-only operation remains supported and
missing mmCIF quality metrics are recorded as null.

Ranking policy version 1 is lexicographic: mapping completeness, mapping
ambiguity, altloc quality, occupancy, RSCC/RSR (higher RSCC and lower RSR),
binding-site completeness, peripheral completeness, contact-environment
completeness, then crystallographic identifier. Remote incompleteness is not a
criterion. This ranking is pre-docking and never reads AutoDock, RMSD, PLIP,
or compound results.

`stages.run_plip_stage` delegates to `plip_runner.run_plip_analysis`, preserving
the recovered Open Babel and PLIP commands. `contacts.detect_contacts` remains
secondary only.
# Persisted global-stage result contracts

Global stage outputs are stored as versioned JSON under `validation/` and are
rehydrated with `sidock_ubuntu.stage_outputs`. The canonical files are:

- `reference_stage.json`
- `receptor_stage.json`
- `grid_stage.json`
- `validation_stage.json`

Each records project-relative artifact paths, SHA-256 hashes, warnings, and
stage-specific metadata. Loaders return `VALID`, `STALE`, or `MISSING` through
the loaded result status fields; they do not rerun scientific work.
