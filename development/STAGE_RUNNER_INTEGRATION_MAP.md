# StageRunner integration map

The production CLI enters `cli.main`, builds a production context, and routes
compound execution through `StageRunner` → `executors.py` → the validated
stage seams. `pipeline.run_all` remains a legacy sequential compatibility API.
The validated scientific functions are:

| Stage | Existing implementation | Inputs | Outputs |
|---|---|---|---|
| INPUT/REFERENCE | `pipeline.prepare_project`, `intake.hetero_groups`, `intake.choose_reference_ligand` | PDB, config | intake and reference files |
| RECEPTOR | `core.clean_receptor_pdb`, `core.prepare_receptor` | PDB, selected ligand | receptor PDB/PDBQT |
| VALIDATION | `validation.prepare_native_ligand`, `validation.run_native_redocking` | receptor, maps, reference ligand | validation JSON/poses |
| GRID | `core.box_from_reference`, `core.write_gpf`, `core.run_command` | receptor and grid | GPF, GLG, maps |
| PREPARATION | `core.split_and_minimize_sdf`, `core.prepare_ligand` | SDF | prepared SDF/PDBQT |
| DOCKING | `core.write_dpf`, `core.run_command` with AutoDock4 | PDBQT, maps, DPF | DPF/DLG |
| POSE_EXTRACTION | `stages.run_pose_extraction_stage` → `analysis.extract_docking_poses` | raw DLG | poses, complexes, mapping metadata |
| CLUSTERING | `stages.run_cluster_stage` → `analysis.analyze_autodock_clusters` | raw DLG | clusters, populations, representatives |
| PLIP | existing complex/export and interaction path in production pipeline | complexes | normalized interactions |
| FIGURES | `figures.pymol_binding_figures`, optional activity figure | analysis outputs | PNGs |
| FINALIZATION | pipeline evidence/results/dashboard adapters | all project outputs | evidence, results, compliance |

`sidock_ubuntu/executors.py` is the adapter seam. It deliberately delegates
to existing functions through production callbacks; it does not duplicate
chemistry, DLG parsing, or interaction algorithms.

Global stage seams persist restart-safe contracts through
`sidock_ubuntu.stage_outputs`; context rehydration loads those contracts and
verifies recorded artifact hashes before reuse.
