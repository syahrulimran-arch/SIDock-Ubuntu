# Codex instructions for SIDock Ubuntu

When the user supplies a receptor PDB and ligand SDF, use the deterministic SIDock Ubuntu pipeline for all chemistry and docking calculations. Do not manually invent or repair scientific results.

1. Activate `.venv` and run `./check_install.sh`.
2. Run `sidock-ubuntu run PROTEIN.pdb COMPOUNDS.sdf --project PROJECT_NAME`.
3. If the pipeline stops, diagnose the logged scientific/tooling failure. Do not bypass QC merely to obtain a result.
4. Read `PROJECT_NAME/analysis/evidence.json` and the CSV/JSON files it references.
5. Follow `PROJECT_NAME/analysis/CODEX_DISCUSSION_PROMPT.md` and write `PROJECT_NAME/analysis/Results_and_Discussion.md`.
6. Numerical claims must come from generated outputs. Do not invent interactions, residues, scores, IC50 values, statistics, or validation claims.
7. Treat geometric contacts as putative contacts, not experimental proof.
8. Preserve original inputs and logs.

Current v0.1 automatic grid rule: a co-crystallized reference ligand must be detectable in the receptor PDB. If no suitable ligand is found, stop rather than performing blind docking.
