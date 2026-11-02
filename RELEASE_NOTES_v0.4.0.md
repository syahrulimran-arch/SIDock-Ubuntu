# SIDock Ubuntu v0.4.0

SIDock Ubuntu v0.4.0 is a restart-safe, provenance-preserving AutoDock4
workflow for Ubuntu.

Highlights:

- StageRunner production orchestration with persisted global and compound state.
- Native-redocking validation gate with fixed 2.0 Å RMSD acceptance.
- Local CCD/mmCIF reference-ligand identification and deterministic ranking.
- Receptor completeness checks and explicit low-confidence E101 handling.
- PLIP as the primary interaction engine; internal contacts remain secondary.
- Disk-only restart-safe finalization, partial-project and BLOCKED-project output.
- Atomic final outputs, artifact hashes, execution traces, and CLI/package support.

Limitations and deferrals:

- AutoDock4 remains CPU-based and external executables are installed separately.
- StageRunner compound execution is serial in v0.4; `--jobs >1` warns and runs serially.
- Remote reference enrichment is deferred.
- Docking scores are predicted energies, not experimental affinities.
- Experimental IC50 metadata is optional and is never inferred from docking.
- Web GUI and parallel StageRunner execution are deferred to v0.5.

Release verification: 155 pytest tests pass. No frozen scientific study was
rerun for this release freeze.
