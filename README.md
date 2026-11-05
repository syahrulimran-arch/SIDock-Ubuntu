# SIDock Ubuntu

SIDock Ubuntu is an automated AutoDock4 molecular-docking workflow for Ubuntu. It performs target validation, CCD-aware reference chemistry reconstruction, receptor and ligand preparation, AutoGrid4/AutoDock4 docking, pose extraction, clustering, PLIP interaction analysis, figures, QC, and provenance tracking.

Receptor and ligand PDBQT preparation uses Open Babel (`obabel`). Meeko is no longer used or required.

The guided workflow now runs post-docking pose/cluster analysis, PLIP, figures
and an evidence-bound draft write-up. Step 5 also analyses existing saved runs
without redocking. See [guided analysis](docs/GUIDED_ANALYSIS.md) for outputs,
dependency recovery and the unvalidated-model status.

## Quick start

```bash
./install.sh
sidock doctor
sidock run protein.pdb compounds.sdf --project my_project
sidock status my_project
```

Native-ligand redocking is a hard gate: a target must achieve fixed-frame heavy-atom RMSD ≤2.0 Å before study-compound docking proceeds. The original PDB/SDF inputs and raw DLG files are retained.

Useful commands:

```bash
sidock config show
sidock config validate
sidock resume my_project
sidock inspect my_project --compound quercetin
sidock compliance my_project
```

AutoDock4 is CPU-based. StageRunner execution is serial in v0.4; `--jobs` is
accepted for compatibility and values above 1 emit a warning before continuing
with one job. Parallel StageRunner execution is deferred to v0.5.

## Reference ligands and safety gates

SIDock scans non-polymer HETATM groups, excludes common solvent/buffer/ion records, writes `validation/reference_ligand_candidates.csv`, and records its deterministic selection. Use `--reference-ligand RES:CHAIN:NUMBER` when an explicit instance is required. Interactive ambiguity handling is supported by the candidate list; noninteractive runs stop with `E101` when selection is not sufficiently supported.

Before receptor preparation, canonical protein residues are checked for missing heavy atoms. Incomplete residues within 5 Å of the reference ligand stop with `E202`. Residues 5–8 Å require `--allow-incomplete-peripheral`; preparation preserves existing heavy atoms and never retries by deleting unmatched residues.

PLIP is required for primary interaction annotations. The internal contact parser is retained only as a secondary cross-check.

Headless Ubuntu automation derived from SIDock AD4 0.7.6. It is designed for Codex to orchestrate a reproducible docking workflow while chemistry and docking remain deterministic.

## v0.4 production workflow

`protein.pdb + compounds.sdf` -> persisted StageRunner intake/reference/receptor/grid/validation -> gated compound stages -> disk-rehydrated finalization.

The workflow deliberately stops when a suitable bound reference ligand cannot be detected or native-ligand redocking validation fails. Blind AI-selected pockets are not enabled. v0.2 also records native-redocking RMSD, AutoDock clusters, and geometry-derived contacts.

**Release-candidate status:** v0.4 orchestration, restart/recovery, provenance, finalization, and local CCD/mmCIF reference selection are verified by the repository test suite. Remote reference enrichment and StageRunner parallel execution are deferred.

## v0.5 GUI foundation

The development branch provides a local FastAPI thin client (`sidock-gui`) for
creating/reopening projects, viewing persisted reference candidates, compound
results, COMPLETE/PARTIAL/BLOCKED project summaries, and downloading canonical
`results.csv`, `evidence.json`, and `report.html` artifacts. Result pages are
read-only and use persisted backend state; they do not run scientific code.

GUI execution is owned by a separate local worker process. Operational job
records and history are stored under `provenance/`; project scientific state
remains authoritative. A worker lost during GUI restart or host reboot is
marked `INTERRUPTED` and requires an explicit Resume. Cancellation and systemd
supervision are deferred; execution remains serial within each project.

### Project execution safety

Scientific writers share one backend-neutral Linux advisory lock at
`provenance/project_execution.lock`. CLI `run`/`resume` and GUI workers acquire
this same lock, so a project cannot be written concurrently by two interfaces.
The lock file's diagnostic JSON is not authoritative; `fcntl` ownership is.
Read-only status, inspection, results, and configuration views do not take the
exclusive lock. GUI job records (`gui_job.json` and
`gui_job_history.jsonl`) remain optional operational metadata and are not used
by the CLI. A lock conflict is an operational error and does not create a
scientific stage failure or execution-trace record.

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
sidock-ubuntu run protein.pdb compounds.sdf --project my_study --ga-runs 20 --jobs 1
```

SDF properties recognized for activity include `IC50`, `IC50_uM`, `IC50_UM`, `IC50 (uM)`, and `IC50 µM`.

## Scientific safeguards

- Original input files are not silently overwritten.
- Ligand preparation uses RDKit ETKDG and MMFF94 where parameters are available, with UFF fallback.
- A co-crystallized ligand defines the grid automatically.
- Docking values are parsed from generated AutoDock files, not invented by the AI.
- Codex receives an `evidence.json` bundle and strict writing instructions.
- Correlation claims are calculated only when experimental IC50 metadata is present; docking scores are predictions, not experimental affinity.
- Failed preparation/docking stages stop or are logged rather than being silently converted into results.

## Origin

The workflow is derived from SIDock AD4 0.7.6 and is being migrated into a GUI-free Ubuntu backend for Codex-controlled automation.
