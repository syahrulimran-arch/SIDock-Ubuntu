# Guided docking analysis

Step 5 now connects the manually reviewed/Open Babel docking route to saved-output
analysis. A successful `dock` worker automatically analyses that run. Existing
completed runs can instead use **Analyse existing results — no redocking**.
The **Other saved runs** links select historical runs independently of current
receptor/grid approvals; analysis never starts a new docking calculation.

## Outputs

Each analysis attempt writes a separate revision below
`manual_docking/runs/<run_id>/analysis/revisions/` and records the current revision
in that run's `analysis_result.json`. Step 5 links to:

- an HTML analysis report with figures and a draft Results and Discussion;
- the editable Markdown draft;
- score, cluster and PLIP interaction counts in CSV;
- residue identities, atom indices and recorded distances in per-compound PLIP CSVs;
- poses, complexes, PLIP XML and command logs;
- evidence/provenance JSON and a ZIP containing the analysis plus its original
  protocol, receptor, prepared ligands, DPF/DLG files and linked SDF when available.

Only completed DLGs with finite, energy-ranked poses enter the ranking. Their
atom identities must match the saved prepared ligand. The receptor used for
analysis is the run's frozen grid receptor, checked against `protocol.json`.
DPF settings are read from that saved run, not the current GUI settings.

PLIP is applied to the selected pose. Its output is filtered to the explicitly
materialized docked ligand `LIG:Z:999`, excluding unrelated receptor components.
A missing or failed PLIP result remains unavailable, with blank CSV counts and
JSON nulls. A successful report of no contacts yields zero. Both current
`report.xml` and older `*_report.xml` filenames are supported.

## Figures and experimental metadata

Figures include an energy ranking, a PLIP residue-contact matrix when contacts
exist, a ligand-centred **2D protein–ligand interaction diagram for each compound
with completed PLIP analysis**, linked input structures, selected-pose coordinate plots and an activity
scatter when at least two exact positive IC50 pairs exist. Ranking/contact
figures display up to 40 compounds and 30 residue columns; detailed CSV/JSON
retain every compound/contact. Structure/pose views cover the top five compounds.
Coordinate plots do not invent molecular bonds or assign contacts from proximity.
The 2D diagrams show recorded residue/chain names, interaction types and PLIP
distances. Ligand atom/group endpoints use PLIP's recorded SMILES-to-PDB mapping
only when its chemical graph matches the supplied structure. Otherwise lines
are explicitly labelled as ligand-level contacts. Dense contact sets are paginated
without dropping records. Completed PLIP reports with zero contacts produce an
explicit no-contact diagram; missing PLIP is never depicted as zero contacts.
Modern PLIP `pi_stacks` records and atom-group indices are retained. Diagrams
and their JSON source records are included in the HTML report and ZIP; Step 5
also provides individual PNG download links. Existing analyses must be rerun
after updating to add these diagrams; docking need not be repeated.
Metal positions are shown from the saved receptor. PyMOL rendering scripts are
always exported for these top poses; molecular PNGs are also rendered when the
`pymol` executable is available. All displayed figures identify the docking as
exploratory/unvalidated.

Experimental comparisons require an SDF linked by the preparation manifest and
prepared-ligand hashes. Numeric values follow the micromolar SDF schema unless
inline units, unit-specific property names or `IC50_unit(s)` metadata specify
otherwise. Recognized nM/mM/M values are converted to µM. Qualified values such
as `>10`, uncertainties and unsupported units are excluded from numerical
correlations. No IC50 is inferred from a docking score.

## Recovery and scientific status

A missing original SDF does not force redocking; experimental comparisons and
input-structure figures are then unavailable. Missing PLIP leaves a partial
analysis with explicit errors while preserving usable poses and scores. Repair
the analysis dependency and click the analysis button again. Original docking
inputs/logs remain unchanged. Failed attempts and earlier report revisions are
retained. A complete report is reused only if its recorded input/output hashes
still match. Edited, added or missing DLGs and modified report artifacts make
the saved report stale.

Successful analysis does not establish scientific docking validation. Guided
reports, summary result CSVs, evidence and plots retain `scientific_approval=false` and the saved
model status, with native redocking marked not performed. The small technical
smoke test is never represented as native RMSD validation.

CLI equivalent for an existing saved run:

```bash
python -m sidock_ubuntu.manual_docking PROJECT analyze --run-id RUN_ID
```

After pulling code changes, restart the existing GUI process with its original
workspace settings. Reloading the browser alone does not reload Python code.
