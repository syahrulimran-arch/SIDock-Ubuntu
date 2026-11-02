# SIDock Ubuntu production rules

These rules apply to controlled real-study use. `READY_FOR_CONTROLLED_REAL_STUDY` does not mean publication-ready or validated for all targets.

1. Original input PDB/SDF files are immutable.
2. Native-ligand redocking must pass at <=2.0 Å before study compounds are docked.
3. A target that fails native redocking blocks downstream study docking.
4. Target-specific parameter tuning to force validation is prohibited.
5. Raw DLG files must be retained.
6. SHA-256 hashes must be recorded for inputs and major derived artifacts.
7. CCD chemistry is authoritative for crystallographic reference ligands.
8. Study SDF chemistry is authoritative for study compounds.
9. Ambiguous study-ligand chemistry must be flagged, not silently repaired.
10. Fixed-frame heavy-atom RMSD must not use ligand superposition.
11. Corrected DLG traversal-to-chemical mapping is mandatory.
12. Lowest docking energy is not synonymous with the experimentally correct pose.
13. Cluster information must be retained.
14. PLIP is the primary interaction annotation engine.
15. Protonation is performed once and recorded.
16. PLIP runs use `--nohydro` on the preserved protonated complex.
17. SIDock internal interaction analysis is experimental/secondary unless explicitly promoted.
18. Docking energy is not experimental binding affinity.
19. IC50/pIC50 correlation is reported only from calculated statistics.
20. No residue interaction may appear in narrative text unless present in evidence files.
21. Failed ligands may not silently disappear from output.
22. Preparation, tautomer and protonation warnings remain visible.
23. Statistical calculations preserve sample count.
24. Results and Discussion generation uses `evidence.json` only.
25. AI-generated scientific statements must be traceable to evidence.
26. Human scientific review is mandatory before publication use.
27. Benchmark failures, including 1M17, remain in validation provenance.
28. This status does not mean publication-ready.
