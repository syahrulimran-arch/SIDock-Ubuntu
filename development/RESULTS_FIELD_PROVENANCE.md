# Results field provenance

| Field | Source |
|---|---|
| `compound_id` | input/preparation |
| `status` | compound stage state |
| `lowest_energy` | AutoDock DLG/ranking |
| `cluster_count` | raw-DLG clustering |
| `largest_cluster_population` | raw-DLG clustering |
| `selected_pose` | pose extraction |
| `hbond_count` | normalized PLIP only |
| `hydrophobic_count` | normalized PLIP only |
| `pi_stack_count` | normalized PLIP only |
| `salt_bridge_count` | normalized PLIP only |
| `warnings` | preparation, validation, and provenance |

`contacts.detect_contacts` is secondary and must not populate primary evidence.
## Reference-selection provenance

- `selected_reference_ligand` -> local PDB/mmCIF instance selected by the CCD-style chemical graph and crystallographic-quality ranking.
- `chemical_group_id` -> canonical chemical graph hash (CCD atoms/bonds, charges, aromaticity, and reliable stereochemistry).
- `mapping_fraction`, `missing_heavy_atoms`, `extra_heavy_atoms` -> observed instance to local CCD mapping.
- `occupancy_summary`, `b_factor_summary`, `rscc`, `rsr` -> local mmCIF/PDB crystallographic fields; unavailable values remain null.
- `selection_confidence`, `selection_method`, `selection_reason` -> `validation/reference_selection.json` and candidate CSV.
