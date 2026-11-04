# General-purpose receptor preparation: design requirements

SIDock must be **PDB-agnostic**. No production decision may depend on a specific PDB ID, chain ID, residue number, or known example such as 4UBP. Use problematic structures as regression cases, not as hidden application configuration.

## Core contract

1. **Discover, do not assume**: inspect chains, sequences, HET groups, modified polymer residues, water, ions, caps, ligands, alternative conformations, missing atoms, occupancy, CONECT/LINK and polymer connectivity. Report source provenance.
2. **Reviewable decisions**: the user selects the reference ligand/site, removal of nonpolymer groups, treatment of waters, ion handling, and structural repair. Do not identify a crystal ligand solely from a residue code or treat a covalent cap as a free inhibitor.
3. **Offline-first chemical templates**: resolve residue templates from the installed toolkit/cache first, then local user-supplied CCD/template files; request internet retrieval only when allowed and available. Templates must match backbone linkage, bond order, stereochemistry and atom naming. A validated KCX template is one example, not an instruction to assume KCX in other proteins.
4. **Poorly determined atoms**: identify missing/zero-occupancy atoms and template mismatches per residue. Report the exact affected atoms, distances to the user-reviewed site, alternative experimental reference coverage, fit quality, and geometry. Never silently delete atoms or invent experimental coordinates. Any repairs are explicitly labeled *modeled*, isolated from the original, and independently approved.
5. **Metal/cofactor handling**: detect metal element and all coordination sites from the input. Support element-specific extensible preparation paths. Existing Ni support is an optional plugin/branch activated only when Ni is present, not a default for all proteins. Do not apply Ni parameters to Zn, Fe, Mg, Mn, Cu, Co or other species. Never assume two metal atoms or fixed chain/residue IDs. Site charges, atom types, coordination treatment and custom AD4 parameters require documented scientific provenance.
6. **Pipeline stages and integrity**: separate source input, reviewed clean structure, Meeko-compatible temporary intermediate, audited final receptor and validation. Use content hashes and provenance to invalidate downstream approvals when chemistry changes. Never automatically start large-scale docking from receptor preparation.
7. **Validation**: syntactic PDBQT and AutoGrid/AutoDock smoke checks are distinct from scientific validation. Native-ligand redocking is appropriate only when an experimental ligand exists; evaluate RMSD, binding-site coverage and metal-coordination geometry as relevant. Clearly report when validation is not possible.

## General fallback for unresolved non-active-site side chains

When experimental reference coordinates are unavailable, SIDock may offer a **modeled side-chain repair**, not a claim to reconstruct experimental coordinates. Apply this policy to *any* residue, chain or PDB ID based on detected missing/zero-occupancy atoms and structural context.

1. Preserve original crystallographic PDB and all reliable occupied coordinates. Record exactly which atoms have occupancy zero, are missing, or are already modeled.
2. Determine whether the affected side chain is part of the binding site, metal coordination, covalent linkage, interface or structural/catalytic interaction. If yes, stop for expert structural review; distance alone is not proof of irrelevance.
3. For distal side chains where modeling is acceptable, sample chemically valid, backbone-conditioned side-chain rotamers using an installed trustworthy offline modeling tool or rotamer library; do not select arbitrary coordinates or rename residue identities.
4. Score candidates by local steric overlap, bond lengths/angles, stereochemistry, hydrogen bonding and contact environment. Optionally apply restrained local optimization with **backbone and experimentally supported atoms fixed or tightly restrained**. Do not permit changes to ions or user-reviewed ligand/site geometry.
5. Retain multiple plausible models where the best rotamer is ambiguous. Record uncertainty and check whether differing rotamers affect pocket electrostatics or docking-site geometry. Never label modeled atoms as experimentally observed.
6. Write a **separate repaired receptor**, along with an atom-level modification map, modeled-atom flags, rotamer selection rationale, clash/geometry measurements, references, software versions and reproducibility information. Require explicit user approval before substituting it into Step 2.
7. Only after review: rerun real Meeko on the modeled copy and compare atom retention, modified residues, protonation and binding-site consistency. A Meeko PASS is a *format/preparation* check, not proof of experimental structure.
8. If no credible modeled configuration can be produced or uncertainty is relevant to the active site, report BLOCKED rather than deleting atoms or suppressing Meeko errors.

**Regression requirement:** The application must detect the same defects even if residue numbers, chain IDs and PDB IDs change. Tests may include 4UBP B:111/B:119/C:326/C:403 as one data set, but no remediation route may depend on their literal identifiers.

## Required implementation approach

- **Metadata-driven preparation**: choose behaviors based on parsed PDB/CCD metadata and explicitly reviewed component roles. No conditions such as `if pdb_id == "4UBP"`, no hardcoded chain/residue lists, and no filename-specific magic.
- **Extensible strategies**: make offline templates and metal handlers discoverable by residue/element identity and supported parameters. If a treatment is unavailable, return an actionable reviewable blocker rather than silently replacing it with an unrelated chemistry model.
- **Existing specializations**: keep current Ni, KCX and ACE handling as clearly labeled compatible strategies where their chemistry applies; do not remove safeguards, generalize an Ni routine to other metals by renaming it, or infer binding-site identity from those residue codes.
- **Repair policy**: use original experimental coordinates unchanged; write every repaired structure separately, with modeled atom IDs, reference structure IDs, alignment and local-validation metrics, plus human approval. Support arbitrary affected residues and chains.
- **Coverage**: test a standard protein without modifications or metals, a protein with a nonstandard covalent residue, an example with metal/cofactor(s), a structure with zero-occupancy side chains, a structure with an unusual terminal cap, and a structure lacking a native ligand. Use 4UBP only as one regression fixture.
- **No uncontrolled dependency retrieval**: prefer installed and user-provided local files; log failures and sources. Do not assume RCSB or any specific external host is reachable.

## Acceptance conditions

- A standard protein is not forced into a Ni or urease-specific preparation path.
- A PDB without HAE or ACE is processed without warnings mentioning them.
- A protein with different chain/residue numbering receives equivalent checks.
- Meeko failures expose the actual unmatched residues/atoms and maintain original files.
- A metal-containing model cannot silently bypass absent element-specific parameters or charges.
- Neither the GUI nor worker hardcodes 4UBP-specific binding-site decisions.
- Unit tests and real host integration tests demonstrate each applicable branch; claims of scientific validation are supported by native-ligand evidence where possible.

## Current status

The existing SIDock Ni restoration workflow is a **Ni-specific strategy**, not a generic metal treatment. The new receptor-repair workflow remains a design/implementation task until Codex's local report shows tested general coverage. Do not claim arbitrary metalloproteins or automatic side-chain modeling are implemented before they are.
