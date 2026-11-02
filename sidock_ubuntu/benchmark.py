"""Deterministic helpers for benchmark classification and comparison."""
from __future__ import annotations

def sampling_result(minimum_rmsd: float, threshold: float = 2.0) -> str:
    return "SAMPLING_PASS" if minimum_rmsd <= threshold else "SAMPLING_FAIL"

def ranking_status(lowest_energy_rmsd: float, minimum_rmsd_rank: int | None, has_native_like: bool, threshold: float = 2.0) -> str:
    if not has_native_like:
        return "NO_NATIVE_LIKE_POSE"
    if lowest_energy_rmsd <= threshold:
        return "NATIVE_LIKE_TOP1"
    if minimum_rmsd_rank is not None and minimum_rmsd_rank <= 3:
        return "NATIVE_LIKE_TOP3"
    return "NATIVE_LIKE_SAMPLED_NOT_TOP3"

def interaction_agreement(sidock: set[tuple[str, str]], reference: set[tuple[str, str]]) -> dict[str, int | float]:
    union = sidock | reference
    agree = len(sidock & reference)
    return {"true_agreement": agree, "sidock_only": len(sidock-reference), "reference_only": len(reference-sidock), "union": len(union), "percent_agreement": 100.0*agree/len(union) if union else 100.0}

def normalized_interaction_key(interaction_type: str, chain: str, residue_name: str, residue_number: str, protein_atom: str = "", ligand_atom: str = "") -> tuple[str, str, str, str, str, str]:
    """Canonical atom-level key; empty atom fields support residue-level comparison."""
    aliases = {"hydrogen_bond": "hydrogen_bond", "hbond": "hydrogen_bond", "hydrophobic": "hydrophobic", "pi-pi": "pi-pi", "pi_stack": "pi-pi", "pi-cation": "pi-cation", "salt_bridge": "salt_bridge", "halogen_bond": "halogen_bond", "metal_coordination": "metal_coordination"}
    return (aliases.get(interaction_type, interaction_type), chain or "", residue_name or "", str(residue_number or ""), protein_atom or "", ligand_atom or "")

def mismatch_classification(sidock_present: bool, plip_present: bool, *, same_residue: bool = False, protonation_difference: bool = False) -> str:
    if sidock_present == plip_present:
        return "AGREEMENT"
    if protonation_difference:
        return "PROTONATION_DIFFERENCE"
    if same_residue:
        return "ATOM_TYPING_DIFFERENCE"
    return "UNRESOLVED"
