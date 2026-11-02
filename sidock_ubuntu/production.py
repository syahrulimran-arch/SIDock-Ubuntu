"""Production-safety policy for controlled SIDock studies."""
from __future__ import annotations

PRIMARY_INTERACTION_ENGINE = "PLIP"
SIDOCK_INTERACTION_STATUS = "EXPERIMENTAL_SECONDARY"


def readiness_gate(*, benchmark_count: int, sampling_pass_count: int,
                   chemistry_ok: bool, rmsd_ok: bool, mapping_ok: bool,
                   ccd_ok: bool, plip_export_ok: bool, protonation_ok: bool,
                   provenance_ok: bool, validation_gate_ok: bool,
                   interaction_primary_ok: bool) -> dict[str, object]:
    checks = {
        "three_benchmarks": benchmark_count >= 3,
        "two_of_three_sampling_pass": sampling_pass_count >= 2,
        "chemistry_preparation": chemistry_ok,
        "fixed_frame_rmsd": rmsd_ok,
        "dlg_mapping": mapping_ok,
        "ccd_chemistry": ccd_ok,
        "plip_export": plip_export_ok,
        "protonation_workflow": protonation_ok,
        "provenance": provenance_ok,
        "target_validation_gate": validation_gate_ok,
        "primary_interaction_engine": interaction_primary_ok,
    }
    return {"status": "READY_FOR_CONTROLLED_REAL_STUDY" if all(checks.values()) else "DEVELOPMENT ONLY", "checks": checks, "primary_interaction_engine": PRIMARY_INTERACTION_ENGINE, "sidock_internal_parser_status": SIDOCK_INTERACTION_STATUS}


def interaction_evidence_allowed(source: str) -> bool:
    return source.upper() == PRIMARY_INTERACTION_ENGINE


def assert_primary_interaction_source(source: str) -> None:
    if not interaction_evidence_allowed(source):
        raise ValueError("Only PLIP interaction records may populate primary scientific evidence.")
