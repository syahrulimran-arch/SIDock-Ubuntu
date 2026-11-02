from dataclasses import dataclass

ERROR_CODES = {
    "INPUT_FILE_MISSING": "E001", "INVALID_SDF": "E002", "DUPLICATE_COMPOUND_ID": "E003",
    "REFERENCE_LIGAND_UNRESOLVED": "E101", "REFERENCE_CHEMISTRY_FAILED": "E102",
    "RECEPTOR_PREPARATION_FAILED": "E201", "BINDING_SITE_RESIDUE_INCOMPLETE": "E202",
    "NATIVE_REDOCKING_FAILED": "E301", "VALIDATION_RMSD_EXCEEDED": "E302",
    "LIGAND_PREPARATION_FAILED": "E401", "AUTOGRID_FAILED": "E501", "AUTODOCK_FAILED": "E502",
    "DLG_PARSE_FAILED": "E601", "POSE_MAPPING_FAILED": "E602", "PLIP_FAILED": "E701",
}

@dataclass
class SIDockError(RuntimeError):
    code_name: str
    message: str
    details: dict | None = None
    def __str__(self):
        return f"{ERROR_CODES.get(self.code_name, 'E999')} {self.code_name}: {self.message}"

def error_payload(error: SIDockError) -> dict:
    return {"code": ERROR_CODES.get(error.code_name, "E999"), "error": error.code_name,
            "message": error.message, "details": error.details or {}}
