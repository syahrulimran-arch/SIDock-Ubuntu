from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

WATER = {"HOH", "WAT", "H2O", "DOD"}
COMMON_IONS = {"NA", "K", "CL", "CA", "MG", "ZN", "MN", "FE", "CU", "CO", "NI", "CD", "HG", "BR", "IOD", "I", "LI", "CS", "SR", "BA", "AL", "AG", "AU", "PB"}
COMMON_BUFFER = {"SO4", "PO4", "NO3", "ACT", "ACE", "FMT", "GOL", "EDO", "PEG", "MPD", "TRS", "MES", "HEP", "BME", "DMS", "DTT", "NAG", "MAN"}


@dataclass(frozen=True)
class HeteroGroup:
    resname: str
    chain: str
    resseq: str
    atoms: tuple[str, ...]

    @property
    def atom_count(self) -> int:
        return len(self.atoms)

    @property
    def key(self) -> str:
        return f"{self.resname}:{self.chain or '-'}:{self.resseq}"


def hetero_groups(pdb: Path) -> list[HeteroGroup]:
    grouped: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("HETATM"):
            continue
        key = (line[17:20].strip().upper(), line[21:22].strip(), line[22:27].strip())
        grouped[key].append(line)
    return [HeteroGroup(r, c, n, tuple(lines)) for (r, c, n), lines in grouped.items()]


def choose_reference_ligand(pdb: Path) -> HeteroGroup | None:
    candidates = []
    for group in hetero_groups(pdb):
        if group.resname in WATER or group.resname in COMMON_IONS or group.resname in COMMON_BUFFER:
            continue
        if group.atom_count < 6:
            continue
        candidates.append(group)
    return max(candidates, key=lambda g: g.atom_count) if candidates else None


def write_group(group: HeteroGroup, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(group.atoms) + "\nEND\n", encoding="utf-8")
    return destination


def protein_intake_summary(pdb: Path) -> dict[str, object]:
    groups = hetero_groups(pdb)
    ref = choose_reference_ligand(pdb)
    return {
        "pdb": str(pdb.resolve()),
        "hetero_groups": [{"key": g.key, "resname": g.resname, "chain": g.chain, "resseq": g.resseq, "atoms": g.atom_count} for g in groups],
        "selected_reference_ligand": ref.key if ref else None,
        "reference_ligand_atoms": ref.atom_count if ref else 0,
    }
