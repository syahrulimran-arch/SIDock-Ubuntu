from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass
class GridBox:
    center_x: float
    center_y: float
    center_z: float
    npts_x: int
    npts_y: int
    npts_z: int
    spacing: float = 0.375


@dataclass
class DockingParameters:
    ga_runs: int = 10
    ga_num_evals: int = 2_500_000
    ga_pop_size: int = 150
    ga_num_generations: int = 27_000
    rmstol: float = 2.0
    parallel_jobs: int = 1


@dataclass
class Project:
    root: str
    receptor_source: str = ""
    receptor_pdbqt: str = ""
    ligand_source: str = ""
    prepared_ligands: list[str] = field(default_factory=list)
    grid: GridBox | None = None
    docking: DockingParameters = field(default_factory=DockingParameters)

    @property
    def path(self) -> Path:
        return Path(self.root)

    def ensure_layout(self) -> None:
        for rel in (
            "receptor", "ligands/source", "ligands/prepared", "grid",
            "docking/dpf", "docking/dlg", "docking/work", "results",
            "analysis", "analysis/poses", "analysis/complexes", "analysis/figures", "logs",
        ):
            (self.path / rel).mkdir(parents=True, exist_ok=True)

    def save(self) -> None:
        self.ensure_layout()
        (self.path / "project.json").write_text(json.dumps(asdict(self), indent=2), encoding="utf-8")
