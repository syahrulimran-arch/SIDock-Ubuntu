"""Verify user-provided AD4 parameter files reach BOTH GPF and DPF.

No AutoDock4 or AutoGrid4 binaries are launched by these tests.
"""
import hashlib
from pathlib import Path

from sidock_ubuntu.model import DockingParameters, GridBox
from sidock_ubuntu.stages import run_compound_docking, run_grid_stage
from sidock_ubuntu import stages


def test_custom_parameters_appear_in_autogrid_gpf(monkeypatch, tmp_path):
    receptor = tmp_path / "input_receptor.pdbqt"
    receptor.write_text(
        "ATOM      1  C   GLY C  20      10.000  20.000  30.000  1.00  0.00    +0.100 C\n"
    )
    parameter = tmp_path / "custom_parameters.dat"
    parameter.write_text("atom_par C 4.00 0.150 33.00 -0.001 0 0 0 -1 -1 0\n")
    monkeypatch.setattr(stages, "run_command", lambda *a, **k: 0)

    box = GridBox(10., 20., 30., 80, 80, 80, 0.375)
    output = run_grid_stage(receptor, tmp_path/"grid", box, ["C"], "autogrid4",
                            parameter_file=parameter)
    gpf = (tmp_path/"grid/receptor.gpf").read_text()
    assert gpf.splitlines()[0] == "parameter_file custom_ad4_parameters.dat"
    assert "gridcenter 10.000 20.000 30.000" in gpf
    assert (tmp_path/"grid/custom_ad4_parameters.dat").read_bytes() == parameter.read_bytes()
    import json
    provenance = json.loads(Path(output["grid_provenance"]).read_text())
    assert provenance["custom_parameter_sha256"] == hashlib.sha256(parameter.read_bytes()).hexdigest()


def test_custom_parameters_appear_in_autodock_dpf(monkeypatch, tmp_path):
    root = tmp_path/"run"
    dpf_dir = root/"docking/dpf"
    dlg_dir = root/"docking/dlg"
    work = root/"docking/work"
    ligand = tmp_path/"A1.pdbqt"
    ligand.write_text(
        "ROOT\nHETATM    1  C1  UNL A   1      10.000  20.000  30.000  1.00  0.00    +0.100 C\n"
        "ENDROOT\nTORSDOF 0\n"
    )
    params = tmp_path/"parameters.dat"
    params.write_text("atom_par C 4.00 0.150 33.00 -0.001 0 0 0 -1 -1 0\n")
    grid = root/"grid/production"
    grid.mkdir(parents=True)
    (grid/"receptor.maps.fld").write_text("fld")
    (grid/"receptor.C.map").write_text("map")
    def fake_run(command, cwd, log):
        Path(command[command.index("-l")+1]).write_text("Estimated Free Energy of Binding = -5.000\n")
        return 0
    monkeypatch.setattr(stages, "run_command", fake_run)

    response = run_compound_docking(
        ligand, "receptor", work, dpf_dir, dlg_dir, ["C"],
        DockingParameters(ga_runs=1), "autodock4", parameter_file=params
    )
    dpf = Path(response["dpf"]).read_text()
    assert dpf.splitlines()[0] == "parameter_file custom_ad4_parameters.dat"
    assert "ga_run 1" in dpf
    assert (work/"A1/custom_ad4_parameters.dat").read_bytes() == params.read_bytes()
