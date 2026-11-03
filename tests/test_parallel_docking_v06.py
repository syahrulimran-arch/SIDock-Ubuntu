import json
import threading
import time
from pathlib import Path

from sidock_ubuntu.config import DEFAULT_CONFIG, validate_config
from sidock_ubuntu.model import DockingParameters
from sidock_ubuntu.orchestrator import StageResult, StageRunner
from sidock_ubuntu.state import COMPOUND_STAGES, ensure_compounds, project_status
from sidock_ubuntu.stages import run_compound_docking


def test_parallel_workers_auto_uses_available_cpus_but_respects_ten_job_ceiling(tmp_path, monkeypatch):
    monkeypatch.setattr('sidock_ubuntu.orchestrator.available_cpu_count', lambda: 32)
    runner=StageRunner(tmp_path, {'validation':{'mode':'skip'}, 'docking':{'parallel_jobs':'auto','max_parallel_jobs':10,'ga_runs':25}}, {})
    assert runner._parallel_workers(10) == 10
    assert runner._parallel_workers(40) == 10


def test_parallel_workers_manual_override_is_capped_by_available_cpu_compounds_and_ceiling(tmp_path, monkeypatch):
    monkeypatch.setattr('sidock_ubuntu.orchestrator.available_cpu_count', lambda: 16)
    runner=StageRunner(tmp_path, {'validation':{'mode':'skip'}, 'docking':{'parallel_jobs':16,'max_parallel_jobs':10,'ga_runs':25}}, {})
    assert runner._parallel_workers(20) == 10
    assert runner._parallel_workers(4) == 4
    runner=StageRunner(tmp_path, {'validation':{'mode':'skip'}, 'docking':{'parallel_jobs':8,'max_parallel_jobs':10,'ga_runs':25}}, {})
    assert runner._parallel_workers(20) == 8


def test_compound_docking_batches_autodock_concurrently_and_persists_ga_runs(tmp_path, monkeypatch):
    monkeypatch.setattr('sidock_ubuntu.orchestrator.available_cpu_count', lambda: 4)
    ensure_compounds(tmp_path, ['A','B','C','D'])

    active=0
    max_active=0
    guard=threading.Lock()

    def prepare(cid):
        return StageResult('PASS', {})

    def dock(cid):
        nonlocal active, max_active
        with guard:
            active += 1
            max_active=max(max_active,active)
        time.sleep(0.04)
        with guard:
            active -= 1
        return StageResult('PASS', {})

    executors={
        'preparation':prepare,
        'docking':dock,
        'pose_extraction':lambda cid: StageResult('PASS', {}),
        'clustering':lambda cid: StageResult('PASS', {}),
        'plip':lambda cid: StageResult('PASS', {}),
        'figures':lambda cid: StageResult('PASS', {}),
    }
    config={'validation':{'mode':'skip'}, 'docking':{'parallel_jobs':'auto','max_parallel_jobs':10,'ga_runs':37}}
    runner=StageRunner(tmp_path,config,executors)
    runner._run_compound_unlocked('A')

    assert max_active >= 2
    state=project_status(tmp_path)
    assert all(state['compounds'][cid]['docking']['status']=='PASS' for cid in ('A','B','C','D'))
    settings=json.loads((tmp_path/'provenance'/'parallel_execution.json').read_text())
    assert settings['detected_cpus']==4
    assert settings['max_parallel_jobs']==10
    assert settings['effective_parallel_workers']==4
    assert settings['ga_runs_per_compound']==37
    assert settings['parallelized_stage']=='docking'


def test_run_compound_docking_uses_independent_ga_setting(tmp_path, monkeypatch):
    ligand=tmp_path/'ligands'/'prepared'/'A.pdbqt'; ligand.parent.mkdir(parents=True)
    ligand.write_text('ATOM      1  C   LIG A   1       1.000   2.000   3.000  1.00  0.00      0.000 C\nTORSDOF 0\n')
    work=tmp_path/'docking'/'work'; work.mkdir(parents=True)
    (work/'receptor.C.map').write_text('map\n')
    (work/'receptor.e.map').write_text('map\n')
    (work/'receptor.d.map').write_text('map\n')
    (work/'receptor.maps.fld').write_text('fld\n')
    prov=tmp_path/'provenance'; prov.mkdir()
    (prov/'parallel_execution.json').write_text(json.dumps({'ga_runs_per_compound':77}))

    captured={}
    def fake_write_dpf(receptor_stem, local_ligand, ligand_types, params, destination):
        captured['ga_runs']=params.ga_runs
        destination.parent.mkdir(parents=True,exist_ok=True)
        destination.write_text(f'ga_run {params.ga_runs}\n')

    def fake_run_command(cmd,cwd,log):
        Path(cmd[4]).write_text('DLG\n')
        return 0

    monkeypatch.setattr('sidock_ubuntu.stages.write_dpf',fake_write_dpf)
    monkeypatch.setattr('sidock_ubuntu.stages.run_command',fake_run_command)
    original=DockingParameters(ga_runs=10)
    result=run_compound_docking(ligand,'receptor',work,tmp_path/'docking'/'dpf',tmp_path/'docking'/'dlg',['C'],original,'autodock4')

    assert captured['ga_runs']==77
    assert original.ga_runs==10
    assert result['ga_runs']==77
    assert (tmp_path/'docking'/'work'/'A').is_dir()


def test_config_keeps_ga_runs_and_parallel_jobs_independent():
    config=json.loads(json.dumps(DEFAULT_CONFIG))
    config['validation']['ga_runs']=7
    config['docking']['ga_runs']=83
    config['docking']['parallel_jobs']='auto'
    config['docking']['max_parallel_jobs']=10
    assert validate_config(config)==[]
    assert config['validation']['ga_runs']==7
    assert config['docking']['ga_runs']==83
    assert config['docking']['max_parallel_jobs']==10
