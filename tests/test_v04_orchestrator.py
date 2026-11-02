from sidock_ubuntu.orchestrator import StageRunner, StageResult, config_change_audit
from sidock_ubuntu.state import ensure_compounds, update_compound, project_status
from sidock_ubuntu.config import DEFAULT_CONFIG, deep_merge
from sidock_ubuntu.executors import ProductionContext, execute_plip, execute_compound_docking

def setup(tmp_path):
    ensure_compounds(tmp_path,['A','B','C','D'])
    for s in ('preparation','docking','pose_extraction','clustering'):
        update_compound(tmp_path,'A',s,'PASS')
    for s in ('preparation','docking','pose_extraction','clustering'):
        update_compound(tmp_path,'B',s,'PASS')
    update_compound(tmp_path,'B','plip','FAIL',error_code='E701')
    update_compound(tmp_path,'C','preparation','PASS'); update_compound(tmp_path,'C','docking','FAIL',error_code='E502')
    return tmp_path

def test_plan_resumes_plip_only(tmp_path):
    p=setup(tmp_path); plan=StageRunner(p,DEFAULT_CONFIG).plan(['B'])
    assert plan['B']==['plip','figures']

def test_failed_only_excludes_pass(tmp_path):
    p=setup(tmp_path); plan=StageRunner(p,DEFAULT_CONFIG).plan(['B'],failed_only=True)
    assert plan['B']==['plip']

def test_compound_filter_excludes_others(tmp_path):
    p=setup(tmp_path); plan=StageRunner(p,DEFAULT_CONFIG).plan(['B'])
    assert 'A' not in plan

def test_force_includes_completed_stages(tmp_path):
    p=setup(tmp_path); plan=StageRunner(p,DEFAULT_CONFIG).plan(['A'],force=True)
    assert plan['A']==['preparation','docking','pose_extraction','clustering','plip','figures']

def test_executor_runs_and_records_attempt(tmp_path):
    p=setup(tmp_path); calls=[]
    def plip(cid): calls.append(cid); return StageResult('PASS')
    runner=StageRunner(p,DEFAULT_CONFIG,{'plip':plip,'figures':lambda cid: StageResult('PASS')})
    runner.run_compound('B',failed_only=True)
    assert calls==['B']; assert project_status(p)['compounds']['B']['plip']['status']=='PASS'

def test_failed_autodock_does_not_run_plip(tmp_path):
    p=setup(tmp_path); calls=[]
    runner=StageRunner(p,DEFAULT_CONFIG,{'docking':lambda cid: calls.append(cid) or StageResult('PASS'), 'pose_extraction':lambda cid: StageResult('PASS'), 'clustering':lambda cid: StageResult('PASS'), 'plip':lambda cid: StageResult('PASS'), 'figures':lambda cid: StageResult('PASS')})
    runner.run_compound('C',failed_only=True); assert calls==['C']

def test_attempt_history_is_preserved_on_retry(tmp_path):
    p=setup(tmp_path); runner=StageRunner(p,DEFAULT_CONFIG,{'plip':lambda cid: StageResult('PASS'), 'figures':lambda cid: StageResult('PASS')})
    runner.run_compound('B',failed_only=True)
    assert len(project_status(p)['compounds']['B']['plip']['attempts'])>=2

def test_config_audit_docking_change(tmp_path):
    a=DEFAULT_CONFIG; b=deep_merge(a,{'docking':{'ga_runs':20}}); out=config_change_audit(tmp_path,a,b)
    assert 'DOCKING' in out['invalidated_stages'] and 'PLIP' in out['invalidated_stages']

def test_config_audit_interaction_change(tmp_path):
    a=DEFAULT_CONFIG; b=deep_merge(a,{'interaction':{'primary_engine':'plip'}}); out=config_change_audit(tmp_path,a,b)
    assert out['changed'] is False

def test_config_audit_figures_change(tmp_path):
    a=DEFAULT_CONFIG; b=deep_merge(a,{'figures':{'dpi':300}}); out=config_change_audit(tmp_path,a,b)
    assert out['invalidated_stages']==['FIGURES','FINALIZATION']

def test_dry_run_audit_written(tmp_path):
    out=config_change_audit(tmp_path,DEFAULT_CONFIG,deep_merge(DEFAULT_CONFIG,{'validation':{'rmsd_threshold':1.0}}))
    assert (tmp_path/'provenance/config_change_audit.json').exists() and 'VALIDATION' in out['invalidated_stages']

def test_pending_compound_plan_starts_preparation(tmp_path):
    setup(tmp_path); plan=StageRunner(tmp_path,DEFAULT_CONFIG).plan(['D'])
    assert plan['D'][0]=='preparation'

def test_same_config_audit_is_empty(tmp_path):
    out=config_change_audit(tmp_path,DEFAULT_CONFIG,DEFAULT_CONFIG)
    assert out['changed'] is False and out['invalidated_stages']==[]

def test_atomic_state_keeps_previous_snapshot(tmp_path):
    from sidock_ubuntu.state import update_state
    update_state(tmp_path,'INPUT','PASS'); update_state(tmp_path,'INPUT','PASS')
    assert (tmp_path/'project_state.previous.json').exists()

def test_production_executor_delegates_plip_callback(tmp_path):
    calls=[]; ctx=ProductionContext(tmp_path,DEFAULT_CONFIG,callbacks={'plip':lambda cid: calls.append(cid) or {'interactions':'x'}})
    assert execute_plip(ctx,'B').status=='PASS' and calls==['B']

def test_production_executor_does_not_invent_docking_backend(tmp_path):
    ctx=ProductionContext(tmp_path,DEFAULT_CONFIG)
    result=execute_compound_docking(ctx,'B')
    assert result.status=='BLOCKED' and result.error_code=='E502'

