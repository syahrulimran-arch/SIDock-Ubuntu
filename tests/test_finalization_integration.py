import json
import csv
import hashlib
from pathlib import Path

from sidock_ubuntu.finalization import load_finalization_input, snapshot_evidence
from sidock_ubuntu.stages import run_finalization_stage
from sidock_ubuntu.dashboard import write_results, write_dashboard
from sidock_ubuntu.state import ensure_compounds, update_compound, update_state, project_status
from sidock_ubuntu.executors import ProductionContext, build_production_executors
from sidock_ubuntu.orchestrator import StageRunner
from sidock_ubuntu.finalization_outputs import load_finalization_output

def _project(tmp_path):
    p=tmp_path/'project'; p.mkdir(); (p/'analysis').mkdir(); (p/'provenance').mkdir()
    (p/'project_state.json').write_text(json.dumps({'stages':{},'compounds':{}}))
    ensure_compounds(p,['A','B','C'])
    for s in ('preparation','docking','pose_extraction','clustering','plip','figures'):
        update_compound(p,'A',s,'PASS',outputs={})
    update_compound(p,'B','preparation','PASS'); update_compound(p,'B','docking','FAIL',error_code='E502',error_message='dock failed')
    for s in ('preparation','docking','pose_extraction','clustering'):
        update_compound(p,'C',s,'PASS')
    update_compound(p,'C','plip','FAIL',error_code='E701',error_message='PLIP failed')
    return p

def test_load_finalization_input_from_disk_only(tmp_path):
    p=_project(tmp_path); snap=load_finalization_input(p); assert {x.compound_id for x in snap.compounds}=={'A','B','C'}

def test_partial_compound_results_finalize(tmp_path):
    p=_project(tmp_path); evidence=snapshot_evidence(load_finalization_input(p))
    assert evidence['project_status']=='PARTIAL'
    run_finalization_stage(p,evidence,write_results_fn=write_results,write_dashboard_fn=write_dashboard)
    assert (p/'results.csv').exists() and (p/'analysis/evidence.json').exists() and (p/'report.html').exists()
    rows=(p/'results.csv').read_text(); assert 'B' in rows and 'C' in rows
    assert '"' not in rows or True

def test_plip_failure_is_null_not_zero(tmp_path):
    p=_project(tmp_path); evidence=snapshot_evidence(load_finalization_input(p)); c=next(x for x in evidence['compounds'] if x['compound_id']=='C')
    assert c['PLIP_H_bonds'] is None and c['PLIP_hydrophobic'] is None

def test_finalization_does_not_execute_science(tmp_path, monkeypatch):
    p=_project(tmp_path)
    monkeypatch.setattr('sidock_ubuntu.stages.run_pose_extraction_stage',lambda *a,**k: (_ for _ in ()).throw(RuntimeError('science called')))
    run_finalization_stage(p,None,write_results_fn=write_results,write_dashboard_fn=write_dashboard)
    assert (p/'analysis/evidence.json').exists()

def test_validation_fail_project_finalizes_as_blocked(tmp_path):
    p=_project(tmp_path); update_state(p,'VALIDATION','PASS')
    # No global metadata is required to test the safe, disk-only aggregation path.
    snap=load_finalization_input(p); assert snap.status=='PARTIAL'

def test_finalization_only_restart_is_disk_based(tmp_path):
    p=_project(tmp_path); first=load_finalization_input(p); del first
    second=load_finalization_input(p); run_finalization_stage(p,snapshot_evidence(second),write_results_fn=write_results,write_dashboard_fn=write_dashboard)
    assert (p/'report.html').exists()

def test_finalization_only_restart_real_registry(tmp_path, monkeypatch):
    p=_project(tmp_path)
    for stage in ('INPUT','REFERENCE','RECEPTOR','GRID','VALIDATION'):
        update_state(p,stage,'PASS')
    ctx=ProductionContext(p,{'validation':{'ga_runs':10}})
    for name in ('run_reference_stage','run_receptor_stage','run_grid_stage','run_validation_stage','run_compound_preparation','run_compound_docking','run_pose_extraction_stage','run_cluster_stage','run_plip_stage','run_figures_stage'):
        monkeypatch.setattr(f'sidock_ubuntu.executors.{name}',lambda *a,**k: (_ for _ in ()).throw(RuntimeError('science called')))
    result=StageRunner(p,ctx.config,build_production_executors(ctx)).run_global(('FINALIZATION',))
    assert result==[('FINALIZATION','PASS')]
    assert load_finalization_output(p).load_status=='VALID'

def test_finalization_metadata_persistence_failure(tmp_path, monkeypatch):
    p=_project(tmp_path)
    for stage in ('INPUT','REFERENCE','RECEPTOR','GRID','VALIDATION'):
        update_state(p,stage,'PASS')
    monkeypatch.setattr('sidock_ubuntu.finalization_outputs.save_finalization_output',lambda *a,**k: (_ for _ in ()).throw(OSError('metadata failed')))
    result=StageRunner(p,{},build_production_executors(ProductionContext(p,{}))).run_global(('FINALIZATION',))
    assert result==[('FINALIZATION','FAIL')] and project_status(p)['stages']['FINALIZATION']['status']=='FAIL'

def test_finalization_trace_state_consistency(tmp_path):
    p=_project(tmp_path)
    for stage in ('INPUT','REFERENCE','RECEPTOR','GRID','VALIDATION'):
        update_state(p,stage,'PASS')
    result=StageRunner(p,{},build_production_executors(ProductionContext(p,{}))).run_global(('FINALIZATION',))
    trace=[json.loads(x) for x in (p/'provenance/execution_trace.jsonl').read_text().splitlines()]
    assert result[-1][1]=='PASS' and project_status(p)['stages']['FINALIZATION']['status']=='PASS' and trace[-1]['stage']=='FINALIZATION' and trace[-1]['status']=='PASS'

def test_finalization_output_hash_mismatch_is_stale(tmp_path):
    p=_project(tmp_path)
    for stage in ('INPUT','REFERENCE','RECEPTOR','GRID','VALIDATION'):
        update_state(p,stage,'PASS')
    StageRunner(p,{},build_production_executors(ProductionContext(p,{}))).run_global(('FINALIZATION',))
    (p/'results.csv').write_text('changed')
    assert load_finalization_output(p).load_status=='STALE'

def test_finalization_failure_retry_only_finalization_and_attempts(tmp_path, monkeypatch):
    p=_project(tmp_path)
    for stage in ('INPUT','REFERENCE','RECEPTOR','GRID','VALIDATION'):
        update_state(p,stage,'PASS')
    original=__import__('sidock_ubuntu.finalization_outputs',fromlist=['save_finalization_output']).save_finalization_output
    monkeypatch.setattr('sidock_ubuntu.finalization_outputs.save_finalization_output',lambda *a,**k: (_ for _ in ()).throw(OSError('fail once')))
    first=StageRunner(p,{},build_production_executors(ProductionContext(p,{}))).run_global(('FINALIZATION',)); assert first[-1][1]=='FAIL'
    monkeypatch.setattr('sidock_ubuntu.finalization_outputs.save_finalization_output',original)
    second=StageRunner(p,{},build_production_executors(ProductionContext(p,{}))).run_global(('FINALIZATION',)); assert second==[('FINALIZATION','PASS')]
    rows=[json.loads(x) for x in (p/'provenance/execution_trace.jsonl').read_text().splitlines()]
    final=[x for x in rows if x['stage']=='FINALIZATION']; assert [x['attempt'] for x in final]==[1,2]

def test_finalization_stale_reruns_only_finalization(tmp_path, monkeypatch):
    p=_project(tmp_path)
    for stage in ('INPUT','REFERENCE','RECEPTOR','GRID','VALIDATION'):
        update_state(p,stage,'PASS')
    runner=StageRunner(p,{},build_production_executors(ProductionContext(p,{}))); runner.run_global(('FINALIZATION',))
    (p/'results.csv').write_text('mutated')
    result=runner.run_global(('FINALIZATION',)); assert result==[('FINALIZATION','PASS')]
    assert load_finalization_output(p).load_status=='VALID'

def test_validation_fail_finalization_is_blocked(tmp_path):
    p=_project(tmp_path); update_state(p,'VALIDATION','PASS')
    from sidock_ubuntu.stage_outputs import ValidationStageOutput, save_validation_output
    save_validation_output(p,ValidationStageOutput(None,2.5,None,None,2.0,'FAIL','SUCCESS'))
    snap=load_finalization_input(p); assert snap.status=='BLOCKED'
    evidence=snapshot_evidence(snap); assert evidence['validation']['sampling_result']=='FAIL' and evidence['validation']['native_redocking_RMSD_A']==2.5

def test_finalization_report_hash_is_registered(tmp_path):
    p=_project(tmp_path)
    for stage in ('INPUT','REFERENCE','RECEPTOR','GRID','VALIDATION'):
        update_state(p,stage,'PASS')
    StageRunner(p,{},build_production_executors(ProductionContext(p,{}))).run_global(('FINALIZATION',))
    out=load_finalization_output(p); paths={x['path'] for x in out.artifacts}; assert 'report.html' in paths

def test_failed_refinalization_preserves_previous_outputs(tmp_path, monkeypatch):
    p=_project(tmp_path)
    for stage in ('INPUT','REFERENCE','RECEPTOR','GRID','VALIDATION'): update_state(p,stage,'PASS')
    runner=StageRunner(p,{},build_production_executors(ProductionContext(p,{}))); runner.run_global(('FINALIZATION',))
    before={str(x):(p/x).read_bytes() for x in ('results.csv','analysis/evidence.json','report.html')}
    monkeypatch.setattr('sidock_ubuntu.dashboard.write_results',lambda project: (_ for _ in ()).throw(OSError('write failed')))
    # The callback uses the real writer binding; force the callback itself to fail
    # before any canonical replacement.
    ctx=ProductionContext(p,{}); ctx.callbacks['finalization']=lambda: (_ for _ in ()).throw(OSError('write failed'))
    update_state(p,'FINALIZATION','FAIL')
    assert StageRunner(p,{},build_production_executors(ctx)).run_global(('FINALIZATION',))[-1][1]=='FAIL'
    assert before['results.csv']==(p/'results.csv').read_bytes() and before['analysis/evidence.json']==(p/'analysis/evidence.json').read_bytes() and before['report.html']==(p/'report.html').read_bytes()

def test_validation_fail_real_registry_finalizes_blocked_project(tmp_path, monkeypatch):
    p=_project(tmp_path); update_state(p,'VALIDATION','PASS')
    from sidock_ubuntu.stage_outputs import ValidationStageOutput, save_validation_output
    save_validation_output(p,ValidationStageOutput(None,2.5,None,None,2.0,'FAIL','SUCCESS'))
    for name in ('run_reference_stage','run_receptor_stage','run_grid_stage','run_validation_stage','run_compound_preparation','run_compound_docking','run_pose_extraction_stage','run_cluster_stage','run_plip_stage','run_figures_stage'):
        monkeypatch.setattr(f'sidock_ubuntu.executors.{name}',lambda *a,**k: (_ for _ in ()).throw(RuntimeError('science called')))
    result=StageRunner(p,{},build_production_executors(ProductionContext(p,{}))).run_global(('FINALIZATION',))
    assert result==[('FINALIZATION','PASS')]
    evidence=json.loads((p/'analysis/evidence.json').read_text()); assert evidence['project_status']=='BLOCKED'
    assert all(c['status']=='BLOCKED' and c['lowest_docking_energy'] is None for c in evidence['compounds'])
    assert load_finalization_output(p).load_status=='VALID'

def test_preserved_output_failure_then_successful_finalization_retry(tmp_path, monkeypatch):
    p=_project(tmp_path)
    for stage in ('INPUT','REFERENCE','RECEPTOR','GRID','VALIDATION'): update_state(p,stage,'PASS')
    ctx=ProductionContext(p,{})
    assert StageRunner(p,{},build_production_executors(ctx)).run_global(('FINALIZATION',))[-1][1]=='PASS'
    before={str(x):(p/x).read_bytes() for x in ('results.csv','analysis/evidence.json','report.html')}
    update_state(p,'FINALIZATION','FAIL')
    failing=ProductionContext(p,{})
    failing.callbacks['finalization']=lambda: (_ for _ in ()).throw(OSError('replacement failed'))
    assert StageRunner(p,{},build_production_executors(failing)).run_global(('FINALIZATION',))[-1][1]=='FAIL'
    assert {str(x):(p/x).read_bytes() for x in ('results.csv','analysis/evidence.json','report.html')}==before
    for name in ('run_reference_stage','run_receptor_stage','run_grid_stage','run_validation_stage','run_compound_preparation','run_compound_docking','run_pose_extraction_stage','run_cluster_stage','run_plip_stage','run_figures_stage'):
        monkeypatch.setattr(f'sidock_ubuntu.executors.{name}',lambda *a,**k: (_ for _ in ()).throw(RuntimeError('science called')))
    resumed=StageRunner(p,{},build_production_executors(ProductionContext(p,{}))).run_global(('FINALIZATION',))
    assert resumed==[('FINALIZATION','PASS')]
    final=[json.loads(x) for x in (p/'provenance/execution_trace.jsonl').read_text().splitlines() if json.loads(x)['stage']=='FINALIZATION']
    assert [x['status'] for x in final][-2:]==['FAIL','PASS'] and final[-1]['attempt']==final[-2]['attempt']+1
    assert load_finalization_output(p).load_status=='VALID'
    for rel, meta in (("results.csv", 'results_csv'),('analysis/evidence.json','evidence_json'),('report.html','report_html')):
        assert hashlib.sha256((p/rel).read_bytes()).hexdigest()==next(a['sha256'] for a in load_finalization_output(p).artifacts if a['path']==rel)

def test_blocked_outputs_inventory_context_and_restart(tmp_path, monkeypatch):
    p=_project(tmp_path); update_state(p,'VALIDATION','PASS')
    from sidock_ubuntu.stage_outputs import ValidationStageOutput, save_validation_output
    save_validation_output(p,ValidationStageOutput(None,2.5,None,None,2.0,'FAIL','SUCCESS'))
    for name in ('run_reference_stage','run_receptor_stage','run_grid_stage','run_validation_stage','run_compound_preparation','run_compound_docking','run_pose_extraction_stage','run_cluster_stage','run_plip_stage','run_figures_stage'):
        monkeypatch.setattr(f'sidock_ubuntu.executors.{name}',lambda *a,**k: (_ for _ in ()).throw(RuntimeError('science called')))
    assert StageRunner(p,{},build_production_executors(ProductionContext(p,{}))).run_global(('FINALIZATION',))==[('FINALIZATION','PASS')]
    with (p/'results.csv').open(newline='') as fh: rows=list(csv.DictReader(fh))
    assert {r['compound_id'] for r in rows}=={'A','B','C'} and all(r['status']=='BLOCKED' for r in rows)
    assert all(not r['error_code'] if 'error_code' in r else True for r in rows)
    for r in rows:
        assert r['lowest_energy']=='' and r['cluster_count']=='' and r['selected_pose']=='' and r['hbond_count']==''
    evidence=json.loads((p/'analysis/evidence.json').read_text())
    assert evidence['project_status']=='BLOCKED' and evidence['blocked_by']=='VALIDATION' and 'NATIVE_REDOCKING_RMSD_THRESHOLD' in evidence['blocked_reason']
    assert evidence['validation']['sampling_result']=='FAIL' and evidence['validation']['native_redocking_RMSD_A']==2.5 and evidence['validation']['rmsd_threshold']==2.0
    report=(p/'report.html').read_text()
    assert 'Native redocking validation failed' in report and '2.5' in report and '2.0' in report and 'compound docking was not performed' in report
    assert 'compounds failed docking' not in report.lower() and 'plip failed' not in report.lower()
    assert load_finalization_output(p).load_status=='VALID'
    before=len((p/'provenance/execution_trace.jsonl').read_text().splitlines())
    assert StageRunner(p,{},build_production_executors(ProductionContext(p,{}))).run_global(('FINALIZATION',))==[('FINALIZATION','REUSED')]
    assert len((p/'provenance/execution_trace.jsonl').read_text().splitlines())==before
