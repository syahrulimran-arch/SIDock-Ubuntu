import json
from pathlib import Path

from sidock_ubuntu.executors import ProductionContext, build_production_executors
from sidock_ubuntu.orchestrator import StageRunner, StageResult
from sidock_ubuntu.stage_outputs import (load_grid_output, load_receptor_output,
    load_reference_output, load_validation_output, StageOutputError)
from sidock_ubuntu.state import update_state, project_status

def _fixture(tmp_path, monkeypatch):
    (tmp_path/'receptor.pdb').write_text('ATOM\n')
    (tmp_path/'reference.pdb').write_text('HETATM\n')
    (tmp_path/'receptor.pdbqt').write_text('ATOM\n')
    (tmp_path/'grid.gpf').write_text('gpf\n'); (tmp_path/'grid.glg').write_text('glg\n')
    (tmp_path/'map.map').write_text('map\n'); (tmp_path/'native.dlg').write_text('dlg\n')
    project=tmp_path/'project'; project.mkdir()
    (project/'project.json').write_text(json.dumps({'receptor_source':str(tmp_path/'receptor.pdb'),'ligand_source':str(tmp_path/'ligands.sdf')}))
    ctx=ProductionContext(project,{'validation':{'ga_runs':10}},tmp_path/'receptor.pdb',tmp_path/'ligands.sdf')
    calls=[]
    class Group:
        key='TEST:LIG:1'; resname='LIG'; chain='A'; resseq='1'
    monkeypatch.setattr('sidock_ubuntu.intake.hetero_groups',lambda p:[Group()])
    monkeypatch.setattr('sidock_ubuntu.intake.choose_reference_ligand',lambda p:Group())
    def reference(*a,**k):
        calls.append('REFERENCE'); return {'reference_pdb':tmp_path/'reference.pdb','selected_ligand':'REF_SENTINEL_84731','resname':'LIG','chain':'A','resseq':'1'}
    def receptor(*a,**k):
        calls.append('RECEPTOR'); assert load_reference_output(project).selected_ligand=='REF_SENTINEL_84731'; return {'clean':tmp_path/'receptor.pdb','receptor_pdbqt':tmp_path/'receptor.pdbqt','fallback_metadata':'REC_SENTINEL_49327'}
    def grid(*a,**k):
        calls.append('GRID'); assert load_receptor_output(project).fallback_metadata=='REC_SENTINEL_49327'; return {'gpf':tmp_path/'grid.gpf','glg':tmp_path/'grid.glg','maps':[tmp_path/'map.map'],'center':[1,2,3],'dimensions':[10,10,10],'spacing':.375,'warnings':['GRID_SENTINEL_82154']}
    def validation(*a,**k):
        calls.append('VALIDATION'); assert 'GRID_SENTINEL_82154' in load_grid_output(project).warnings; return {'native_dlg':None,'minimum_rmsd':1.0,'rmsd_threshold':2.0,'scientific_status':'PASS','execution_status':'SUCCESS'}
    monkeypatch.setattr('sidock_ubuntu.executors.run_reference_stage',reference)
    monkeypatch.setattr('sidock_ubuntu.executors.run_receptor_stage',receptor)
    monkeypatch.setattr('sidock_ubuntu.executors.run_grid_stage',grid)
    monkeypatch.setattr('sidock_ubuntu.executors.run_validation_stage',validation)
    monkeypatch.setattr('sidock_ubuntu.validation.prepare_native_ligand',lambda *a,**k:{'pdbqt':tmp_path/'native.pdbqt'})
    monkeypatch.setattr('sidock_ubuntu.executors.box_from_reference',lambda p: object())
    return project,ctx,calls

def test_global_real_registry_persists_and_routes_in_order(tmp_path, monkeypatch):
    project,ctx,calls=_fixture(tmp_path,monkeypatch)
    for stage in ('INPUT','REFERENCE','RECEPTOR','GRID','VALIDATION','FINALIZATION'):
        update_state(project,stage,'PASS' if stage=='INPUT' else 'PENDING')
    result=StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    assert calls==['REFERENCE','RECEPTOR','GRID','VALIDATION'], (result, project_status(project))
    assert [x[1] for x in result]==['PASS']*4, project_status(project)
    state=project_status(project)
    for stage,loader in [('REFERENCE',load_reference_output),('RECEPTOR',load_receptor_output),('GRID',load_grid_output),('VALIDATION',load_validation_output)]:
        assert state['stages'][stage]['result_metadata']
        assert loader(project).load_status=='VALID'

def test_global_restart_after_receptor_runs_only_downstream(tmp_path, monkeypatch):
    project,ctx,calls=_fixture(tmp_path,monkeypatch)
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR'))
    calls.clear()
    from sidock_ubuntu.executors import production_context_from_project
    fresh=production_context_from_project(project,ctx.config)
    StageRunner(project,fresh.config,build_production_executors(fresh)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    assert calls==['GRID','VALIDATION']

def test_global_restart_after_reference_runs_remaining_stages(tmp_path, monkeypatch):
    project,ctx,calls=_fixture(tmp_path,monkeypatch)
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE',))
    calls.clear()
    from sidock_ubuntu.executors import production_context_from_project
    fresh=production_context_from_project(project,ctx.config)
    StageRunner(project,fresh.config,build_production_executors(fresh)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    assert calls==['RECEPTOR','GRID','VALIDATION']

def test_global_restart_after_grid_runs_validation_only(tmp_path, monkeypatch):
    project,ctx,calls=_fixture(tmp_path,monkeypatch)
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID'))
    calls.clear()
    from sidock_ubuntu.executors import production_context_from_project
    fresh=production_context_from_project(project,ctx.config)
    StageRunner(project,fresh.config,build_production_executors(fresh)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    assert calls==['VALIDATION']

def test_global_execution_trace_records_order(tmp_path, monkeypatch):
    project,ctx,_=_fixture(tmp_path,monkeypatch)
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    rows=[json.loads(x) for x in (project/'provenance/execution_trace.jsonl').read_text().splitlines()]
    assert [(r['scope'],r['stage'],r['status']) for r in rows]==[('GLOBAL',s,'PASS') for s in ('REFERENCE','RECEPTOR','GRID','VALIDATION')]

def test_validation_scientific_failure_is_persisted(tmp_path, monkeypatch):
    project,ctx,_=_fixture(tmp_path,monkeypatch)
    def fail(*a,**k): return {'native_dlg':None,'minimum_rmsd':2.5,'rmsd_threshold':2.0,'scientific_status':'FAIL','execution_status':'SUCCESS'}
    monkeypatch.setattr('sidock_ubuntu.executors.run_validation_stage',fail)
    result=StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    assert result[-1][1]=='PASS'
    out=load_validation_output(project); assert out.load_status=='VALID' and out.scientific_status=='FAIL' and out.minimum_rmsd==2.5

def test_validation_scientific_fail_survives_restart_without_rerun(tmp_path, monkeypatch):
    project,ctx,calls=_fixture(tmp_path,monkeypatch)
    def fail(*a,**k): return {'native_dlg':None,'minimum_rmsd':2.5,'rmsd_threshold':2.0,'scientific_status':'FAIL','execution_status':'SUCCESS'}
    monkeypatch.setattr('sidock_ubuntu.executors.run_validation_stage',fail)
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    calls.clear()
    from sidock_ubuntu.executors import production_context_from_project
    fresh=production_context_from_project(project,ctx.config)
    result=StageRunner(project,fresh.config,build_production_executors(fresh)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    assert calls==[] and result[-1][1]=='REUSED'
    out=load_validation_output(project); assert out.load_status=='VALID' and out.scientific_status=='FAIL'

def _persistence_failure_test(tmp_path, monkeypatch, stage, save_name, expected_next, call_index):
    project,ctx,calls=_fixture(tmp_path,monkeypatch)
    monkeypatch.setattr(f'sidock_ubuntu.stage_outputs.{save_name}',lambda *a,**k: (_ for _ in ()).throw(OSError('metadata write failed')))
    order={'REFERENCE':('REFERENCE',),'RECEPTOR':('REFERENCE','RECEPTOR'),'GRID':('REFERENCE','RECEPTOR','GRID'),'VALIDATION':('REFERENCE','RECEPTOR','GRID','VALIDATION')}[stage]
    result=StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(order)
    assert result[-1][1]=='FAIL'
    assert expected_next not in calls
    trace=[json.loads(x) for x in (project/'provenance/execution_trace.jsonl').read_text().splitlines()]
    assert trace[-1]['scope']=='GLOBAL' and trace[-1]['stage']==stage and trace[-1]['status']=='FAIL'
    assert project_status(project)['stages'][stage]['status']=='FAIL'

def test_reference_persistence_failure_blocks_receptor(tmp_path, monkeypatch):
    _persistence_failure_test(tmp_path,monkeypatch,'REFERENCE','save_reference_output','RECEPTOR',0)

def test_receptor_persistence_failure_blocks_grid(tmp_path, monkeypatch):
    _persistence_failure_test(tmp_path,monkeypatch,'RECEPTOR','save_receptor_output','GRID',1)

def test_grid_persistence_failure_blocks_validation(tmp_path, monkeypatch):
    _persistence_failure_test(tmp_path,monkeypatch,'GRID','save_grid_output','VALIDATION',2)

def test_validation_persistence_failure_is_not_reusable(tmp_path, monkeypatch):
    project,ctx,calls=_fixture(tmp_path,monkeypatch)
    monkeypatch.setattr('sidock_ubuntu.stage_outputs.save_validation_output',lambda *a,**k: (_ for _ in ()).throw(OSError('metadata write failed')))
    result=StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    assert result[-1][1]=='FAIL' and 'VALIDATION' in calls
    assert not (project/'validation/validation_stage.json').exists()

def test_validation_pass_unblocks_compound_preparation(tmp_path, monkeypatch):
    project,ctx,_=_fixture(tmp_path,monkeypatch)
    from sidock_ubuntu.state import ensure_compounds
    ensure_compounds(project,['A'])
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    from sidock_ubuntu.orchestrator import StageResult
    calls=[]; runner=StageRunner(project,ctx.config,{'preparation':lambda cid: calls.append(cid) or StageResult('PASS')})
    assert runner.run_compound('A')[0]==('preparation','PASS') and calls==['A']

def test_validation_fail_blocks_compound_preparation_and_restart(tmp_path, monkeypatch):
    project,ctx,_=_fixture(tmp_path,monkeypatch)
    from sidock_ubuntu.state import ensure_compounds
    ensure_compounds(project,['A'])
    monkeypatch.setattr('sidock_ubuntu.executors.run_validation_stage',lambda *a,**k:{'native_dlg':None,'minimum_rmsd':2.5,'rmsd_threshold':2.0,'scientific_status':'FAIL','execution_status':'SUCCESS'})
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    calls=[]; runner=StageRunner(project,ctx.config,{'preparation':lambda cid:calls.append(cid) or StageResult('PASS')})
    assert runner.run_compound('A')==[('preparation','BLOCKED')] and calls==[]
    from sidock_ubuntu.executors import production_context_from_project
    fresh=production_context_from_project(project,ctx.config); calls=[]
    assert StageRunner(project,fresh.config,{'preparation':lambda cid:calls.append(cid) or StageResult('PASS')}).run_compound('A')==[('preparation','BLOCKED')] and calls==[]

def test_stale_grid_blocks_validation_callback(tmp_path, monkeypatch):
    project,ctx,calls=_fixture(tmp_path,monkeypatch)
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID'))
    (tmp_path/'map.map').write_text('changed')
    calls.clear(); result=StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    assert calls==['GRID','VALIDATION']

def _routing_fixture(tmp_path, monkeypatch):
    project,ctx,calls=_fixture(tmp_path,monkeypatch)
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    calls.clear()
    return project,ctx,calls

def test_missing_reference_result_metadata_routing(tmp_path, monkeypatch):
    project,ctx,calls=_routing_fixture(tmp_path,monkeypatch)
    state=project_status(project); state['stages']['REFERENCE'].pop('result_metadata',None)
    from sidock_ubuntu.state import _atomic_write
    _atomic_write(project/'project_state.json',state)
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    assert calls==['REFERENCE','RECEPTOR','GRID','VALIDATION']

def test_missing_reference_metadata_file_routing(tmp_path, monkeypatch):
    project,ctx,calls=_routing_fixture(tmp_path,monkeypatch)
    meta=project/project_status(project)['stages']['REFERENCE']['result_metadata']; (project/meta).unlink()
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    assert calls==['REFERENCE','RECEPTOR','GRID','VALIDATION']

def test_stale_reference_routing_actual_callbacks(tmp_path, monkeypatch):
    project,ctx,calls=_routing_fixture(tmp_path,monkeypatch)
    (tmp_path/'reference.pdb').write_text('changed')
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    assert calls==['REFERENCE','RECEPTOR','GRID','VALIDATION']

def test_stale_receptor_routing_does_not_rerun_reference(tmp_path, monkeypatch):
    project,ctx,calls=_routing_fixture(tmp_path,monkeypatch)
    (tmp_path/'receptor.pdbqt').write_text('changed')
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    assert calls==['RECEPTOR','GRID','VALIDATION']

def test_global_grid_artifact_mutation_is_stale(tmp_path, monkeypatch):
    project,ctx,_=_fixture(tmp_path,monkeypatch)
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID'))
    (tmp_path/'map.map').write_text('changed\n')
    assert load_grid_output(project).load_status=='STALE'

def test_validation_operational_error_is_distinct_and_traced(tmp_path, monkeypatch):
    project,ctx,calls=_fixture(tmp_path,monkeypatch)
    def broken(*a,**k): raise RuntimeError('validation DLG unavailable')
    monkeypatch.setattr('sidock_ubuntu.executors.run_validation_stage',broken)
    result=StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    assert result[-1][1]=='FAIL'
    rows=[json.loads(x) for x in (project/'provenance/execution_trace.jsonl').read_text().splitlines()]
    assert rows[-1]['stage']=='VALIDATION' and rows[-1]['status']=='FAIL'

def test_validation_artifact_hash_mismatch_is_detected(tmp_path, monkeypatch):
    project,ctx,_=_fixture(tmp_path,monkeypatch)
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    # The fixture has no native DLG artifact; add a persisted required artifact
    # through the existing metadata contract and then mutate it.
    from sidock_ubuntu.stage_outputs import ValidationStageOutput, artifact, save_validation_output
    native=tmp_path/'native.dlg'; native.write_text('original')
    save_validation_output(project,ValidationStageOutput(str(native),1.0,None,None,2.0,'PASS','SUCCESS',artifacts=[artifact(project,native)]))
    native.write_text('changed')
    assert load_validation_output(project).load_status=='STALE'

def test_incomplete_metadata_is_controlled(tmp_path):
    p=tmp_path/'validation'; p.mkdir(parents=True)
    (p/'reference_stage.json').write_text(json.dumps({'schema_version':'0.4-global-stage-1','result':{'selected_ligand':'x'}}))
    with __import__('pytest').raises(StageOutputError): load_reference_output(tmp_path)

def test_rehydration_has_no_scientific_side_effect(tmp_path, monkeypatch):
    project,ctx,_=_fixture(tmp_path,monkeypatch)
    StageRunner(project,ctx.config,build_production_executors(ctx)).run_global(('REFERENCE','RECEPTOR','GRID','VALIDATION'))
    for name in ('run_reference_stage','run_receptor_stage','run_grid_stage','run_validation_stage'):
        monkeypatch.setattr(f'sidock_ubuntu.executors.{name}',lambda *a,**k: (_ for _ in ()).throw(RuntimeError('science called')))
    from sidock_ubuntu.executors import production_context_from_project
    production_context_from_project(project,ctx.config).rehydrate_global_outputs()
