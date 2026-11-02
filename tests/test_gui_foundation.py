import json
from pathlib import Path
import pytest
import multiprocessing, time

from sidock_ubuntu.gui_service import create_project, project_view, safe_project_name, workspace_path
from sidock_ubuntu.gui_service import get_compounds, get_compound_detail, get_project_results, resolve_project_artifact, get_doctor_status, get_config, get_project_config, validate_config_values, start_run, reconcile_jobs
from sidock_ubuntu.job_store import create_job, load_job, update_job, project_lock, process_identity

def _lock_worker(path, ready, hold):
    from sidock_ubuntu.job_store import project_lock
    with project_lock(Path(path),'child'):
        ready.set(); hold.wait()
from sidock_ubuntu.state import update_compound, update_state, ensure_compounds

def fixture_inputs(tmp_path):
    pdb=tmp_path/'protein.pdb'; pdb.write_text('ATOM\n')
    sdf=tmp_path/'compounds.sdf'; sdf.write_text('minimal\n')
    return pdb,sdf

def test_project_name_and_workspace_are_confined(tmp_path):
    assert safe_project_name('My Project')=='My-Project'
    with pytest.raises(ValueError): safe_project_name('../escape')
    assert workspace_path(tmp_path,'ok').parent==tmp_path.resolve()

def test_project_creation_persists_canonical_inputs_and_state(tmp_path):
    pdb,sdf=fixture_inputs(tmp_path); project=create_project(tmp_path/'workspace','demo',pdb,sdf)
    assert (project/'input/protein.pdb').exists() and (project/'input/compounds.sdf').exists()
    state=json.loads((project/'project_state.json').read_text())
    assert set(('INPUT','REFERENCE','FINALIZATION')).issubset(state['stages'])

def test_project_view_reopens_from_disk(tmp_path):
    pdb,sdf=fixture_inputs(tmp_path); project=create_project(tmp_path/'workspace','demo',pdb,sdf)
    reopened=project_view(project)
    assert reopened['project']==str(project.resolve()) and reopened['state']['stages']['REFERENCE']['status']=='PENDING'

def test_existing_project_cannot_be_overwritten(tmp_path):
    pdb,sdf=fixture_inputs(tmp_path); create_project(tmp_path/'workspace','demo',pdb,sdf)
    with pytest.raises(FileExistsError): create_project(tmp_path/'workspace','demo',pdb,sdf)

def test_compound_results_are_reconstructed_from_persisted_state(tmp_path):
    pdb,sdf=fixture_inputs(tmp_path); project=create_project(tmp_path/'workspace','demo',pdb,sdf)
    ensure_compounds(project,['A']); update_compound(project,'A','preparation','PASS'); update_compound(project,'A','docking','FAIL',error_code='E502',error_message='AutoDock failed')
    rows=get_compounds(project); assert rows[0]['status']=='FAIL' and rows[0]['error_code']=='E502'
    detail=get_compound_detail(project,'A'); assert detail['stages']['docking']['error_code']=='E502'

def test_artifact_resolver_rejects_traversal_and_results_are_read_only(tmp_path):
    pdb,sdf=fixture_inputs(tmp_path); project=create_project(tmp_path/'workspace','demo',pdb,sdf)
    (project/'results.csv').write_text('compound_id,status\nA,BLOCKED\n')
    assert resolve_project_artifact(project,'results.csv')==(project/'results.csv').resolve()
    with pytest.raises(ValueError): resolve_project_artifact(project,'../protein.pdb')

def _write_evidence(project, status, compounds, validation=None):
    (project/'analysis').mkdir(exist_ok=True)
    (project/'analysis/evidence.json').write_text(json.dumps({'project_status':status,'validation':validation or {},'compounds':compounds}),encoding='utf-8')

def test_results_view_preserves_complete_and_partial_statuses(tmp_path):
    pdb,sdf=fixture_inputs(tmp_path); project=create_project(tmp_path/'workspace','demo',pdb,sdf); ensure_compounds(project,['A','B'])
    update_compound(project,'A','preparation','PASS'); update_compound(project,'B','preparation','FAIL',error_code='E502')
    _write_evidence(project,'PARTIAL',[{'compound_id':'A','lowest_docking_energy':-7.2},{'compound_id':'B','lowest_docking_energy':None}])
    view=get_project_results(project); assert view['project_status']=='PARTIAL' and view['compound_counts']['FAIL']==1

def test_results_view_preserves_blocked_validation_context(tmp_path):
    pdb,sdf=fixture_inputs(tmp_path); project=create_project(tmp_path/'workspace','demo',pdb,sdf); ensure_compounds(project,['A'])
    update_state(project,'VALIDATION','PASS')
    _write_evidence(project,'BLOCKED',[{'compound_id':'A','status':'BLOCKED','lowest_docking_energy':None}],{'sampling_result':'FAIL','native_redocking_RMSD_A':2.5,'rmsd_threshold':2.0})
    view=get_project_results(project); assert view['project_status']=='BLOCKED' and view['validation']['native_redocking_RMSD_A']==2.5

def test_zero_interaction_and_unavailable_values_remain_distinct(tmp_path):
    pdb,sdf=fixture_inputs(tmp_path); project=create_project(tmp_path/'workspace','demo',pdb,sdf); ensure_compounds(project,['A','B'])
    update_compound(project,'A','preparation','PASS'); update_compound(project,'B','preparation','PASS')
    _write_evidence(project,'PARTIAL',[{'compound_id':'A','PLIP_H_bonds':0},{'compound_id':'B','PLIP_H_bonds':None}])
    rows={x['compound_id']:x for x in get_compounds(project)}; evidence=json.loads((project/'analysis/evidence.json').read_text())['compounds']
    assert evidence[0]['PLIP_H_bonds']==0 and evidence[1]['PLIP_H_bonds'] is None

def test_doctor_wrapper_returns_structured_dependency_status():
    data=get_doctor_status(); assert data['dependencies'] and {'name','status','detail','required'} <= set(data['dependencies'][0])

def test_config_wrapper_exposes_scientific_defaults():
    data=get_config(); assert data['validation']['rmsd_threshold']==2.0 and data['docking']['population_size']==150

def test_config_validation_accepts_valid_and_rejects_invalid_values():
    assert validate_config_values({'validation':{'rmsd_threshold':2.0}})['valid']
    result=validate_config_values({'validation':{'rmsd_threshold':-1}}); assert result['valid'] is False and result['errors']

def test_project_resolved_config_is_read_from_disk(tmp_path):
    pdb,sdf=fixture_inputs(tmp_path); project=create_project(tmp_path/'workspace','demo',pdb,sdf)
    (project/'provenance').mkdir(exist_ok=True); (project/'provenance/resolved_config.yaml').write_text('validation:\n  rmsd_threshold: 1.75\n')
    data=get_project_config(project); assert data['config']['validation']['rmsd_threshold']==1.75 and data['source']=='project resolved configuration'

def test_gui_routes_include_doctor_and_config(tmp_path):
    pytest.importorskip('fastapi')
    from sidock_ubuntu.gui import create_app
    routes={getattr(route,'path',None) for route in create_app(tmp_path/'workspace').routes}
    assert {'/doctor','/api/system/doctor','/config','/api/config','/projects/{project_id}/config'} <= routes

def test_job_record_is_atomic_and_has_history(tmp_path):
    project=tmp_path/'project'; project.mkdir(); (project/'provenance').mkdir()
    job=create_job(project); update_job(project,status='FAILED',last_error={'message':'test'})
    assert load_job(project)['status']=='FAILED' and len((project/'provenance/gui_job_history.jsonl').read_text().splitlines())>=2

def test_live_project_lock_rejects_second_owner(tmp_path):
    project=tmp_path/'project'; project.mkdir()
    with project_lock(project,'one'):
        with pytest.raises(BlockingIOError):
            with project_lock(project,'two'): pass

def test_dead_running_job_reconciles_without_science(tmp_path):
    project=tmp_path/'project'; project.mkdir(); (project/'provenance').mkdir(); (project/'project_state.json').write_text(json.dumps({'stages':{'INPUT':{'status':'PENDING'}}}))
    job=create_job(project); update_job(project,status='RUNNING',pid=999999,process_start_identity='missing')
    out=reconcile_jobs(tmp_path); assert out[0]['status']=='INTERRUPTED'

def test_job_run_returns_before_backend_finishes(tmp_path, monkeypatch):
    pdb,sdf=fixture_inputs(tmp_path); project=create_project(tmp_path/'workspace','demo',pdb,sdf)
    import sidock_ubuntu.gui_service as service
    def backend(*args,**kwargs):
        import time; time.sleep(.15)
    monkeypatch.setattr(service,'run_stage_runner',backend)
    result=start_run(project); assert result['status']=='STARTING'
    assert load_job(project)['status'] in {'QUEUED','STARTING','RUNNING','SUCCEEDED'}
    import time; time.sleep(.3); assert load_job(project)['status']=='SUCCEEDED'

def test_duplicate_run_is_rejected(tmp_path, monkeypatch):
    pdb,sdf=fixture_inputs(tmp_path); project=create_project(tmp_path/'workspace','demo',pdb,sdf)
    import sidock_ubuntu.gui_service as service, time
    monkeypatch.setattr(service,'run_stage_runner',lambda *a,**k: time.sleep(.3))
    start_run(project)
    with pytest.raises(RuntimeError,match='PROJECT_ALREADY_RUNNING'): start_run(project)

def test_real_cross_process_lock_rejects_and_then_recovers(tmp_path):
    project=tmp_path/'project'; project.mkdir(); ready=multiprocessing.Event(); hold=multiprocessing.Event()
    child=multiprocessing.Process(target=_lock_worker,args=(str(project),ready,hold)); child.start(); assert ready.wait(3)
    with pytest.raises(BlockingIOError):
        with project_lock(project,'parent'): pass
    hold.set(); child.join(3); assert child.exitcode==0
    with project_lock(project,'parent'): pass

def test_reconciliation_is_idempotent_and_corrupt_job_is_browseable(tmp_path):
    project=tmp_path/'project'; project.mkdir(); (project/'provenance').mkdir(); (project/'project_state.json').write_text(json.dumps({'stages':{'INPUT':{'status':'PENDING'}}}))
    (project/'provenance/gui_job.json').write_text('{broken')
    assert load_job(project) is None
    first=reconcile_jobs(tmp_path); second=reconcile_jobs(tmp_path); assert first==second==[]

def test_gui_app_routes_when_fastapi_is_installed(tmp_path):
    pytest.importorskip('fastapi')
    from sidock_ubuntu.gui import create_app
    app=create_app(tmp_path/'workspace')
    routes={getattr(route,'path',None) for route in app.routes}
    assert '/projects' in routes and '/projects/{project_id}/reference' in routes and '/projects/{project_id}/resume' in routes

def test_multipart_minimal_project_creation_redirects_and_preserves_files(tmp_path):
    pytest.importorskip('fastapi')
    from fastapi.testclient import TestClient
    from sidock_ubuntu.gui import create_app
    client=TestClient(create_app(tmp_path/'workspace'))
    sdf=b'\n     RDKit          2D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5981   -0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  END\n$$$$\n'
    response=client.post('/projects',files={'protein':('native.pdb',b'ATOM\n','chemical/x-pdb'),'compounds':('library.sdf',sdf,'chemical/octet-stream')},data={'name':'Tryptase Screening'},follow_redirects=False)
    assert response.status_code==303 and response.headers['location'].endswith('/Tryptase-Screening')
    project=tmp_path/'workspace'/'Tryptase-Screening'
    assert (project/'input/protein.pdb').read_bytes()==b'ATOM\n'
    assert (project/'input/compounds.sdf').read_bytes()==sdf

def test_multipart_optional_sources_and_duplicate_protection(tmp_path):
    pytest.importorskip('fastapi')
    from fastapi.testclient import TestClient
    from sidock_ubuntu.gui import create_app
    client=TestClient(create_app(tmp_path/'workspace'))
    sdf=b'\n     RDKit          2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n'
    files={'protein':('native.pdb',b'ATOM','text/plain'),'compounds':('library.sdf',sdf,'text/plain'),'mmcif':('structure.cif',b'data_x','text/plain'),'ccd':('ccd.json',b'{}','application/json')}
    assert client.post('/projects',files=files,data={'name':'full project'},follow_redirects=False).status_code==303
    duplicate=client.post('/projects',files=files,data={'name':'full project'},follow_redirects=False)
    assert duplicate.status_code==200 and 'already exists' in duplicate.text

def test_compound_page_filters_and_paginates_deterministically(tmp_path):
    project=tmp_path/'project'; project.mkdir(); ensure_compounds(project,[f'C{i:02d}' for i in range(30)])
    from sidock_ubuntu.gui_service import get_compound_page
    first=get_compound_page(project,page=1,page_size=25); second=get_compound_page(project,page=2,page_size=25)
    assert first['total']==30 and len(first['items'])==25 and [x['compound_id'] for x in second['items']]==['C25','C26','C27','C28','C29']

def test_multipart_oversize_and_invalid_sdf_are_rejected_without_project(tmp_path, monkeypatch):
    pytest.importorskip('fastapi')
    from fastapi.testclient import TestClient
    import sidock_ubuntu.gui as gui
    from sidock_ubuntu.gui import create_app
    monkeypatch.setattr(gui,'MAX_UPLOAD_BYTES',5)
    client=TestClient(create_app(tmp_path/'workspace'))
    oversized=client.post('/projects',files={'protein':('x.pdb',b'ATOM','text/plain'),'compounds':('x.sdf',b'123456','text/plain')},data={'name':'too-big'},follow_redirects=False)
    assert oversized.status_code==200 and 'exceeds' in oversized.text and not (tmp_path/'workspace'/'too-big').exists()
    monkeypatch.setattr(gui,'MAX_UPLOAD_BYTES',100000)
    invalid=client.post('/projects',files={'protein':('x.pdb',b'ATOM','text/plain'),'compounds':('x.sdf',b'not an sdf','text/plain')},data={'name':'bad-sdf'},follow_redirects=False)
    assert invalid.status_code==200 and 'E002' in invalid.text and not (tmp_path/'workspace'/'bad-sdf').exists()

def test_dashboard_action_matrix_is_centralized(tmp_path, monkeypatch):
    import sidock_ubuntu.gui_service as service
    project=tmp_path/'project'; project.mkdir(); (project/'analysis').mkdir(); (project/'project_state.json').write_text('{}')
    cases=[('PENDING',None,'RUN'),('PENDING',{'status':'INTERRUPTED'},'RESUME'),('STALE',None,'RESUME'),('COMPLETE',None,'VIEW_RESULTS'),('COMPLETE_WITH_WARNINGS',None,'VIEW_RESULTS'),('PARTIAL',None,'VIEW_RESULTS'),('BLOCKED',None,'VIEW_VALIDATION'),('PENDING',{'status':'FAILED'},'RESUME')]
    for status,job,expected in cases:
        monkeypatch.setattr(service,'project_status',lambda p,s=status:{'stages':{'FINALIZATION':{'status':s},'VALIDATION':{'status':'PASS'}}})
        monkeypatch.setattr(service,'load_job',lambda p,j=job:j)
        monkeypatch.setattr(service,'get_compounds',lambda p:[])
        monkeypatch.setattr(service,'reference_view',lambda p:{'selection':{'selected_instance':'A','confidence':'HIGH'}})
        monkeypatch.setattr(service,'_evidence',lambda p:{})
        monkeypatch.setattr(service,'project_execution_lock',lambda *a,**k: __import__('contextlib').nullcontext())
        assert service.dashboard_view(project)['primary_action']==expected
    monkeypatch.setattr(service,'reference_view',lambda p:{'selection':{'confidence':'LOW'}})
    monkeypatch.setattr(service,'load_job',lambda p:None)
    assert service.dashboard_view(project)['primary_action']=='SELECT_REFERENCE'
