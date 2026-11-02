"""FastAPI thin client for persisted SIDock projects."""
from __future__ import annotations
import os
from pathlib import Path
from .gui_service import (create_project, get_compound_detail, get_compounds, get_compound_page, get_project_results,
                           get_config, get_doctor_status, get_job, get_job_history, get_project_config, get_worker_log, resolve_project_artifact, validate_config_values,
                           get_job, project_view, reconcile_jobs, resume_project, safe_project_name,
                           set_reference_override, start_run, workspace_path, validate_uploaded_inputs, MAX_UPLOAD_BYTES)
from .state import project_status, derived_project_status

try:
    from fastapi import FastAPI, File, Form, HTTPException, UploadFile, Query
    from fastapi.staticfiles import StaticFiles
    from fastapi.responses import HTMLResponse, FileResponse
except ImportError:  # pragma: no cover
    FastAPI = None

def _workspace(value=None):
    root=Path(value or os.environ.get('SIDOCK_WORKSPACE','sidock_projects')).resolve(); root.mkdir(parents=True,exist_ok=True); return root

def create_app(workspace=None):
    if FastAPI is None: raise RuntimeError('FastAPI is required for the v0.5 GUI')
    root=_workspace(workspace); reconcile_jobs(root); app=FastAPI(title='SIDock Ubuntu',version='0.5.1')
    static_dir=Path(__file__).with_name('static'); static_dir.mkdir(exist_ok=True); app.mount('/static', StaticFiles(directory=str(static_dir)), name='static')
    def page(title, body):
        import html as h
        return HTMLResponse(f"<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>{h.escape(title)} - SIDock</title><style>body{{font:16px system-ui;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#17202a}}nav{{display:flex;gap:1rem;margin-bottom:2rem;flex-wrap:wrap}}table{{border-collapse:collapse;width:100%;margin:1rem 0}}th,td{{padding:.55rem;border-bottom:1px solid #ccd;text-align:left}}.badge{{display:inline-block;padding:.2rem .5rem;border-radius:.35rem;background:#e8eef5}}.error{{color:#a00;background:#fee;padding:.7rem}}label{{display:block;margin:.8rem 0}}input,select,button{{font:inherit;padding:.35rem}}@media(max-width:650px){{table{{font-size:.85rem}}th,td{{padding:.3rem}}}}</style></head><body><nav><a href='/projects'>Projects</a><a href='/projects/new'>New Project</a><a href='/doctor'>Doctor</a><a href='/config'>Configuration</a></nav>{body}</body></html>")
    @app.get('/',response_class=HTMLResponse)
    def index(): return page('Projects', "<h1>SIDock Ubuntu</h1><p>Local thin client for persisted projects.</p><p><a href='/projects'>Open projects</a> or <a href='/projects/new'>create a project</a>.</p>")
    @app.get('/projects/new',response_class=HTMLResponse)
    def new_project():
        return page('New Project', "<h1>New project</h1><form action='/projects' method='post' enctype='multipart/form-data'><label>Project name <input name='name' required></label><label>Protein PDB <input name='protein' type='file' accept='.pdb,.ent' required></label><small>Required crystallographic receptor structure.</small><label>Compound SDF <input name='compounds' type='file' accept='.sdf' required></label><small>Required compound library.</small><label>mmCIF (optional) <input name='mmcif' type='file' accept='.cif,.mmcif'></label><label>Local CCD (optional) <input name='ccd' type='file' accept='.json,.cif'></label><button>Create project</button></form>")
    @app.get('/projects',response_class=HTMLResponse)
    def projects():
        import html as h
        items=[]
        for p in sorted((x for x in root.iterdir() if (x/'project_state.json').exists()), key=lambda x:x.stat().st_mtime, reverse=True):
            state=project_status(p); final=derived_project_status(state)
            items.append(f"<tr><td><a href='/projects/{h.escape(p.name)}'>{h.escape(p.name)}</a></td><td><span class='badge'>{h.escape(final)}</span></td><td><a href='/projects/{h.escape(p.name)}'>Open</a></td></tr>")
        body="<h1>Projects</h1>"+("<table><tr><th>Project</th><th>Status</th><th></th></tr>"+''.join(items)+"</table>" if items else "<p>No projects yet.</p><p><a href='/projects/new'>Create your first project</a></p>")
        return page('Projects',body)

    @app.get('/doctor',response_class=HTMLResponse)
    def doctor():
        data=get_doctor_status(); import html as _html
        rows=''.join(f"<tr><td>{_html.escape(x['name'])}</td><td>{_html.escape(x['status'])}</td><td>{_html.escape(str(x['detail']))}</td></tr>" for x in data['dependencies'])
        return HTMLResponse(f"<h1>System diagnostics</h1><p>Overall: {_html.escape(data['status'])}</p><table>{rows}</table><a href='/doctor'>Refresh</a>")

    @app.get('/api/system/doctor')
    def doctor_api(): return get_doctor_status()

    @app.get('/config',response_class=HTMLResponse)
    def config_page():
        import html as _html, json as _json
        data=get_config(); return HTMLResponse(f"<h1>Configuration</h1><pre>{_html.escape(_json.dumps(data,indent=2))}</pre><p>StageRunner execution remains serial; jobs above 1 are not enabled.</p>")

    @app.get('/api/config')
    def config_api(): return get_config()

    @app.post('/config/validate')
    async def config_validate(values: dict): return validate_config_values(values)

    @app.post('/api/config/validate')
    async def config_validate_api(values: dict): return validate_config_values(values)
    @app.post('/projects')
    async def create(name:str=Form(...),protein:UploadFile=File(...),compounds:UploadFile=File(...),mmcif:UploadFile|None=File(None),ccd:UploadFile|None=File(None)):
        upload=root/'.uploads'; upload.mkdir(exist_ok=True); token=safe_project_name(name); files=[]
        try:
            for item,suffix in ((protein,'.pdb'),(compounds,'.sdf'),(mmcif,'.cif'),(ccd,'.json')):
                if item is None: continue
                allowed={'.pdb':('.pdb','.ent'),'.sdf':('.sdf',),'.cif':('.cif','.mmcif'),'.json':('.json','.cif')}[suffix]
                if Path(item.filename or '').suffix.lower() not in allowed: raise ValueError('E001 INVALID_INPUT_FILE_EXTENSION')
                path=upload/f'{token}-{len(files)}{suffix}'; total=0
                with path.open('wb') as fh:
                    while True:
                        chunk=await item.read(1024*1024)
                        if not chunk: break
                        total += len(chunk)
                        if total > MAX_UPLOAD_BYTES: raise ValueError('E001 UPLOAD_TOO_LARGE: Uploaded file exceeds the configured size limit.')
                        fh.write(chunk)
                files.append(path)
            mmcif_file=files[2] if mmcif else None
            ccd_file=files[3] if (mmcif and ccd) else (files[2] if ccd else None)
            validate_uploaded_inputs(files[0],files[1],mmcif_file,ccd_file)
            project=create_project(root,name,files[0],files[1],mmcif_file,ccd_file); from fastapi.responses import RedirectResponse; return RedirectResponse(f'/projects/{project.name}',status_code=303)
        except (ValueError,FileExistsError) as exc:
            import html as h
            label='A project with this identifier already exists.' if str(exc).startswith('E003') else str(exc)
            return page('New Project',f"<h1>New project</h1><p class='error'>{h.escape(label)}</p><p><a href='/projects/new'>Return to the project form</a></p>")
    def get_project(project_id):
        try: project=workspace_path(root,project_id)
        except ValueError as exc: raise HTTPException(status_code=400,detail=str(exc))
        if not (project/'project_state.json').exists(): raise HTTPException(status_code=404,detail='E001 PROJECT_NOT_FOUND')
        return project
    @app.get('/projects/{project_id}')
    def detail(project_id):
        data=project_view(get_project(project_id)); import html as h
        state=data['state']; final=derived_project_status(state); counts=data.get('reference',{}); dash=data['dashboard']
        body=f"<h1>{h.escape(Path(data['project']).name)}</h1><p>Scientific status: <span class='badge'>{h.escape(final)}</span></p><p>Operational job: <span class='badge'>{h.escape(str(dash.get('operational_job_status') or 'No GUI job recorded'))}</span></p><p>Reference: {h.escape(str(dash['reference_status']))}; validation: {h.escape(str(dash['validation_status']))}</p><p>Compounds: {dash['compound_total']} total, {dash['compound_success']} successful, {dash['compound_failed']} failed, {dash['compound_blocked']} blocked.</p><p>{h.escape(dash['explanatory_message'])}</p><p><a href='{h.escape(dash['primary_action_url'] or '#')}'>{h.escape(dash['primary_action'])}</a></p><p>Reference artifacts: {'available' if counts else 'not selected yet'}</p><nav><a href='/projects/{project_id}/reference'>Reference</a><a href='/projects/{project_id}/compounds'>Compounds</a><a href='/projects/{project_id}/results'>Results</a><a href='/projects/{project_id}/jobs'>Jobs</a><a href='/projects/{project_id}/config'>Config</a></nav>"
        if final=='BLOCKED': body += '<p class="error">Native redocking validation failed; compound docking was not performed.</p>'
        if data.get('active_job'): body += '<p>Project currently executing.</p>'
        return page(Path(data['project']).name,body)
    @app.get('/projects/{project_id}/status')
    def status(project_id): return project_status(get_project(project_id))

    @app.get('/api/projects/{project_id}/job')
    def job(project_id): return get_job(get_project(project_id))

    @app.get('/projects/{project_id}/jobs',response_class=HTMLResponse)
    def jobs_page(project_id):
        import html as _html
        history=get_job_history(get_project(project_id)); rows=''.join(f"<tr><td>{_html.escape(str(x.get('attempt')))}</td><td>{_html.escape(str(x.get('action')))}</td><td>{_html.escape(str(x.get('status')))}</td><td>{_html.escape(str(x.get('last_error') or ''))}</td></tr>" for x in history)
        body='<h1>Job history</h1>'+ (f"<table><tr><th>Attempt</th><th>Action</th><th>Status</th><th>Error</th></tr>{rows}</table>" if rows else '<p>No GUI execution jobs have been recorded for this project.</p>')
        return page('Jobs - '+project_id,body)

    @app.get('/api/projects/{project_id}/jobs')
    def jobs_api(project_id): return get_job_history(get_project(project_id))

    @app.get('/projects/{project_id}/worker-log',response_class=HTMLResponse)
    def worker_log(project_id):
        import html as _html
        return HTMLResponse(f"<h1>Worker log</h1><pre>{_html.escape(get_worker_log(get_project(project_id)))}</pre>")

    @app.get('/projects/{project_id}/config',response_class=HTMLResponse)
    def project_config(project_id):
        import html as _html, json as _json
        data=get_project_config(get_project(project_id)); return HTMLResponse(f"<h1>Project configuration</h1><p>{_html.escape(str(data['source']))}</p><pre>{_html.escape(_json.dumps(data.get('config'),indent=2,default=str))}</pre>")

    @app.get('/api/projects/{project_id}/config')
    def project_config_api(project_id): return get_project_config(get_project(project_id))
    @app.get('/projects/{project_id}/reference')
    def reference(project_id): return project_view(get_project(project_id))['reference']
    @app.post('/projects/{project_id}/reference')
    def choose_reference(project_id, instance: str):
        try: return set_reference_override(get_project(project_id), instance)
        except ValueError as exc: raise HTTPException(status_code=400, detail=str(exc))
    @app.post('/projects/{project_id}/run')
    def run(project_id):
        try: return start_run(get_project(project_id))
        except RuntimeError as exc: raise HTTPException(status_code=409,detail=str(exc))
    @app.post('/projects/{project_id}/resume')
    def resume(project_id):
        try: return resume_project(get_project(project_id))
        except RuntimeError as exc: raise HTTPException(status_code=409,detail=str(exc))
    @app.get('/projects/{project_id}/compounds',response_class=HTMLResponse)
    def compounds(project_id, page_num:int=Query(1, alias='page', ge=1), page_size:int=25, status:str='ALL', search:str='', sort:str='compound_id'):
        data=get_compound_page(get_project(project_id),page=page_num,page_size=page_size,status=status,search=search,sort=sort); body=[f"<h1>Compounds</h1><form><label>Search <input name='search' value='{__import__('html').escape(search)}'></label><label>Status <select name='status'><option>ALL</option><option>PASS</option><option>FAIL</option><option>BLOCKED</option><option>STALE</option></select></label><input type='hidden' name='page_size' value='{data['page_size']}'><button>Filter</button></form><p>Showing {min((data['page']-1)*data['page_size']+1,data['total']) if data['total'] else 0}–{min(data['page']*data['page_size'],data['total'])} of {data['total']}</p><table><tr><th>ID</th><th>Status</th><th>Docking</th><th>Details</th></tr>"]
        for row in data['items']:
            cid=__import__('html').escape(str(row['compound_id'])); body.append(f"<tr><td>{cid}</td><td>{__import__('html').escape(row['status'])}</td><td>{__import__('html').escape(str(row.get('docking_energy') if row.get('docking_energy') is not None else 'Unavailable'))}</td><td><a href='/projects/{project_id}/compounds/{cid}'>Details</a></td></tr>")
        if not data['total']:
            return page('Compounds - '+project_id,'<h1>Compounds</h1><p>No compounds are available in this project.</p>')
        links=[]
        for n in range(1,data['pages']+1): links.append(f"<a href='?page={n}&page_size={data['page_size']}&status={status}&search={__import__('urllib.parse').parse.quote(search)}&sort={sort}'>{n}</a>")
        return page('Compounds',''.join(body)+f"</table><p>{' '.join(links)}</p>")

    @app.get('/api/projects/{project_id}/compounds')
    def compounds_api(project_id): return get_compounds(get_project(project_id))

    @app.get('/projects/{project_id}/compounds/{compound_id}')
    def compound_detail(project_id, compound_id):
        try: detail=get_compound_detail(get_project(project_id),compound_id)
        except KeyError: raise HTTPException(status_code=404,detail='COMPOUND_NOT_FOUND')
        import html as _html
        status=_html.escape(detail['status']); evidence=detail.get('evidence',{})
        def val(key): return evidence.get(key) if evidence.get(key) is not None else 'Unavailable'
        return HTMLResponse(f"<h1>{_html.escape(compound_id)}</h1><p>Status: <b>{status}</b></p><p>Docking energy (predicted): {_html.escape(str(val('lowest_docking_energy')))}</p><p>Pose: {_html.escape(str(val('selected_pose')))}</p><p>H-bonds: {_html.escape(str(val('PLIP_H_bonds') if evidence.get('PLIP_H_bonds') is not None else 'Unavailable'))}</p><p>Error: {_html.escape(str(detail.get('error_code') or 'None'))}</p>")

    @app.get('/api/projects/{project_id}/compounds/{compound_id}')
    def compound_detail_api(project_id, compound_id):
        try: return get_compound_detail(get_project(project_id),compound_id)
        except KeyError: raise HTTPException(status_code=404,detail='COMPOUND_NOT_FOUND')

    @app.get('/projects/{project_id}/results', response_class=HTMLResponse)
    def results_page(project_id):
        data=get_project_results(get_project(project_id)); import html as _html
        validation=data['validation']; counts=data['compound_counts']; blocked='';
        if data['project_status']=='BLOCKED': blocked='<p><b>Native redocking validation failed; compound docking was not performed.</b></p>'
        links=[]
        for rel,label in (('results.csv','results.csv'),('analysis/evidence.json','evidence.json'),('report.html','report.html')):
            if (get_project(project_id)/rel).is_file(): links.append(f"<a href='/projects/{project_id}/artifacts/{rel}'>{label}</a>")
        empty='' if links else '<p>Final results are not available yet.</p>'
        return page('Results - '+project_id,f"<h1>Results</h1><p>Project status: <b>{_html.escape(str(data['project_status']))}</b></p>{blocked}<p>Validation: {_html.escape(str(validation.get('sampling_result','Not run yet')))}; RMSD {_html.escape(str(validation.get('native_redocking_RMSD_A','Unavailable')))}; threshold {_html.escape(str(validation.get('rmsd_threshold','Unavailable')))}</p><p>Counts: {_html.escape(str(counts))}</p>{empty}{' | '.join(links)}")

    @app.get('/api/projects/{project_id}/results')
    def results_api(project_id): return get_project_results(get_project(project_id))

    @app.get('/projects/{project_id}/artifacts/{artifact_path:path}')
    def artifact(project_id, artifact_path):
        try: path=resolve_project_artifact(get_project(project_id),artifact_path)
        except (ValueError,FileNotFoundError) as exc: raise HTTPException(status_code=400,detail=str(exc))
        return FileResponse(path)
    @app.get('/projects/{project_id}/evidence')
    def evidence(project_id):
        path=get_project(project_id)/'analysis'/'evidence.json'
        if not path.exists(): raise HTTPException(status_code=404,detail='EVIDENCE_NOT_AVAILABLE')
        return FileResponse(path,media_type='application/json',filename='evidence.json')
    return app

app=create_app() if FastAPI is not None else None

def main():
    if FastAPI is None: raise SystemExit('FastAPI is required for the v0.5 GUI')
    import uvicorn; uvicorn.run('sidock_ubuntu.gui:app',host='127.0.0.1',port=8000,reload=False)
