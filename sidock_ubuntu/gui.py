"""FastAPI thin client for persisted SIDock projects."""
from __future__ import annotations
import os
from pathlib import Path
from .gui_service import (create_project, get_compound_detail, get_compounds, get_compound_page, get_project_results,
                           get_config, get_doctor_status, get_job, get_job_history, get_project_config, get_worker_log, resolve_project_artifact, validate_config_values,
                           get_job, project_view, reconcile_jobs, resume_project, safe_project_name, input_compound_count,
                           set_reference_override, start_run, workspace_path, validate_uploaded_inputs, attach_project_inputs, MAX_UPLOAD_BYTES)
from .state import project_status, derived_project_status

try:
    from fastapi import FastAPI, File, Form, HTTPException, UploadFile, Query, Request
    from fastapi.staticfiles import StaticFiles
    from fastapi.responses import HTMLResponse, FileResponse
except ImportError:  # pragma: no cover
    FastAPI = None

def _workspace(value=None):
    root=Path(value or os.environ.get('SIDOCK_WORKSPACE','sidock_projects')).resolve(); root.mkdir(parents=True,exist_ok=True); return root

def create_app(workspace=None):
    if FastAPI is None: raise RuntimeError('FastAPI is required for the v0.5 GUI')
    root=_workspace(workspace); reconcile_jobs(root); app=FastAPI(title='SIDock Ubuntu',version='0.6.0.dev0')
    static_dir=Path(__file__).with_name('static'); static_dir.mkdir(exist_ok=True); app.mount('/static', StaticFiles(directory=str(static_dir)), name='static')
    def page(title, body):
        import html as h
        return HTMLResponse(f"<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>{h.escape(title)} - SIDock</title><style>body{{font:16px system-ui;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#17202a}}nav{{display:flex;gap:1rem;margin-bottom:2rem;flex-wrap:wrap}}table{{border-collapse:collapse;width:100%;margin:1rem 0}}th,td{{padding:.55rem;border-bottom:1px solid #ccd;text-align:left}}.badge{{display:inline-block;padding:.2rem .5rem;border-radius:.35rem;background:#e8eef5}}.error{{color:#a00;background:#fee;padding:.7rem}}label{{display:block;margin:.8rem 0}}input,select,button{{font:inherit;padding:.35rem}}@media(max-width:650px){{table{{font-size:.85rem}}th,td{{padding:.3rem}}}}</style></head><body><nav><a href='/projects'>Projects</a><a href='/projects/new'>New Project</a><a href='/doctor'>Doctor</a><a href='/config'>Configuration</a></nav>{body}</body></html>")
    @app.get('/',response_class=HTMLResponse)
    def index(): return page('Projects', "<h1>SIDock Ubuntu</h1><p>Local thin client for persisted projects.</p><p><a href='/projects'>Open projects</a> or <a href='/projects/new'>create a project</a>.</p>")
    @app.get('/projects/new', response_class=HTMLResponse)
    def guided_new_redirect():
        from fastapi.responses import RedirectResponse
        return RedirectResponse('/guided/new', status_code=303)

    @app.get('/guided/new', response_class=HTMLResponse)
    def guided_new_project():
        return page('Create project', "<h1>Step 1 — Create docking project</h1>"
                    "<p>Enter the project name. Protein selection and compound upload happen in their own review steps.</p>"
                    "<form action='/guided/projects' method='post'><label>Project name "
                    "<input name='name' required maxlength='80'></label>"
                    "<button>Create project and continue</button></form>"
                    "<p><a href='/projects/new/automatic'>Legacy fully automated creation (optional)</a></p>")

    @app.post('/guided/projects')
    def guided_create(name: str=Form(...)):
        from .guided_workflow import create_guided_project
        from fastapi.responses import RedirectResponse
        try:
            project = create_guided_project(root, name)
        except (ValueError, FileExistsError) as exc:
            raise HTTPException(status_code=409, detail=str(exc))
        return RedirectResponse(f'/projects/{project.name}/guided/1', status_code=303)

    @app.get('/projects/new/automatic',response_class=HTMLResponse)
    def new_project():
        return page('New Project', "<h1>New project</h1><form action='/projects' method='post' enctype='multipart/form-data'><label>Project name <input name='name' required></label><label>Protein PDB <input name='protein' type='file' accept='.pdb,.ent' required></label><small>Required crystallographic or prepared receptor structure.</small><label>Compound SDF <input name='compounds' type='file' accept='.sdf' required></label><small>Required compound library.</small><label>mmCIF (optional) <input name='mmcif' type='file' accept='.cif,.mmcif'></label><label>Local CCD (optional) <input name='ccd' type='file' accept='.json,.cif'></label><label>Validation structure (optional) <input name='validation_structure' type='file' accept='.pdb,.ent,.cif,.mmcif'></label><label>PDB ID (optional) <input name='pdb_id'></label><button>Create project</button></form>")
    @app.get('/projects',response_class=HTMLResponse)
    def projects():
        import html as h
        items=[]
        for p in sorted((x for x in root.iterdir() if (x/'project_state.json').exists()), key=lambda x:x.stat().st_mtime, reverse=True):
            state=project_status(p); final=derived_project_status(state)
            guided=(p/'manual_docking/workflow.json').is_file()
            if guided:
                from .guided_workflow import approved
                if approved(p,4): final='READY FOR DOCKING'
                elif approved(p,3): final='STEP 4: SETUP REVIEW'
                elif approved(p,2): final='STEP 3: COMPOUND REVIEW'
                else: final='STEP 2: PROTEIN REVIEW'
            action=(f'/projects/{h.escape(p.name)}/guided/1' if guided
                    else f'/projects/{h.escape(p.name)}')
            items.append(f"<tr><td><a href='{action}'>{h.escape(p.name)}</a></td><td><span class='badge'>{h.escape(final)}</span></td><td><a href='{action}'>Open</a></td></tr>")
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
        data=get_config(); return HTMLResponse(f"<h1>Configuration</h1><pre>{_html.escape(_json.dumps(data,indent=2))}</pre><p>AutoDock compound jobs can run concurrently; parallel_jobs=auto uses the CPUs available to the SIDock process.</p>")

    @app.get('/api/config')
    def config_api(): return get_config()

    @app.post('/config/validate')
    async def config_validate(values: dict): return validate_config_values(values)

    @app.post('/api/config/validate')
    async def config_validate_api(values: dict): return validate_config_values(values)
    @app.post('/projects')
    async def create(name:str=Form(...),protein:UploadFile=File(...),compounds:UploadFile=File(...),mmcif:UploadFile|None=File(None),ccd:UploadFile|None=File(None),validation_structure:UploadFile|None=File(None),pdb_id:str|None=Form(None)):
        upload=root/'.uploads'; upload.mkdir(exist_ok=True); token=safe_project_name(name); files=[]
        try:
            for item,suffix in ((protein,'.pdb'),(compounds,'.sdf'),(mmcif,'.cif'),(ccd,'.json'),(validation_structure,'.pdb')):
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
            validation_file=files[4] if validation_structure else None
            project=create_project(root,name,files[0],files[1],mmcif_file,ccd_file,validation_file,pdb_id); from fastapi.responses import RedirectResponse; return RedirectResponse(f'/projects/{project.name}',status_code=303)
        except (ValueError,FileExistsError) as exc:
            import html as h
            label='A project with this identifier already exists.' if str(exc).startswith('E003') else str(exc)
            return page('New Project',f"<h1>New project</h1><p class='error'>{h.escape(label)}</p><p><a href='/projects/new'>Return to the project form</a></p>")
    def get_project(project_id):
        try: project=workspace_path(root,project_id)
        except ValueError as exc: raise HTTPException(status_code=400,detail=str(exc))
        if not (project/'project_state.json').exists(): raise HTTPException(status_code=404,detail='E001 PROJECT_NOT_FOUND')
        return project

    @app.get('/projects/{project_id}/guided/{step}', response_class=HTMLResponse)
    def guided_step(project_id: str, step: int):
        from .guided_ui import render
        if step not in (1, 2, 3, 4, 5):
            raise HTTPException(status_code=404, detail='UNKNOWN_GUIDED_STEP')
        project = get_project(project_id)
        return page(f'Step {step} - {project.name}', render(project, step))

    @app.post('/projects/{project_id}/guided/protein-upload')
    async def guided_protein_upload(project_id: str, protein: UploadFile=File(...)):
        from .guided_workflow import save_input
        from fastapi.responses import RedirectResponse
        project = get_project(project_id)
        if Path(protein.filename or '').suffix.lower() not in ('.pdb', '.ent'):
            raise HTTPException(status_code=400, detail='Expected a .pdb or .ent file')
        target = project / 'input' / 'guided_protein_upload.ent'
        target.parent.mkdir(parents=True, exist_ok=True)
        total = 0
        with target.open('wb') as handle:
            while True:
                chunk = await protein.read(1024*1024)
                if not chunk: break
                total += len(chunk)
                if total > MAX_UPLOAD_BYTES:
                    target.unlink(missing_ok=True)
                    raise HTTPException(status_code=413, detail='Protein file is too large')
                handle.write(chunk)
        try:
            save_input(project, target, 'protein')
        except (RuntimeError, ValueError) as exc:
            raise HTTPException(status_code=400, detail=str(exc))
        finally:
            target.unlink(missing_ok=True)
        return RedirectResponse(f'/projects/{project_id}/guided/2',status_code=303)

    @app.post('/projects/{project_id}/guided/protein-review')
    def guided_protein_review(project_id:str, ligand_key:str=Form(''),
                              remove_groups:list[str]=Form([]),
                              remove_waters:str|None=Form(None),
                              padding:float=Form(8.0), center_x:str=Form(''),
                              center_y:str=Form(''), center_z:str=Form(''),
                              exploratory_metal_free:str|None=Form(None),
                              acknowledge_non_ligand_component:str|None=Form(None)):
        from .guided_workflow import select_protein
        from .manual_docking import settings as manual_settings, save_settings
        from fastapi.responses import RedirectResponse
        project=get_project(project_id)
        try:
            center=None if ligand_key else (float(center_x),float(center_y),float(center_z))
            select_protein(project, ligand_key=ligand_key or None,
                           remove_groups=remove_groups, remove_waters=bool(remove_waters),
                           center=center, padding=padding,
                           acknowledge_non_ligand_component=bool(acknowledge_non_ligand_component))
            values=manual_settings(project)
            values["allow_missing_metal"]=bool(exploratory_metal_free)
            save_settings(project, values)
        except (ValueError,RuntimeError) as exc:
            raise HTTPException(status_code=400,detail=str(exc))
        return RedirectResponse(f'/projects/{project_id}/guided/2',status_code=303)

    @app.post('/projects/{project_id}/guided/receptor-import')
    async def guided_receptor_import(project_id:str, upload: UploadFile|None=File(None),
                                     existing_path:str|None=Form(None)):
        from .manual_docking import root as manual_root, job_status
        from .guided_workflow import _json, install_reviewed_receptor
        from fastapi.responses import RedirectResponse
        project=get_project(project_id)
        if job_status(project).get('status') == 'RUNNING':
            raise HTTPException(status_code=409,detail='Manual task is running')
        if not _json(manual_root(project)/'inputs/protein_choice.json'):
            raise HTTPException(status_code=409,detail='Save Step 2 protein review first')
        try:
            if upload is not None and upload.filename:
                if Path(upload.filename).suffix.lower() != '.pdbqt':
                    raise ValueError('Only a prepared .pdbqt is accepted')
                destination = manual_root(project)/'uploaded_receptor.pdbqt'
                destination.parent.mkdir(parents=True,exist_ok=True)
                total=0
                with destination.open('wb') as out:
                    while True:
                        chunk=await upload.read(1024*1024)
                        if not chunk: break
                        total += len(chunk)
                        if total > MAX_UPLOAD_BYTES:
                            destination.unlink(missing_ok=True)
                            raise ValueError('Uploaded receptor exceeds size limit')
                        out.write(chunk)
                install_reviewed_receptor(project,destination)
            elif existing_path:
                path=(project/existing_path.strip()).resolve()
                if not path.is_relative_to(project.resolve()):
                    raise ValueError('PDBQT must be inside this project')
                install_reviewed_receptor(project,path)
            else:
                raise ValueError('Select a prepared PDBQT or provide an existing project-relative path')
        except (ValueError,RuntimeError,OSError) as exc:
            raise HTTPException(status_code=400,detail=str(exc))
        return RedirectResponse(f'/projects/{project_id}/guided/2',status_code=303)

    @app.post('/projects/{project_id}/guided/compound-upload')
    async def guided_compound_upload(project_id:str, compounds:UploadFile=File(...)):
        from .guided_workflow import save_input
        from fastapi.responses import RedirectResponse
        project=get_project(project_id)
        if Path(compounds.filename or '').suffix.lower() != '.sdf':
            raise HTTPException(status_code=400,detail='Expected a .sdf file')
        temp=project/'input/guided_compounds_upload.sdf'
        temp.parent.mkdir(parents=True,exist_ok=True)
        total=0
        with temp.open('wb') as out:
            while True:
                chunk=await compounds.read(1024*1024)
                if not chunk: break
                total+=len(chunk)
                if total>MAX_UPLOAD_BYTES:
                    temp.unlink(missing_ok=True)
                    raise HTTPException(status_code=413,detail='SDF file is too large')
                out.write(chunk)
        try:
            from rdkit import Chem
            if not any(x is not None for x in Chem.SDMolSupplier(str(temp),sanitize=True,removeHs=False)):
                raise ValueError('SDF contains no valid molecules')
            save_input(project,temp,'compounds')
        except (RuntimeError,ValueError) as exc:
            raise HTTPException(status_code=400,detail=str(exc))
        finally:
            temp.unlink(missing_ok=True)
        return RedirectResponse(f'/projects/{project_id}/guided/3',status_code=303)

    @app.post('/projects/{project_id}/guided/settings')
    def guided_save_settings(project_id: str, center_x:float=Form(...), center_y:float=Form(...),
                             center_z:float=Form(...), npts_x:int=Form(...), npts_y:int=Form(...),
                             npts_z:int=Form(...), spacing:float=Form(...), ga_runs:int=Form(...),
                             parallel_jobs:int=Form(...), ga_num_evals:int=Form(...),
                             ga_pop_size:int=Form(...), ga_num_generations:int=Form(...),
                             allow_missing_metal:str|None=Form(None)):
        from .manual_docking import save_settings
        from fastapi.responses import RedirectResponse
        project=get_project(project_id)
        try:
            save_settings(project,locals())
        except (ValueError,TypeError) as exc:
            raise HTTPException(status_code=400,detail=str(exc))
        return RedirectResponse(f'/projects/{project_id}/guided/4',status_code=303)

    @app.post('/projects/{project_id}/guided/parameters')
    async def guided_custom_parameters(project_id: str, upload: UploadFile|None=File(None),
                                       existing_path: str|None=Form(None),
                                       next_step:int=Form(4)):
        from .manual_docking import install_parameter_file, root as manual_root, job_status
        from fastapi.responses import RedirectResponse
        project=get_project(project_id)
        if job_status(project).get('status') == 'RUNNING':
            raise HTTPException(status_code=409, detail='Manual task running; parameters cannot be replaced')
        try:
            if upload is not None and upload.filename:
                if Path(upload.filename).suffix.lower() not in ('.dat','.txt'):
                    raise ValueError('Expected an AutoDock4 .dat or .txt force-field library')
                path=manual_root(project)/'uploaded_custom_parameters.dat'
                path.parent.mkdir(parents=True,exist_ok=True)
                total=0
                with path.open('wb') as out:
                    while True:
                        chunk=await upload.read(1024*1024)
                        if not chunk: break
                        total+=len(chunk)
                        if total>2_000_000:
                            path.unlink(missing_ok=True)
                            raise ValueError('Parameter file exceeds 2 MB')
                        out.write(chunk)
                install_parameter_file(project,path)
            elif existing_path and existing_path.strip():
                source=(project/existing_path.strip()).resolve()
                if not source.is_relative_to(project.resolve()):
                    raise ValueError('Parameter file must be inside this project')
                install_parameter_file(project,source)
            else:
                raise ValueError('Select a complete AutoDock4 parameter library')
        except (ValueError,RuntimeError,OSError,UnicodeError) as exc:
            raise HTTPException(status_code=400,detail=str(exc))
        return RedirectResponse(f'/projects/{project_id}/guided/{2 if next_step == 2 else 4}',status_code=303)

    @app.post('/projects/{project_id}/guided/ni-reconstruct')
    def guided_reconstruct_nickel(project_id: str, atom_type:str=Form(...),
                                  site_charges:str=Form(...),
                                  parameter_source:str=Form(...),
                                  charge_source:str=Form(...),
                                  confirm_nickel_model:str|None=Form(None)):
        from .nickel_reconstruction import reconstruct, parse_charges
        from .manual_docking import job_status
        from fastapi.responses import RedirectResponse
        project=get_project(project_id)
        try:
            if job_status(project).get('status') == 'RUNNING':
                raise RuntimeError('A manual worker is still running; reconstruction is blocked')
            result=reconstruct(
                project, atom_type=atom_type.strip(),
                charges=parse_charges(site_charges),
                parameter_source=parameter_source,
                charge_source=charge_source,
                nonbonded_model_acknowledged=(confirm_nickel_model == 'I_ACKNOWLEDGE_NI_MODEL'),
            )
        except (RuntimeError, ValueError, OSError) as exc:
            raise HTTPException(status_code=409, detail=str(exc))
        return RedirectResponse(f'/projects/{project_id}/guided/2', status_code=303)

    @app.post('/projects/{project_id}/guided/approve')
    def guided_approve(project_id:str, step:int=Form(...)):
        from .guided_workflow import approve
        from fastapi.responses import RedirectResponse
        project=get_project(project_id)
        try:
            approve(project,step)
        except (ValueError,RuntimeError) as exc:
            raise HTTPException(status_code=409,detail=str(exc))
        return RedirectResponse(f'/projects/{project_id}/guided/{min(step+1,5)}',status_code=303)

    @app.post('/projects/{project_id}/guided/approve-exploratory')
    def guided_approve_exploratory(project_id: str,
                                   acknowledgement: str = Form(...)):
        from .guided_workflow import approve_exploratory
        from fastapi.responses import RedirectResponse
        project = get_project(project_id)
        try:
            approve_exploratory(project, acknowledgement=acknowledgement)
        except (RuntimeError, ValueError) as exc:
            raise HTTPException(status_code=409, detail=str(exc))
        return RedirectResponse(f'/projects/{project_id}/guided/3', status_code=303)

    @app.post('/projects/{project_id}/guided/action')
    def guided_run_action(project_id:str, action:str=Form(...),
                          confirm_metal_free:str|None=Form(None)):
        from .guided_workflow import approved, compound_profile, ready_to_dock, import_existing_compounds, _json
        from .manual_docking import launch, root as manual_root, job_status, settings as manual_settings, save_settings
        from fastapi.responses import RedirectResponse
        project=get_project(project_id)
        if action not in ('receptor','receptor_metal_free','receptor_ni_intermediate','receptor_openbabel','receptor_openbabel_pipeline','exploratory_smoke','prepare','reuse_compounds','grid','dock'):
            raise HTTPException(status_code=400,detail='Unknown step action')
        try:
            if action in ('receptor', 'receptor_metal_free', 'receptor_ni_intermediate', 'receptor_openbabel', 'receptor_openbabel_pipeline', 'exploratory_smoke'):
                if not _json(manual_root(project)/'inputs/protein_choice.json'):
                    raise RuntimeError('First save Step 2 protein and crystal ligand choices.')
                if action == 'receptor_metal_free':
                    if confirm_metal_free != 'I_ACKNOWLEDGE_METAL_FREE':
                        raise RuntimeError(
                            'Explicitly acknowledge metal-free exploratory receptor chemistry in Step 2.'
                        )
                    if job_status(project).get('status') == 'RUNNING':
                        raise RuntimeError('A manual preparation action is still active')
                    values=manual_settings(project)
                    values['allow_missing_metal']=True
                    save_settings(project,values)
                step=2
            elif action == 'prepare':
                # Ligand chemistry is independent of receptor preparation.
                # Only Step 3 approval and subsequent docking need approved Step 2.
                profile=compound_profile(project)
                if profile.get('errors') or not profile.get('available'):
                    raise RuntimeError('Compound SDF must contain valid, uniquely named molecules: '+ '; '.join(profile.get('errors',[])))
                step=3
            elif action == 'reuse_compounds':
                # Reuse is also independent of protein approval; verify the SDF
                # and filenames before copying in guided_workflow.
                if job_status(project).get('status') == 'RUNNING':
                    raise RuntimeError('Wait until the active manual task is finished')
                import_existing_compounds(project)
                return RedirectResponse(f'/projects/{project_id}/guided/3',status_code=303)
            elif action == 'grid':
                if not approved(project,2) or not approved(project,3):
                    raise RuntimeError('Approve both protein and compound preparation first')
                step=4
            else:
                ready_to_dock(project)
                step=5
            launch(project,action)
        except (RuntimeError,ValueError) as exc:
            raise HTTPException(status_code=409,detail=str(exc))
        return RedirectResponse(f'/projects/{project_id}/guided/{step}',status_code=303)

    @app.get('/projects/{project_id}')
    def detail(project_id):
        project=get_project(project_id)
        if (project/'manual_docking/workflow.json').is_file():
            from .guided_workflow import approved
            from .manual_docking import job_status
            import html as _h
            progress=next((n for n in (4,3,2) if approved(project,n)),1)
            goal=min(progress+1,5)
            current_job=job_status(project)
            run_message=(f"<p>Last manual job: {_h.escape(str(current_job.get('action')))} "
                         f"— {_h.escape(str(current_job.get('status')))}</p>"
                         if current_job else "")
            return page(project.name,
                f"<h1>{_h.escape(project.name)}</h1>"
                "<p>SIDock guided workflow. Every scientific stage is manually reviewed. "
                "The previous automatic pipeline state is independent of these approvals.</p>"
                f"<p><b>Approved through Step {progress}</b> (Step 1 means project created).</p>"
                + run_message
                + f"<p><a href='/projects/{_h.escape(project_id)}/guided/{goal}'>"
                  f"Continue to Step {goal}</a></p>"
                + f"<p><a href='/projects/{_h.escape(project_id)}/guided/1'>All five steps</a> "
                  f"| <a href='/projects/{_h.escape(project_id)}/results'>Legacy results</a></p>")
        data=project_view(project); import html as h
        state=data['state']; final=derived_project_status(state); counts=data.get('reference',{}); dash=data['dashboard']
        if dash['primary_action']=='RUN':
            policy=dash.get('validation_policy','required')
            options=''.join(f"<option value='{value}'{' selected' if value==policy else ''}>{label}</option>" for value,label in (('if_available','Validate if available (continue as unvalidated if RMSD fails)'),('required','Require native validation (RMSD ≤ 2 Å)'),('skip','Skip native redocking validation')))
            parallel=dash.get('parallel_jobs','auto'); detected=dash.get('detected_cpus',1); ga=dash.get('ga_runs',10)
            action=(f"<form method='post' action='{h.escape(dash['primary_action_url'])}'>"
                    f"<label>Validation policy <select name='validation_mode'>{options}</select></label>"
                    f"<label>GA runs per compound <input name='ga_runs' type='number' min='1' max='10000' value='{int(ga)}' required></label>"
                    f"<small>This controls search repetitions inside each compound docking job. It is independent of CPU count.</small>"
                    f"<label>Parallel compound jobs <input name='parallel_jobs' value='{h.escape(str(parallel))}'></label>"
                    f"<small>Use <b>auto</b> to use up to the {int(detected)} CPUs available to SIDock, one AutoDock compound job per worker.</small><br>"
                    f"<button type='submit'>RUN</button></form>")
        elif dash['primary_action']=='RESUME':
            action=f"<form method='post' action='{h.escape(dash['primary_action_url'])}'><button type='submit'>RESUME</button></form>"
        else:
            action=f"<a href='{h.escape(dash['primary_action_url'] or '#')}'>{h.escape(dash['primary_action'])}</a>"
        if (get_project(project_id)/'manual_docking/workflow.json').exists():
            action=(f"<p><a href='/projects/{h.escape(project_id)}/guided/1'>"
                    "<b>Continue manually reviewed five-step docking workflow</b></a></p>")
        else:
            action=(f"<p><a href='/projects/{h.escape(project_id)}/guided/1'>"
                    "Switch to reviewed five-step docking workflow</a></p>" + action)
        failed_stages = [
            (stage, row) for stage, row in state.get('stages', {}).items()
            if isinstance(row, dict) and row.get('status') in {'FAIL', 'BLOCKED'}
        ]
        errors_html = ''.join(
            f"<tr><td>{h.escape(str(stage))}</td><td>{h.escape(str(row.get('status')))}</td>"
            f"<td>{h.escape('; '.join(map(str, row.get('errors') or [])) or 'No error detail recorded')}</td></tr>"
            for stage, row in failed_stages
        )
        diagnostics = (
            "<h2>Why did the execution stop?</h2><table><tr><th>Stage</th><th>Status</th><th>Recorded error</th></tr>"
            + errors_html + "</table>"
        ) if failed_stages else ''
        reference_resolved = bool(counts.get('selection', {}).get('selected_instance')) if isinstance(counts.get('selection'), dict) else False
        reference_artifacts_status = 'reference ligand selected' if reference_resolved else 'no valid reference selection recorded'
        if not reference_resolved:
            diagnostics += (
                f"<p><b>Reference unresolved.</b> Check the recorded reference-stage error above and "
                f"<a href='/projects/{h.escape(project_id)}/reference'>reference candidates</a> "
                "before retrying. Intake metadata alone is not a selected ligand.</p>"
            )
        body=f"<h1>{h.escape(Path(data['project']).name)}</h1><p>Scientific status: <span class='badge'>{h.escape(final)}</span></p><p>Operational job: <span class='badge'>{h.escape(str(dash.get('operational_job_status') or 'No GUI job recorded'))}</span></p><p>Reference: {h.escape(str(dash['reference_status']))}; validation stage: {h.escape(str(dash['validation_status']))}; validation policy: <b>{h.escape(str(dash.get('validation_policy','required')))}</b></p><p>Docking setup: <b>{h.escape(str(dash.get('ga_runs',10)))}</b> GA runs per compound; parallel jobs: <b>{h.escape(str(dash.get('parallel_jobs','auto')))}</b>; detected CPUs: <b>{h.escape(str(dash.get('detected_cpus',1)))}</b>.</p><p>Input compounds: {dash['input_compound_count']}; execution records: {dash['compound_total']} ({dash['compound_success']} successful, {dash['compound_failed']} failed, {dash['compound_blocked']} blocked).</p><p>{h.escape(dash['explanatory_message'])}</p>{action}<p>Reference: {reference_artifacts_status}</p>{diagnostics}<nav><a href='/projects/{project_id}/manual'>Manual Docking</a><a href='/projects/{project_id}/reference'>Reference</a><a href='/projects/{project_id}/compounds'>Compounds</a><a href='/projects/{project_id}/results'>Results</a><a href='/projects/{project_id}/jobs'>Jobs</a><a href='/projects/{project_id}/config'>Config</a></nav>"
        if final=='BLOCKED': body += '<p class="error">Strict native redocking validation blocked compound docking.</p>'
        if data.get('active_job'): body += '<p>Project currently executing.</p>'
        return page(Path(data['project']).name,body)
    @app.get('/projects/{project_id}/manual', response_class=HTMLResponse)
    def manual_docking_page(project_id):
        from .manual_docking import (
            settings as manual_settings, root as manual_root, job_status,
            active_run, prepared_ligands, _input_metals, _pdbqt_metals,
        )
        import html as h
        p = get_project(project_id)
        t = manual_settings(p)
        m = manual_root(p)
        receptor = m / 'receptor.pdbqt'
        original_metals = _input_metals(p)
        receptor_metals = _pdbqt_metals(p)
        job = job_status(p)
        run = active_run(p)
        ligands = prepared_ligands(p)
        exp = lambda value: h.escape(str(value), quote=True)
        fields = ''.join(
            f"<label>{label} <input type='number' step='{step}' name='{name}' "
            f"value='{exp(t[name])}' required></label>"
            for name, label, step in (
                ('center_x', 'Center X (Å)', '0.001'),
                ('center_y', 'Center Y (Å)', '0.001'),
                ('center_z', 'Center Z (Å)', '0.001'),
                ('npts_x', 'Grid points X (20–126, even)', '2'),
                ('npts_y', 'Grid points Y (20–126, even)', '2'),
                ('npts_z', 'Grid points Z (20–126, even)', '2'),
                ('spacing', 'Grid spacing (Å)', '0.025'),
                ('ga_runs', 'GA runs per compound', '1'),
                ('parallel_jobs', 'Simultaneous AutoDock jobs', '1'),
                ('ga_num_evals', 'Maximum evaluations', '1'),
                ('ga_pop_size', 'GA population size', '1'),
                ('ga_num_generations', 'Maximum generations', '1'),
            )
        )
        metal_warning = (
            f"<p class='error'><b>Metal audit:</b> the original structure contains {original_metals} Ni ions; "
            f"the uploaded receptor contains {receptor_metals}. "
            "A metal-free model is not equivalent to the crystallographic urease active site. "
            "Any such docking must be explicitly labeled exploratory and cannot validate nickel coordination.</p>"
            if original_metals and receptor_metals != original_metals else ''
        )
        user_warning = (
            "<p><b>Important:</b> A valid receptor PDBQT must already have chemically appropriate "
            "charges and atom types. SIDock does not invent nickel parameters. "
            "This manual workflow does not perform native redocking or certify the receptor.</p>"
        )
        action_status = (
            f"<p><b>Last action:</b> {exp(job.get('action'))}: "
            f"{exp(job.get('status'))}. {exp(job.get('error',''))}</p>"
            if job else '<p>No manual action has been launched.</p>'
        )
        results = ''
        if run and run.is_dir():
            import json as _json
            grid_file = run / 'grid/production/grid_provenance.json'
            dlg_files = sorted((run/'docking/dlg').glob('*.dlg'))
            root_rel = run.relative_to(p).as_posix()
            results = (
                f"<p><b>Run ID:</b> {exp(run.name)} "
                f"| Grid: {'CREATED' if grid_file.exists() else 'NOT CREATED'} "
                f"| DLG files: {len(dlg_files)}/{len(ligands)}</p>"
            )
            if grid_file.exists():
                results += f"<p><a href='/projects/{project_id}/artifacts/{root_rel}/grid/production/receptor.gpf'>View grid parameters (.gpf)</a></p>"
            if (run/'results.csv').is_file():
                results += f"<p><a href='/projects/{project_id}/artifacts/{root_rel}/results.csv'>Download manual results (.csv)</a></p>"
            if dlg_files:
                results += '<details><summary>Docking logs</summary><ul>'
                for f in dlg_files:
                    results += (f"<li><a href='/projects/{project_id}/artifacts/"
                                f"{root_rel}/docking/dlg/{exp(f.name)}'>{exp(f.name)}</a></li>")
                results += '</ul></details>'
        sizes = [t[k]*t['spacing'] for k in ('npts_x','npts_y','npts_z')]
        body = (
            f"<h1>Manual AutoDock4: {exp(project_id)}</h1>"
            "<p>This is a separate manual workflow. It does not restart the automatic pipeline "
            "or discard any existing project docking outputs.</p>"
            "<h2>1. Prepared receptor PDBQT</h2>"
            f"<p>Current receptor: {exp(receptor) if receptor.is_file() else 'Not imported'}</p>"
            f"<form method='post' action='/projects/{project_id}/manual/receptor' enctype='multipart/form-data'>"
            "<label>Upload prepared receptor .pdbqt <input type='file' name='upload' accept='.pdbqt'></label>"
            "<label>Or import existing project-relative .pdbqt path "
            "<input name='existing_path' placeholder='provenance/4ubp_protein_diagnostic.pdbqt'></label>"
            "<button>Use this receptor</button></form>"
            + metal_warning + user_warning
            + "<h2>2. Grid box and docking parameters</h2>"
            + f"<p>Approximate box dimensions: {sizes[0]:.2f} × {sizes[1]:.2f} × {sizes[2]:.2f} Å</p>"
            + "<p>Change the center and grid-point counts yourself. The initial center is estimated "
              "from the crystal HAE coordinates when available.</p>"
            + f"<form method='post' action='/projects/{project_id}/manual/settings'>"
            + fields
            + "<label><input type='checkbox' name='allow_missing_metal' value='on' "
              + ("checked" if t['allow_missing_metal'] else "")
              + "> I acknowledge this is exploratory metal-free docking if the catalytic Ni ions are absent</label>"
            + "<button>Save manual parameters</button></form>"
            + "<h2>3. Execute each step separately</h2>"
            + f"<p>Prepared ligand files: {len(ligands)}</p>"
        )
        for action, label in (
            ('prepare','Prepare ligands from existing project SDF (skip completed preparations)'),
            ('grid','Generate AutoGrid4 maps for the saved box'),
            ('dock','Run AutoDock4 (skip completed compounds)'),
        ):
            body += (f"<form method='post' action='/projects/{project_id}/manual/action'>"
                     f"<input type='hidden' name='action' value='{action}'>"
                     f"<button type='submit'>{label}</button></form>")
        body += action_status + results
        body += (f"<p><a href='/projects/{project_id}/manual'>Refresh manual status</a> "
                 f"| <a href='/projects/{project_id}/manual/log'>View manual worker log</a> "
                 f"| <a href='/projects/{project_id}'>Back to project</a></p>")
        return page('Manual Docking - '+project_id, body)

    @app.post('/projects/{project_id}/manual/settings')
    def manual_settings_save(project_id, center_x:float=Form(...), center_y:float=Form(...),
                             center_z:float=Form(...), npts_x:int=Form(...), npts_y:int=Form(...),
                             npts_z:int=Form(...), spacing:float=Form(...), ga_runs:int=Form(...),
                             parallel_jobs:int=Form(...), ga_num_evals:int=Form(...),
                             ga_pop_size:int=Form(...), ga_num_generations:int=Form(...),
                             allow_missing_metal:str|None=Form(None)):
        from .manual_docking import save_settings
        try:
            save_settings(get_project(project_id), locals())
        except (ValueError, TypeError) as exc:
            raise HTTPException(status_code=400, detail=str(exc))
        from fastapi.responses import RedirectResponse
        return RedirectResponse(f'/projects/{project_id}/manual', status_code=303)

    @app.post('/projects/{project_id}/manual/receptor')
    async def manual_receptor_select(project_id, upload: UploadFile|None=File(None),
                                     existing_path:str|None=Form(None)):
        from .manual_docking import install_receptor, root as manual_root, job_status
        p = get_project(project_id)
        if job_status(p).get('status') == 'RUNNING':
            raise HTTPException(status_code=409, detail='Manual action is running')
        try:
            if upload is not None and upload.filename:
                if not upload.filename.lower().endswith('.pdbqt'):
                    raise ValueError('Only .pdbqt files are supported')
                dest = manual_root(p) / 'uploaded_receptor.pdbqt'
                dest.parent.mkdir(parents=True, exist_ok=True)
                total = 0
                with dest.open('wb') as f:
                    while True:
                        chunk = await upload.read(1024 * 1024)
                        if not chunk: break
                        total += len(chunk)
                        if total > MAX_UPLOAD_BYTES: raise ValueError('PDBQT exceeds allowed upload size')
                        f.write(chunk)
                install_receptor(p, dest)
            elif existing_path and existing_path.strip():
                target = (p / existing_path.strip()).resolve()
                if not target.is_relative_to(p.resolve()):
                    raise ValueError('Only project-local receptor files are allowed')
                install_receptor(p, target)
            else:
                raise ValueError('Upload a prepared PDBQT or enter an existing path')
        except (OSError, ValueError) as exc:
            raise HTTPException(status_code=400, detail=str(exc))
        from fastapi.responses import RedirectResponse
        return RedirectResponse(f'/projects/{project_id}/manual', status_code=303)

    @app.post('/projects/{project_id}/manual/action')
    def manual_action(project_id, action:str=Form(...)):
        from .manual_docking import launch
        try:
            launch(get_project(project_id), action)
        except (ValueError, RuntimeError) as exc:
            raise HTTPException(status_code=409, detail=str(exc))
        from fastapi.responses import RedirectResponse
        return RedirectResponse(f'/projects/{project_id}/manual', status_code=303)

    @app.get('/projects/{project_id}/manual/log', response_class=HTMLResponse)
    def manual_log(project_id):
        from .manual_docking import root as manual_root
        import html as h
        p=get_project(project_id)
        path=manual_root(p)/'manual_worker.log'
        data=path.read_text(encoding='utf-8', errors='replace')[-16000:] if path.exists() else 'No manual worker log yet.'
        return page('Manual worker log', f"<h1>Manual Docking log</h1><pre>{h.escape(data)}</pre>"
                    f"<p><a href='/projects/{project_id}/manual'>Back to manual docking</a></p>")

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
    @app.get('/projects/{project_id}/reference-inputs', response_class=HTMLResponse)
    def reference_inputs(project_id):
        return page('Reference inputs', f"<h1>Reference inputs</h1><p>Attach or replace validation structure, mmCIF metadata, or CCD chemistry.</p><form method='post' enctype='multipart/form-data'><label>Validation structure <input type='file' name='validation_structure' accept='.pdb,.ent,.cif,.mmcif'></label><label>mmCIF <input type='file' name='mmcif' accept='.cif,.mmcif'></label><label>CCD JSON <input type='file' name='ccd' accept='.json,.cif'></label><button>Update reference inputs</button></form>")
    @app.post('/projects/{project_id}/reference-inputs')
    async def reference_inputs_update(project_id, validation_structure: UploadFile|None=File(None), mmcif: UploadFile|None=File(None), ccd: UploadFile|None=File(None)):
        p=get_project(project_id); upload=root/'.uploads'; upload.mkdir(exist_ok=True); files=[]
        try:
            for item,suffix in ((validation_structure,'.pdb'),(mmcif,'.cif'),(ccd,'.json')):
                if item is None: files.append(None); continue
                path=upload/f'{p.name}-update-{len(files)}{suffix}'; total=0
                with path.open('wb') as fh:
                    while True:
                        chunk=await item.read(1024*1024)
                        if not chunk: break
                        total+=len(chunk)
                        if total>MAX_UPLOAD_BYTES: raise ValueError('E001 UPLOAD_TOO_LARGE')
                        fh.write(chunk)
                files.append(path)
            attach_project_inputs(p,validation_structure=files[0],mmcif=files[1],ccd=files[2])
            from fastapi.responses import RedirectResponse
            return RedirectResponse(f'/projects/{project_id}',status_code=303)
        except (ValueError,RuntimeError) as exc:
            raise HTTPException(status_code=400,detail=str(exc))
    @app.post('/projects/{project_id}/reference')
    def choose_reference(project_id, instance: str):
        try: return set_reference_override(get_project(project_id), instance)
        except ValueError as exc: raise HTTPException(status_code=400, detail=str(exc))
    @app.post('/projects/{project_id}/run')
    def run(project_id, request: Request, validation_mode: str=Form('required'), ga_runs: int=Form(10), parallel_jobs: str=Form('auto')):
        try:
            result=start_run(get_project(project_id), validation_mode=validation_mode, ga_runs=ga_runs, parallel_jobs=parallel_jobs)
            if 'text/html' in request.headers.get('accept','') or request.headers.get('content-type','').startswith('application/x-www-form-urlencoded'):
                from fastapi.responses import RedirectResponse
                return RedirectResponse(f'/projects/{project_id}',status_code=303)
            return result
        except RuntimeError as exc: raise HTTPException(status_code=409,detail=str(exc))
    @app.post('/projects/{project_id}/resume')
    def resume(project_id, request: Request):
        try:
            result=resume_project(get_project(project_id))
            if 'text/html' in request.headers.get('accept','') or request.headers.get('content-type','').startswith('application/x-www-form-urlencoded'):
                from fastapi.responses import RedirectResponse
                return RedirectResponse(f'/projects/{project_id}',status_code=303)
            return result
        except RuntimeError as exc: raise HTTPException(status_code=409,detail=str(exc))
    @app.get('/projects/{project_id}/compounds',response_class=HTMLResponse)
    def compounds(project_id, page_num:int=Query(1, alias='page', ge=1), page_size:int=25, status:str='ALL', search:str='', sort:str='compound_id'):
        data=get_compound_page(get_project(project_id),page=page_num,page_size=page_size,status=status,search=search,sort=sort); body=[f"<h1>Compounds</h1><form><label>Search <input name='search' value='{__import__('html').escape(search)}'></label><label>Status <select name='status'><option>ALL</option><option>PASS</option><option>FAIL</option><option>BLOCKED</option><option>STALE</option></select></label><input type='hidden' name='page_size' value='{data['page_size']}'><button>Filter</button></form><p>Showing {min((data['page']-1)*data['page_size']+1,data['total']) if data['total'] else 0}–{min(data['page']*data['page_size'],data['total'])} of {data['total']}</p><table><tr><th>ID</th><th>Status</th><th>Docking</th><th>Details</th></tr>"]
        for row in data['items']:
            cid=__import__('html').escape(str(row['compound_id'])); body.append(f"<tr><td>{cid}</td><td>{__import__('html').escape(row['status'])}</td><td>{__import__('html').escape(str(row.get('docking_energy') if row.get('docking_energy') is not None else 'Unavailable'))}</td><td><a href='/projects/{project_id}/compounds/{cid}'>Details</a></td></tr>")
        if not data['total']:
            count,_=input_compound_count(get_project(project_id))
            return page('Compounds - '+project_id,f'<h1>Compounds</h1><p>Compound execution has not started yet. {count} compounds were found in the input SDF.</p>')
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
        if data['project_status']=='BLOCKED': blocked='<p><b>Strict native redocking validation failed; compound docking was not performed.</b></p>'
        elif validation.get('docking_validation_status')=='UNVALIDATED': blocked='<p><b>UNVALIDATED DOCKING:</b> native redocking was skipped or did not satisfy the validation criterion. Interpret docking poses and energies as exploratory predictions.</p>'
        links=[]
        for rel,label in (('results.csv','results.csv'),('analysis/evidence.json','evidence.json'),('report.html','report.html'),('analysis/scientific_report.html','Scientific analysis and write-up'),('analysis/scientific_report.md','Editable write-up'),('analysis/scientific_report.json','Analysis JSON')):
            if (get_project(project_id)/rel).is_file(): links.append(f"<a href='/projects/{project_id}/artifacts/{rel}'>{label}</a>")
        empty='' if links else '<p>Final results are not available yet.</p>'
        return page('Results - '+project_id,f"<h1>Results</h1><p>Project status: <b>{_html.escape(str(data['project_status']))}</b></p>{blocked}<p>Validation policy: {_html.escape(str(validation.get('validation_policy','Unavailable')))}; validation status: {_html.escape(str(validation.get('docking_validation_status','Not run yet')))}; native result: {_html.escape(str(validation.get('sampling_result','Not run yet')))}; RMSD {_html.escape(str(validation.get('native_redocking_RMSD_A','Unavailable')))}; threshold {_html.escape(str(validation.get('rmsd_threshold','Unavailable')))}</p><p>Counts: {_html.escape(str(counts))}</p>{empty}{' | '.join(links)}")

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
