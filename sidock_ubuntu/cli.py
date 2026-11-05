from __future__ import annotations
import argparse
import json, shutil, sys, importlib.util
from pathlib import Path
from .pipeline import run_stage_runner
from .config import load_config, summary, validate_config
from .dashboard import write_dashboard, write_results
from .state import project_status, STAGES, invalidate_for_config_change, _atomic_write, derived_stage_statuses, derived_project_status
from .compliance import evaluate
from .orchestrator import StageRunner, config_change_audit
from .executors import production_context_from_project, build_production_executors
from .project_lock import ProjectExecutionLocked
from . import __version__
from .dependencies import resolve_plip_executable


def parser() -> argparse.ArgumentParser:
    p=argparse.ArgumentParser(prog='sidock',description='Automated AutoDock4 molecular docking workflow')
    sub=p.add_subparsers(dest='command',required=True)
    r=sub.add_parser('run',help='Run receptor + SDF through the complete docking workflow')
    r.add_argument('receptor',type=Path)
    r.add_argument('ligands',type=Path)
    r.add_argument('--project',type=Path,default=Path('sidock_project'))
    r.add_argument('--padding',type=float,default=8.0)
    r.add_argument('--ga-runs',type=int,default=None)
    r.add_argument('--jobs',type=int,default=1)
    r.add_argument('--top-n',type=int,default=20)
    r.add_argument('--seed',type=int,default=20250301)
    r.add_argument('--allow-failed-redocking-validation', action='store_true', help='Developer/testing override; continue after failed native redocking validation')
    r.add_argument('--config', type=Path, default=None)
    r.add_argument('--reference-ligand', default=None, help='Explicit ligand instance, e.g. P43:D:1263')
    r.add_argument('--allow-incomplete-peripheral', action='store_true', help='Allow audited incomplete residues 5–8 Å from reference ligand')
    sub.add_parser('doctor', help='Check the SIDock runtime and docking executables')
    sr=sub.add_parser('report', help='Generate evidence-grounded analysis and write-up without redocking')
    sr.add_argument('project', type=Path)
    for name in ('status', 'resume'):
        q = sub.add_parser(name, help=f'{name.title()} a project')
        q.add_argument('project', type=Path)
        if name == 'resume':
            q.add_argument('--force', action='store_true')
            q.add_argument('--failed-only', action='store_true', help='Retry only recorded failed compound stages')
            q.add_argument('--compound', default=None, help='Retry one compound stage')
            q.add_argument('--dry-run', action='store_true', help='Show invalidation/resume plan without executing')
    c = sub.add_parser('config', help='Show or validate configuration')
    c.add_argument('action', choices=('show','validate'))
    c.add_argument('--file', type=Path, default=None)
    v = sub.add_parser('validate', help='Validate input availability and configuration')
    v.add_argument('receptor', type=Path); v.add_argument('--config', type=Path, default=None)
    i = sub.add_parser('inspect', help='Inspect a project or compound')
    i.add_argument('project', type=Path); i.add_argument('--compound', default=None)
    c = sub.add_parser('compliance', help='Evaluate machine-checkable production rules')
    c.add_argument('project', type=Path)
    return p

def doctor_status() -> dict:
    checks = [("Python", sys.executable), ("RDKit", importlib.util.find_spec("rdkit")),
              ("AutoGrid4", shutil.which("autogrid4")), ("AutoDock4", shutil.which("autodock4")),
              ("Open Babel", shutil.which("obabel")), ("PyMOL", shutil.which("pymol"))]
    items=[]
    for label,value in checks:
        available=bool(value); detail=value if isinstance(value,str) else ('available' if value else 'not detected')
        items.append({'name':label,'status':'AVAILABLE' if available else 'MISSING','detail':detail,'required':label not in {'PLIP','PyMOL'}})
    module=importlib.util.find_spec('plip'); items.append({'name':'PLIP Python module','status':'AVAILABLE' if module else 'MISSING','detail':str(module) if module else 'not detected','required':True})
    plip=resolve_plip_executable(); items.append({'name':'PLIP executable','status':'AVAILABLE' if plip['available'] else 'MISSING','detail':plip.get('path') or plip['message'],'version':plip.get('version'),'required':True})
    return {'status':'READY' if all(x['status']=='AVAILABLE' for x in items if x['required']) else 'NOT_READY','dependencies':items}

def _doctor() -> int:
    report = doctor_status()
    for item in report['dependencies']:
        detail = item.get('detail') or ''
        version = f" ({item['version']})" if item.get('version') else ''
        print(f"[{'PASS' if item['status']=='AVAILABLE' else 'FAIL'}] {item['name']}{version}{' ' + str(detail) if detail else ''}")
    return 0 if report['status'] == 'READY' else 1

def _print_status(project: Path) -> None:
    data = project_status(project)
    stages = derived_stage_statuses(data)
    print(f"Project status: {derived_project_status(data)}")
    for stage in STAGES: print(f"[{stages.get(stage,{}).get('status','PENDING')}] {stage}")


def main() -> None:
    args=parser().parse_args()
    if args.command=='run':
        if args.jobs > 1:
            print('StageRunner parallel execution is not enabled in v0.4; continuing serially with --jobs 1.', file=sys.stderr)
            args.jobs = 1
        config = load_config(args.config)
        if args.ga_runs is not None: config['validation']['ga_runs'] = args.ga_runs
        errors = validate_config(config)
        if errors: raise SystemExit('Invalid configuration: ' + '; '.join(errors))
        print(summary(config))
        try:
            result = run_stage_runner(args.project,args.receptor,args.ligands,config=config,padding=args.padding,parallel_jobs=args.jobs,top_n=args.top_n,seed=args.seed,reference_ligand=args.reference_ligand,allow_incomplete_peripheral=args.allow_incomplete_peripheral)
        except ProjectExecutionLocked as exc:
            raise SystemExit(str(exc))
        try:
            import yaml
            p = args.project / 'provenance'; p.mkdir(parents=True, exist_ok=True)
            (p / 'resolved_config.yaml').write_text(yaml.safe_dump(config, sort_keys=False), encoding='utf-8')
        except ImportError:
            pass
        write_results(args.project); write_dashboard(args.project)
    elif args.command == 'doctor':
        raise SystemExit(_doctor())
    elif args.command == 'report':
        from .scientific_report import write_scientific_report
        print(json.dumps(write_scientific_report(args.project),indent=2))
    elif args.command == 'status':
        _print_status(args.project)
    elif args.command == 'resume':
        state=project_status(args.project)
        current_config=load_config(); previous_config=current_config
        resolved=args.project/'provenance/resolved_config.yaml'
        if resolved.exists():
            try:
                import yaml
                previous_config=yaml.safe_load(resolved.read_text(encoding='utf-8')) or current_config
            except Exception:
                previous_config=current_config
        audit=config_change_audit(args.project,previous_config,current_config)
        if audit.get('changed'):
            state=invalidate_for_config_change(state,previous_config,current_config); _atomic_write(args.project/'project_state.json',state)
        compounds=list(state.get('compounds',{})); ctx=production_context_from_project(args.project, current_config); runner=StageRunner(args.project, current_config, build_production_executors(ctx))
        if compounds:
            selected=[args.compound] if args.compound else compounds
            plan=runner.plan(selected,failed_only=args.failed_only,force=args.force)
            print(json.dumps(plan,indent=2))
            if args.dry_run: return
            for cid in selected:
                if all(v.get('status')=='PASS' for v in state.get('compounds',{}).get(cid,{}).values() if isinstance(v,dict) and 'status' in v) and not args.force:
                    print(f'Compound {cid} is already complete.'); continue
                try:
                    runner.run_compound(cid,failed_only=args.failed_only,force=args.force)
                except ProjectExecutionLocked as exc:
                    raise SystemExit(str(exc))
            evaluate(args.project); return
        if (args.project / 'project_state.json').exists() and all(v.get('status') == 'PASS' for v in state.get('stages',{}).values()):
            _print_status(args.project); print('Nothing to resume; completed outputs were not rerun.')
        else:
            project = json.loads((args.project/'project.json').read_text())
            if not args.force: raise SystemExit('Resume requires stage-aware checkpoints; use --force to rerun the legacy pipeline explicitly.')
            run_stage_runner(args.project, Path(project['receptor_source']), Path(project['ligand_source']), config=load_config(), parallel_jobs=project.get('docking',{}).get('parallel_jobs',1))
    elif args.command == 'config':
        config = load_config(args.file)
        if args.action == 'show': print(json.dumps(config, indent=2))
        else:
            errors = validate_config(config); print('VALID' if not errors else '\n'.join(errors)); raise SystemExit(0 if not errors else 1)
    elif args.command == 'validate':
        config = load_config(args.config); errors = validate_config(config)
        if not args.receptor.exists(): errors.append('E001 INPUT_FILE_MISSING: receptor does not exist')
        print('VALID' if not errors else '\n'.join(errors)); raise SystemExit(0 if not errors else 1)
    elif args.command == 'inspect':
        print(json.dumps(project_status(args.project), indent=2))
    elif args.command == 'compliance':
        result=evaluate(args.project); print(json.dumps(result, indent=2)); raise SystemExit(0 if result['status']=='PASS' else 1)


if __name__=='__main__':
    main()
