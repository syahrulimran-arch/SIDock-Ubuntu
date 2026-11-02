"""Adapters between the stage orchestrator and the validated scientific backend.

These adapters intentionally contain no chemistry or docking algorithms.  They
provide a stable seam for the existing production functions and are exercised
with injected callbacks in orchestration tests.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Any
from .orchestrator import StageResult
from .stages import (run_reference_stage, run_receptor_stage, run_grid_stage,
                     run_validation_stage, run_compound_preparation,
                     run_compound_docking, run_pose_extraction_stage,
                     run_cluster_stage, run_plip_stage, run_figures_stage,
                     run_finalization_stage)
from .core import box_from_reference

@dataclass
class ProductionContext:
    project: Path
    config: dict
    receptor: Path | None = None
    ligand_sdf: Path | None = None
    callbacks: dict[str, Callable[..., Any]] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def rehydrate_global_outputs(self):
        from .stage_outputs import (load_reference_output, load_receptor_output,
                                     load_grid_output, load_validation_output)
        loaders=(('reference',load_reference_output),('receptor',load_receptor_output),
                 ('grid',load_grid_output),('validation',load_validation_output))
        for name, loader in loaders:
            try: self.metadata[name+'_output']=loader(self.project)
            except FileNotFoundError: pass
        return self.metadata

    def callback(self, name: str, *args, **kwargs):
        """Call an existing validated operation supplied by the production adapter."""
        fn=self.callbacks.get(name)
        if fn is None:
            return None
        return fn(*args, **kwargs)

def execute_compound_preparation(ctx: ProductionContext, compound: str) -> StageResult:
    result=ctx.callback('prepare_ligand', compound)
    if result is None: return StageResult('BLOCKED', error_code='E401', error_message='No production ligand-preparation binding supplied')
    return result if isinstance(result,StageResult) else StageResult('PASS', outputs=result)

def execute_compound_docking(ctx: ProductionContext, compound: str) -> StageResult:
    result=ctx.callback('dock_compound', compound)
    if result is None: return StageResult('BLOCKED', error_code='E502', error_message='No production AutoDock binding supplied')
    return result if isinstance(result,StageResult) else StageResult('PASS', outputs=result)

def execute_pose_extraction(ctx: ProductionContext, compound: str) -> StageResult:
    result=ctx.callback('extract_poses', compound)
    if result is None: return StageResult('BLOCKED', error_code='E602', error_message='No corrected DLG traversal extractor binding supplied')
    return result if isinstance(result,StageResult) else StageResult('PASS', outputs=result)

def execute_clustering(ctx: ProductionContext, compound: str) -> StageResult:
    result=ctx.callback('cluster_dlg', compound)
    if result is None: return StageResult('BLOCKED', error_code='E601', error_message='No raw-DLG cluster parser binding supplied')
    return result if isinstance(result,StageResult) else StageResult('PASS', outputs=result)

def execute_plip(ctx: ProductionContext, compound: str) -> StageResult:
    result=ctx.callback('plip', compound)
    if result is None: return StageResult('BLOCKED', error_code='E701', error_message='No PLIP production binding supplied')
    return result if isinstance(result,StageResult) else StageResult('PASS', outputs=result)

def execute_figures(ctx: ProductionContext, compound: str) -> StageResult:
    result=ctx.callback('figures', compound)
    if result is None: return StageResult('BLOCKED', error_code='E999', error_message='No figure production binding supplied')
    return result if isinstance(result,StageResult) else StageResult('PASS', outputs=result)

def bind_production_executors(ctx: ProductionContext) -> dict[str, Callable[[str], StageResult]]:
    """Bind StageRunner to the existing validated stage seams.

    This adapter contains path/context translation only. Scientific work remains
    implemented by ``stages.py`` and its existing backend functions.
    """
    def call(name, *args, **kwargs):
        fn = ctx.callbacks.get(name)
        if fn is None:
            return StageResult('BLOCKED', error_code='E999', error_message=f'Missing production binding: {name}')
        try:
            value = fn(*args, **kwargs)
            return value if isinstance(value, StageResult) else StageResult('PASS', outputs=value if isinstance(value, dict) else {'result': value})
        except Exception as exc:
            code = {'prepare': 'E401', 'dock': 'E502', 'pose': 'E602', 'cluster': 'E601', 'plip': 'E701'}.get(name, 'E999')
            return StageResult('FAIL', error_code=code, error_message=str(exc))

    return {
        'preparation': lambda cid: call('prepare', cid),
        'docking': lambda cid: call('dock', cid),
        'pose_extraction': lambda cid: call('pose', cid),
        'clustering': lambda cid: call('cluster', cid),
        'plip': lambda cid: call('plip', cid),
        'figures': lambda cid: call('figures', cid),
    }

def build_production_executors(ctx: ProductionContext) -> dict[str, Callable[..., StageResult]]:
    """Return the single registry shared by fresh runs and resume."""
    registry = bind_production_executors(ctx)
    for name in ('reference','receptor','grid','validation','finalization'):
        callback_name = name
        registry[name] = lambda *args, _name=callback_name, **kwargs: _global_result(ctx, _name, *args, **kwargs)
    return registry

def _global_result(ctx: ProductionContext, name: str, *args, **kwargs) -> StageResult:
    fn=ctx.callbacks.get(name) or _default_global_callback(ctx,name)
    if fn is None:
        return StageResult('BLOCKED', error_code='E999', error_message=f'Missing production binding: {name}')
    try:
        value=fn(*args, **kwargs)
        result=value if isinstance(value,StageResult) else StageResult('PASS', outputs=value if isinstance(value,dict) else {'result':value})
        if result.status in {'PASS','FAIL'} and isinstance(result.outputs,dict):
            _persist_global_result(ctx,name,result)
        return result
    except Exception as exc:
        return StageResult('FAIL', error_code='E999', error_message=str(exc))

def _default_global_callback(ctx: ProductionContext, name: str):
    """Construct the real global seam adapter from persisted/bootstrap context."""
    from .intake import hetero_groups, choose_reference_ligand
    from .core import pdbqt_atom_types, resolve_executable
    from .validation import prepare_native_ligand
    if name=='reference':
        def reference():
            groups=hetero_groups(ctx.receptor)
            # Persist the chemistry-aware candidate audit when the local PDB
            # contains discoverable candidates.  Minimal seam fixtures may not
            # have HETATM records, so legacy backend tests remain compatible.
            from .intake import reference_candidates, receptor_completeness
            from .reference import write_reference_selection, sha256_file, load_local_ccd, load_local_mmcif, reference_config_fingerprint, map_instance_to_ccd
            try:
                candidates=reference_candidates(ctx.receptor) if ctx.receptor and ctx.receptor.exists() else []
            except (AttributeError, ValueError):
                # Test/backend adapters may expose only the minimal group
                # identity needed by the existing reference seam.
                candidates=[]
            selection_artifacts={}; ref_cfg=ctx.config.get('reference',{})
            selected=next((g for g in groups if g.key==ref_cfg.get('ligand')),None) or choose_reference_ligand(ctx.receptor)
            ccd_source=load_local_ccd(ref_cfg.get('ccd_path'), selected.resname if selected else None) if ref_cfg.get('ccd_path') else {'available':False,'parse_status':'MISSING','warnings':['CCD_SOURCE_NOT_CONFIGURED']}
            mmcif_source=load_local_mmcif(ref_cfg.get('mmcif_path')) if ref_cfg.get('mmcif_path') else {'available':False,'parse_status':'MISSING','warnings':['MMCIF_SOURCE_NOT_CONFIGURED']}
            if candidates:
                for row in candidates:
                    if ccd_source.get('available'):
                        row.update(map_instance_to_ccd(row,ccd_source.get('chemistry')))
                    row['ccd_source_hash']=ccd_source.get('source_hash')
                    mm_key=':'.join(str(row.get(k,'')) for k in ('ligand_code','chain','residue_number','insertion_code'))
                    mm=mmcif_source.get('instances',{}).get(mm_key)
                    if mm:
                        row.update({k:mm.get(k) for k in ('occupancy_summary','b_factor_summary','rscc','rsr','model','altloc')})
                    group=next((g for g in groups if getattr(g,'key',None)==row.get('key')),None)
                    if group is not None:
                        local=receptor_completeness(ctx.receptor,group)
                        binding=[x for x in local if x.get('within_5A')]
                        peripheral=[x for x in local if x.get('within_8A') and not x.get('within_5A')]
                        contact=[x for x in local if x.get('within_5A') or x.get('within_8A')]
                        bad=lambda xs: sum(1 for x in xs if str(x.get('classification','')).startswith('INCOMPLETE'))
                        row.update({'binding_site_complete':bad(binding)==0,'binding_site_incomplete_residue_count':bad(binding),'binding_site_incomplete_atom_count':sum(len(str(x.get('missing_heavy_atoms','')).split(';')) for x in binding if x.get('missing_heavy_atoms')),'peripheral_incomplete_residue_count':bad(peripheral),'peripheral_incomplete_atom_count':sum(len(str(x.get('missing_heavy_atoms','')).split(';')) for x in peripheral if x.get('missing_heavy_atoms')),'contact_environment_complete':bad(contact)==0,'contact_environment_complete_fraction':(len(contact)-bad(contact))/len(contact) if contact else 1.0})
                selection_artifacts=write_reference_selection(ctx.project,candidates,override=ref_cfg.get('ligand'),input_hashes={'receptor':sha256_file(ctx.receptor),'reference_config':reference_config_fingerprint(ctx.config),'ccd':ccd_source.get('source_hash'),'mmcif':mmcif_source.get('source_hash')},provenance={'source':'local_pdb','ccd':ccd_source,'mmcif':mmcif_source})
            if selected is None: raise RuntimeError('E101 REFERENCE_LIGAND_UNRESOLVED')
            out=run_reference_stage(ctx.receptor,ctx.project,selected)
            out.setdefault('selected_ligand',selected.key); out.setdefault('resname',selected.resname); out.setdefault('chain',selected.chain); out.setdefault('resseq',selected.resseq); out.setdefault('reference_pdb',out.get('reference_pdb'))
            out.update({'candidate_csv':selection_artifacts.get('candidate_csv'),'selection_json':selection_artifacts.get('selection_json'),'source_metadata':{'ccd':ccd_source,'mmcif':mmcif_source},'input_hashes':{str(ctx.receptor):sha256_file(ctx.receptor),'reference_config':reference_config_fingerprint(ctx.config),**({str(ref_cfg['ccd_path']):ccd_source.get('source_hash')} if ccd_source.get('source_hash') else {}),**({str(ref_cfg['mmcif_path']):mmcif_source.get('source_hash')} if mmcif_source.get('source_hash') else {})}})
            return out
        return reference
    if name=='receptor':
        def receptor():
            ctx.rehydrate_global_outputs(); ref=ctx.metadata.get('reference_output')
            if not ref or getattr(ref,'load_status',None)!='VALID': raise RuntimeError('REFERENCE_OUTPUT_INVALID')
            class G: pass
            g=G(); g.resname=ref.residue_name; g.chain=ref.chain; g.resseq=ref.residue_number
            return run_receptor_stage(ctx.receptor,ctx.project,g,allow_bad_res=False)
        return receptor
    if name=='grid':
        def grid():
            ctx.rehydrate_global_outputs(); ref=ctx.metadata.get('reference_output'); rec=ctx.metadata.get('receptor_output')
            if not ref or getattr(ref,'load_status',None)!='VALID' or not rec or getattr(rec,'load_status',None)!='VALID': raise RuntimeError('GLOBAL_PREREQUISITE_INVALID')
            from .intake import hetero_groups
            ref_path=ctx.project/ref.reference_structure; grid_box=box_from_reference(ref_path)
            return run_grid_stage(ctx.project/rec.prepared_receptor_pdbqt,ctx.project/'grid',grid_box,['C','N','O','S'],resolve_executable('autogrid4','AutoGrid4'))
        return grid
    if name=='validation':
        def validation():
            ctx.rehydrate_global_outputs(); ref=ctx.metadata.get('reference_output'); rec=ctx.metadata.get('receptor_output'); grid=ctx.metadata.get('grid_output')
            if not all(x and getattr(x,'load_status',None)=='VALID' for x in (ref,rec,grid)): raise RuntimeError('GLOBAL_PREREQUISITE_INVALID')
            native=ctx.project/ref.reference_structure; info={'resname':ref.residue_name,'chain':ref.chain,'resseq':ref.residue_number}
            native_prepared=prepare_native_ligand(native,ctx.project/'analysis/validation')
            box=type('Grid',(),{'center_x':(grid.center or [0,0,0])[0],'center_y':(grid.center or [0,0,0])[1],'center_z':(grid.center or [0,0,0])[2],'npts_x':(grid.dimensions or [0,0,0])[0],'npts_y':(grid.dimensions or [0,0,0])[1],'npts_z':(grid.dimensions or [0,0,0])[2],'spacing':grid.spacing or 0.375})()
            from .model import DockingParameters
            return run_validation_stage(native,ctx.project/rec.prepared_receptor_pdbqt,ctx.project/'grid',ctx.project/'docking/work'/'native_validation',ctx.project/'analysis/validation',box,DockingParameters(ga_runs=ctx.config.get('validation',{}).get('ga_runs',10)),info)
        return validation
    if name=='finalization':
        def finalization():
            from .finalization import load_finalization_input, snapshot_evidence
            from .dashboard import write_results, write_dashboard
            from .compliance import evaluate
            run_finalization_stage(ctx.project,None,write_results_fn=write_results,write_dashboard_fn=write_dashboard,compliance_fn=evaluate)
            return {'results_csv':ctx.project/'results.csv','evidence_json':ctx.project/'analysis/evidence.json','report_html':ctx.project/'report.html','project_status':snapshot_evidence(load_finalization_input(ctx.project)).get('project_status','PENDING')}
        return finalization
    return None

def _persist_global_result(ctx: ProductionContext, name: str, result: StageResult) -> None:
    """Persist a global seam result and link it from project state."""
    from .stage_outputs import (ReferenceStageOutput, ReceptorStageOutput, GridStageOutput,
        ValidationStageOutput, artifact, save_reference_output, save_receptor_output,
        save_grid_output, save_validation_output)
    from .state import update_state
    out=result.outputs or {}; p=ctx.project
    if name=='reference':
        ref=out.get('reference_pdb') or out.get('reference_structure')
        if ref is None: return
        ref=Path(ref); selection_meta=out.get('selection_metadata'); extra=[artifact(p,ref)]
        for candidate in (out.get('candidate_csv'), out.get('selection_json')):
            if candidate: extra.append(artifact(p,candidate))
        value=ReferenceStageOutput(out.get('selected_ligand',out.get('reference_ligand','')),out.get('resname',''),out.get('chain',''),str(out.get('resseq','')),str(ref),out.get('reference_chemistry'),out.get('ccd_id'),selection_meta,out.get('warnings',[]),extra,out.get('source_metadata',{}),out.get('input_hashes',{}))
        meta=save_reference_output(p,value); update_state(p,'REFERENCE',result.status,inputs=out.get('input_hashes',{}),outputs={str(meta):artifact(p,meta).sha256},result_metadata=str(meta.relative_to(p)))
    elif name=='receptor':
        pdb=Path(out.get('clean') or out.get('prepared_receptor_pdb')); pdbqt=Path(out.get('receptor_pdbqt') or out.get('prepared_receptor_pdbqt'))
        value=ReceptorStageOutput(str(pdb),str(pdbqt),out.get('completeness_report'),out.get('fallback_metadata'),out.get('retained_chains',[]),out.get('warnings',[]),[artifact(p,pdb),artifact(p,pdbqt)])
        meta=save_receptor_output(p,value); update_state(p,'RECEPTOR',result.status,outputs={str(meta):artifact(p,meta).sha256},result_metadata=str(meta.relative_to(p)))
    elif name=='grid':
        maps=[str(x) for x in out.get('maps',[])]; value=GridStageOutput(str(out.get('gpf','')),str(out.get('glg','')),maps,out.get('center',[]),out.get('dimensions',[]),out.get('spacing'),out.get('warnings',[]),[artifact(p,out[k]) for k in ('gpf','glg') if out.get(k)] + [artifact(p,x) for x in maps])
        meta=save_grid_output(p,value); update_state(p,'GRID',result.status,outputs={str(meta):artifact(p,meta).sha256},result_metadata=str(meta.relative_to(p)))
    elif name=='validation':
        native_dlg=out.get('native_dlg'); value=ValidationStageOutput(native_dlg,out.get('minimum_rmsd',out.get('rmsd_A')),out.get('median_rmsd'),out.get('lowest_energy_rmsd'),out.get('rmsd_threshold'),out.get('scientific_status',out.get('status','PASS')),out.get('execution_status','SUCCESS'),out.get('native_like_pose_count'),[artifact(p,native_dlg)] if native_dlg else [],out.get('warnings',[]))
        meta=save_validation_output(p,value); update_state(p,'VALIDATION',result.status,outputs={str(meta):artifact(p,meta).sha256},result_metadata=str(meta.relative_to(p)))
    elif name=='finalization':
        from .finalization_outputs import FinalizationStageOutput, save_finalization_output
        paths=[Path(out['results_csv']),Path(out['evidence_json']),Path(out['report_html'])]
        arts=[{'path':str(x.relative_to(p)),'sha256':__import__('hashlib').sha256(x.read_bytes()).hexdigest()} for x in paths]
        value=FinalizationStageOutput(str(paths[0].relative_to(p)),str(paths[1].relative_to(p)),str(paths[2].relative_to(p)),out.get('project_status','PENDING'),[],arts)
        meta=save_finalization_output(p,value); update_state(p,'FINALIZATION',result.status,outputs={str(meta):artifact(p,meta).sha256},result_metadata=str(meta.relative_to(p)))

def production_context_from_project(project: Path, config: dict) -> ProductionContext:
    """Build a resume context from recorded project artifacts."""
    import json
    data=json.loads((project/'project.json').read_text(encoding='utf-8'))
    from .model import DockingParameters
    params=DockingParameters(**{k:v for k,v in data.get('docking',{}).items() if k in DockingParameters.__dataclass_fields__})
    ctx=ProductionContext(project, config, Path(data.get('receptor_source','')),
                          Path(data.get('ligand_source','')))
    ctx.rehydrate_global_outputs()
    def prepared(cid):
        return {'prepared': project/'ligands/prepared'/f'{cid}.pdbqt'}
    def dock(cid):
        lig=project/'ligands/prepared'/f'{cid}.pdbqt'
        return run_compound_docking(lig, 'receptor', project/'docking/work', project/'docking/dpf', project/'docking/dlg', [], params, 'autodock4')
    def pose(cid): return run_pose_extraction_stage(project/'docking/dlg'/f'{cid}.dlg', project/'analysis')
    def cluster(cid): return run_cluster_stage(project/'docking/dlg'/f'{cid}.dlg', project/'analysis')
    def plip(cid):
        pose=run_pose_extraction_stage(project/'docking/dlg'/f'{cid}.dlg', project/'analysis')
        selected=next((p['complex'] for p in pose.poses if p.get('model')==pose.selected_model), None)
        if selected is None: raise RuntimeError('E701 PLIP_FAILED')
        return run_plip_stage(selected, project/'analysis'/f'plip_{cid}', project/'analysis'/f'interactions_plip_{cid}.csv', compound_id=cid, pose_id=str(pose.selected_model))
    def figures(cid): return run_figures_stage(project/'analysis', project/'analysis/figures')
    ctx.callbacks.update({'prepare':prepared,'dock':dock,'pose':pose,'cluster':cluster,'plip':plip,'figures':figures})
    return ctx
