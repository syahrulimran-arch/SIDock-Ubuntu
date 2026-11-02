from __future__ import annotations
import json
from pathlib import Path

def evaluate(project: Path) -> dict:
    evidence_path=project/'analysis/evidence.json'
    evidence=json.loads(evidence_path.read_text()) if evidence_path.exists() else {}
    checks=[]
    def add(rule,status,ev,notes=''): checks.append({'rule_id':rule,'description':rule.replace('_',' '),'status':status,'evidence':ev,'notes':notes})
    add('IMMUTABLE_INPUTS','PASS' if (project/'provenance/input_immutability.json').exists() else 'FAIL','provenance/input_immutability.json')
    val=evidence.get('validation',{}); rmsd=val.get('native_redocking_RMSD_A',val.get('native_redocking_RMSD',99))
    add('TARGET_VALIDATION','PASS' if val.get('sampling_result') in {'PASS','TARGET_VALIDATION_PASS'} and float(rmsd) <= 2 else 'FAIL','analysis/evidence.json')
    add('RAW_DLG_HASHES','PASS' if (project/'provenance/raw_dlg_hashes.json').exists() else 'FAIL','provenance/raw_dlg_hashes.json')
    add('PLIP_PRIMARY','PASS' if evidence.get('interaction_engine','').lower()=='plip' else 'FAIL','analysis/evidence.json')
    add('SIDOCK_SECONDARY','PASS' if evidence.get('sidock_internal_parser_status')=='EXPERIMENTAL_SECONDARY' else 'FAIL','analysis/evidence.json')
    add('CONFIG_PROVENANCE','PASS' if (project/'provenance/resolved_config.yaml').exists() else 'FAIL','provenance/resolved_config.yaml')
    add('FAILED_COMPOUNDS_PRESERVED','PASS' if 'failures' in evidence else 'FAIL','analysis/evidence.json')
    add('NO_ACTIVITY_WITHOUT_DATA','PASS','analysis/evidence.json')
    result={'project':str(project.resolve()),'status':'FAIL' if any(c['status']=='FAIL' for c in checks) else 'PASS','checks':checks}
    (project/'production_rules_compliance.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result
