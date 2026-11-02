import json
from pathlib import Path
import pytest

from sidock_ubuntu.reference import (
    canonical_chemical_graph, chemical_graph_id, confidence_for_candidates,
    extract_mmcif_quality, group_equivalent_instances, map_instance_to_ccd,
    select_reference_instance, write_reference_selection, load_local_ccd,
    load_local_mmcif, reference_config_fingerprint, rank_instances,
)

def graph(order=1):
    return {'atoms':[{'id':'C1','element':'C'},{'id':'O1','element':'O'}], 'bonds':[('C1','O1',order)]}

def candidate(key, gid, **kw):
    return {'key':key,'ligand_code':'LIG','chemical_group_id':gid,'mapping_fraction':kw.pop('mapping_fraction',1.0), **kw}

def test_chemical_graph_distinguishes_bond_topology():
    assert chemical_graph_id(graph(1)) != chemical_graph_id(graph(2))
    assert canonical_chemical_graph(graph(1)) == canonical_chemical_graph(graph(1))

def test_same_ccd_instances_group_across_chains():
    rows=[candidate('LIG:A:1','CCD-X'),candidate('LIG:B:1','CCD-X'),candidate('LIG:C:1','CCD-X')]
    assert list(group_equivalent_instances(rows)) == ['CCD-X']
    assert len(group_equivalent_instances(rows)['CCD-X']) == 3

def test_mapping_quality_is_explicit():
    out=map_instance_to_ccd({'key':'A','atoms':['C1','O1']},{'id':'X','atoms':['C1','O1','N1'],'bonds':[]})
    assert out['mapped_heavy_atoms']==2 and out['total_ccd_heavy_atoms']==3
    assert out['mapping_fraction']==pytest.approx(2/3) and out['missing_heavy_atoms']==1

def test_complete_mapping_ranks_before_incomplete_and_altloc():
    rows=[candidate('B','X',mapping_fraction=.8,altloc_present=True),candidate('A','X',mapping_fraction=1.0,occupancy_summary=.7)]
    assert select_reference_instance(rows)['selected']['key']=='A'

def test_distinct_groups_are_low_confidence_and_require_override():
    rows=[candidate('A','X'),candidate('B','Y')]
    assert confidence_for_candidates(rows)=='LOW'
    with pytest.raises(RuntimeError,match='E101'): select_reference_instance(rows)
    assert select_reference_instance(rows,'B')['selected']['key']=='B'
    assert select_reference_instance(rows,'B')['selected']['selection_method']=='EXPLICIT_OVERRIDE'

def test_invalid_override_is_controlled():
    with pytest.raises(RuntimeError,match='E101'): select_reference_instance([candidate('A','X')],'missing')

def test_selection_csv_and_json_are_persisted(tmp_path):
    rows=[candidate('LIG:A:1','X',occupancy_summary=.8),candidate('LIG:B:1','X',occupancy_summary=.6)]
    result=write_reference_selection(tmp_path,rows,input_hashes={'pdb':'abc'},provenance={'source':'local'})
    assert Path(result['candidate_csv']).exists() and Path(result['selection_json']).exists()
    data=json.loads(Path(result['selection_json']).read_text())
    assert data['selection_method']=='AUTOMATIC' and data['input_hashes']['pdb']=='abc'

def test_pdb_only_quality_fallback_has_no_fabricated_mmcif_values(tmp_path):
    missing=tmp_path/'missing.cif'
    assert extract_mmcif_quality(missing)=={}

def test_local_ccd_loader_records_hash_and_chemistry(tmp_path):
    ccd=tmp_path/'P43.json'; ccd.write_text(json.dumps({'component_id':'P43','atoms':[{'id':'C1','element':'C'}],'bonds':[]}))
    result=load_local_ccd(ccd,'P43')
    assert result['available'] and result['parse_status']=='PASS' and result['source_hash'] and result['chemistry']['component_id']=='P43'

def test_invalid_local_ccd_is_controlled(tmp_path):
    ccd=tmp_path/'bad.json'; ccd.write_text('{not json')
    result=load_local_ccd(ccd,'P43')
    assert result['available'] is False and result['parse_status']=='INVALID' and result['warnings']

def test_missing_local_sources_are_nonfatal_and_provenanced(tmp_path):
    ccd=load_local_ccd(tmp_path/'none.json','P43'); cif=load_local_mmcif(tmp_path/'none.cif')
    assert ccd['parse_status']=='MISSING' and cif['parse_status']=='MISSING' and cif['available'] is False

def test_reference_fingerprint_excludes_unrelated_settings():
    a=reference_config_fingerprint({'reference':{'ligand':'A'},'docking':{'ga_runs':1},'interaction':{'x':1}})
    b=reference_config_fingerprint({'reference':{'ligand':'A'},'docking':{'ga_runs':999},'interaction':{'x':2}})
    assert a==b

def test_reference_fingerprint_changes_override():
    a=reference_config_fingerprint({'reference':{'ligand':'A'}}); b=reference_config_fingerprint({'reference':{'ligand':'B'}})
    assert a!=b

def test_local_mmcif_fields_are_extracted_when_present(tmp_path):
    cif=tmp_path/'x.cif'
    cif.write_text('''data_x\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.label_alt_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.pdbx_PDB_model_num\nHETATM 1 C C1 P43 A 1 A 1263 ? . 1 2 3 0.75 18.0 1\n''')
    result=load_local_mmcif(cif)
    assert result['available'] and result['parse_status']=='PASS'
    assert next(iter(result['instances'].values()))['occupancy_summary']==pytest.approx(.75)

def test_mapping_outranks_better_rscc():
    rows=[candidate('A','X',mapping_fraction=1,rscc=.80),candidate('B','X',mapping_fraction=.8,rscc=.99)]
    assert rank_instances(rows)[0]['key']=='A'

def test_altloc_outranks_lower_priority_metrics():
    rows=[candidate('A','X',altloc_present=False,rscc=.80),candidate('B','X',altloc_present=True,rscc=.99)]
    assert rank_instances(rows)[0]['key']=='A'

def test_occupancy_is_deterministic_priority():
    rows=[candidate('A','X',occupancy_summary=1.0),candidate('B','X',occupancy_summary=.5)]
    assert rank_instances(rows)[0]['key']=='A'

def test_rscc_then_rsr_ranking():
    assert rank_instances([candidate('A','X',rscc=.94),candidate('B','X',rscc=.88)])[0]['key']=='A'
    assert rank_instances([candidate('A','X',rsr=.12),candidate('B','X',rsr=.21)])[0]['key']=='A'

def test_missing_quality_metrics_remain_rankable():
    rows=[candidate('A','X'),candidate('B','X',rscc=None,rsr=None)]
    assert rank_instances(rows)[0]['key']=='A'

def test_mixed_quality_metrics_order_is_repeatable():
    rows=[candidate('A','X',rscc=.9,rsr=.2),candidate('B','X',rscc=.9),candidate('C','X',rsr=.1),candidate('D','X')]
    assert [[r['key'] for r in rank_instances(rows)] for _ in range(3)]==[['A','B','C','D']]*3

def test_binding_site_completeness_outranks_peripheral_and_tie_break():
    rows=[candidate('A','X',binding_site_incomplete_residue_count=0,peripheral_incomplete_residue_count=1),candidate('B','X',binding_site_incomplete_residue_count=1,peripheral_incomplete_residue_count=0)]
    assert rank_instances(rows)[0]['key']=='A'

def test_contact_environment_is_used_after_local_completeness():
    rows=[candidate('A','X',contact_environment_complete_fraction=1),candidate('B','X',contact_environment_complete_fraction=.5)]
    assert rank_instances(rows)[0]['key']=='A'

def test_remote_incompleteness_does_not_affect_local_rank():
    rows=[candidate('A','X',remote_incomplete_residue_count=99),candidate('B','X',remote_incomplete_residue_count=0)]
    assert rank_instances(rows)[0]['key']=='A'  # stable identifier tie-break

def test_selection_json_records_policy_and_reason(tmp_path):
    result=write_reference_selection(tmp_path,[candidate('A','X',binding_site_incomplete_residue_count=0)])
    data=json.loads(Path(result['selection_json']).read_text())
    assert data['reference_ranking_policy_version']=='1' and data['ranking_policy'] and data['selection_reason']
