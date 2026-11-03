from sidock_ubuntu.structure_intake import classify_structure


def _pdb(path, *, groups, remarks='HEADER    TEST\n'):
    lines=[remarks]
    serial=1
    for record, resname, chain, resseq, atoms in groups:
        for atom, element in atoms:
            lines.append(f'{record:<6}{serial:5d} {atom:^4s} {resname:>3s} {chain}{int(resseq):4d}    {serial:7.3f}{serial+1:8.3f}{serial+2:8.3f}  1.00 20.00           {element:>2s}')
            serial+=1
    path.write_text('\n'.join(lines)+'\n')


def test_holo_carbohydrate_and_ion_are_separated(tmp_path):
    p=tmp_path/'holo.pdb'
    _pdb(p, groups=[('ATOM','ALA','A',1,[('N','N'),('CA','C'),('C','C'),('O','O')]),('HETATM','GLC','A',2,[(f'C{i}','C') for i in range(1,7)]),('HETATM','CA','A',3,[('CA','CA')])])
    result=classify_structure(p)
    assert result.classification=='HOLO_NATIVE'
    assert [x['component'] for x in result.candidates]==['GLC']
    assert [x['component'] for x in result.ions]==['CA']


def test_processed_protein_only_and_ion_only_are_not_native(tmp_path):
    processed=tmp_path/'3AJ7ncbi.pdb'
    _pdb(processed, remarks='REMARK Accelrys Discovery Studio processed receptor\n', groups=[('ATOM','ALA','A',1,[('N','N'),('CA','C'),('C','C'),('O','O'),('H','H')])])
    assert classify_structure(processed).classification=='PROCESSED_STRIPPED'
    ion=tmp_path/'3AJ7.pdb'
    _pdb(ion, groups=[('ATOM','ALA','A',1,[('N','N'),('CA','C'),('C','C'),('O','O')]),('HETATM','CA','A',2,[('CA','CA')])])
    assert classify_structure(ion).classification=='ION_ONLY'


def test_non_hetatm_organic_component_is_candidate(tmp_path):
    p=tmp_path/'converted.pdb'
    _pdb(p, groups=[('ATOM','LIG','A',2,[(f'C{i}','C') for i in range(1,7)])])
    result=classify_structure(p)
    assert result.classification=='HOLO_NATIVE'
    assert result.candidates[0]['component']=='LIG'


def test_water_only_is_apo_no_reference(tmp_path):
    p=tmp_path/'water.pdb'
    _pdb(p, groups=[('ATOM','ALA','A',1,[('N','N'),('CA','C'),('C','C'),('O','O')]),('HETATM','HOH','A',2,[('O','O')])])
    result=classify_structure(p)
    assert result.classification in {'APO_NO_REFERENCE','PROCESSED_STRIPPED'}
    assert result.waters[0]['component']=='HOH'

def test_local_mmcif_metadata_precedes_ambiguous_pdb(tmp_path):
    p=tmp_path/'ambiguous.pdb'; _pdb(p, groups=[('ATOM','GLC','A',2,[(f'C{i}','C') for i in range(1,7)])])
    cif=tmp_path/'structure.cif'
    cif.write_text('''data_3A4A\n\nloop_\n_chem_comp.id\n_chem_comp.type\nGLC\n"non-polymer"\nCA\n"non-polymer"\nHOH\n"non-polymer"\n''')
    result=classify_structure(p, mmcif=cif)
    assert result.classification=='HOLO_NATIVE'
    assert any(x['component']=='GLC' for x in result.candidates)
    assert any(x['component']=='CA' for x in result.ions)
    assert any(x['component']=='HOH' for x in result.waters)
