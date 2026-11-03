from sidock_ubuntu.structure_context import compare_validation_receptor

def _p(path, mutation=None, shift=0.0):
    lines=[]
    for i,name in enumerate(('ALA','GLY','SER','VAL'),1):
        if mutation and i==mutation: name='ASP'
        lines += [f'ATOM  {i:5d}  CA  {name:>3} A{i:4d}    {1.0*i+shift:8.3f}{2.0:8.3f}{3.0:8.3f}  1.00 20.00           C']
    path.write_text('\n'.join(lines)+'\n')

def test_stripped_derivative_equivalence(tmp_path):
    a,b=tmp_path/'production.pdb',tmp_path/'validation.pdb'; _p(a); _p(b)
    c=compare_validation_receptor(a,b)
    assert c.classification=='PREPARED_DERIVATIVE_OF_SAME_STRUCTURE' and c.mapping_coverage==1.0 and c.sequence_identity==1.0

def test_different_structure_and_binding_mutation_rejected(tmp_path):
    a,b=tmp_path/'production.pdb',tmp_path/'validation.pdb'; _p(a); _p(b,shift=3.0)
    assert compare_validation_receptor(a,b).classification=='RELATED_BUT_DIFFERENT_STRUCTURE'
    _p(b,mutation=2); assert compare_validation_receptor(a,b,binding_site={('A','2')}).classification=='INCOMPATIBLE'
