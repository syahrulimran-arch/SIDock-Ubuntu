import json
import pytest

from sidock_ubuntu.gui_service import attach_project_inputs, create_project


def fixture_inputs(tmp_path):
    """Create the minimal receptor/SDF pair needed by these attachment tests.

    Keep this helper local so the test module does not depend on ``tests`` being
    installed or importable as a Python package.
    """
    pdb = tmp_path/'protein.pdb'
    pdb.write_text('ATOM\n')
    sdf = tmp_path/'compounds.sdf'
    sdf.write_text(
        'minimal\n'
        '  SIDock\n'
        '\n'
        '  1  0  0  0  0  0  0  0  0  0999 V2000\n'
        '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n'
        'M  END\n'
        '$$$$\n'
    )
    return pdb, sdf


def test_existing_project_inputs_are_replaced_and_dependents_staled(tmp_path):
    pdb, sdf = fixture_inputs(tmp_path)
    project=create_project(tmp_path/'workspace','demo',pdb,sdf)
    validation=tmp_path/'holo.pdb'; validation.write_text('HEADER HOLO\nATOM\n')
    mmcif=tmp_path/'structure.cif'; mmcif.write_text('data_x\n')
    ccd=tmp_path/'ccd.json'; ccd.write_text(json.dumps({'components': []}))
    result=attach_project_inputs(project,validation_structure=validation,mmcif=mmcif,ccd=ccd)
    assert 'validation_structure' in result['updated']
    assert (project/'input/validation_structure.pdb').read_text()==validation.read_text()
    assert json.loads((project/'gui_inputs.json').read_text())['ccd'].endswith('input/ccd.json')
    assert json.loads((project/'project_state.json').read_text())['stages']['REFERENCE']['status']=='STALE'


def test_invalid_attachment_does_not_replace_existing_input(tmp_path):
    pdb, sdf = fixture_inputs(tmp_path)
    project=create_project(tmp_path/'workspace','demo',pdb,sdf)
    old=(project/'input/ccd.json')
    bad=tmp_path/'bad.ccd'; bad.write_text('{not-json')
    with pytest.raises(ValueError): attach_project_inputs(project,ccd=bad)
    assert not old.exists()


def test_completed_project_requires_explicit_input_change_confirmation(tmp_path):
    pdb, sdf = fixture_inputs(tmp_path)
    project=create_project(tmp_path/'workspace','demo',pdb,sdf)
    state=json.loads((project/'project_state.json').read_text())
    state['stages']['FINALIZATION']['status']='PASS'
    (project/'project_state.json').write_text(json.dumps(state))
    validation=tmp_path/'holo.pdb'; validation.write_text('HEADER HOLO\nATOM\n')
    with pytest.raises(ValueError, match='CONFIRMATION'):
        attach_project_inputs(project,validation_structure=validation)
