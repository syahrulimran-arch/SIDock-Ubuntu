from pathlib import Path
from types import SimpleNamespace

from sidock_ubuntu.executors import ProductionContext, _default_global_callback


def test_automatic_reference_selection_uses_validation_structure(tmp_path, monkeypatch):
    production = tmp_path/'input'/'protein.pdb'
    validation = tmp_path/'input'/'validation_holo.pdb'
    production.parent.mkdir(parents=True, exist_ok=True)
    production.write_text('PRODUCTION STRIPPED\n', encoding='utf-8')
    validation.write_text('VALIDATION HOLO\n', encoding='utf-8')

    group = SimpleNamespace(key='GLC:A:1', resname='GLC', chain='A', resseq='1')

    import sidock_ubuntu.intake as intake
    import sidock_ubuntu.executors as executors

    monkeypatch.setattr(intake, 'semantic_reference_groups', lambda path, mmcif=None: [group])
    monkeypatch.setattr(intake, 'hetero_groups', lambda path: [group])
    monkeypatch.setattr(intake, 'reference_candidates', lambda path: [])

    chosen_from = []

    def choose(path):
        chosen_from.append(Path(path))
        return group

    monkeypatch.setattr(intake, 'choose_reference_ligand', choose)

    def fake_reference_stage(pdb, project, selected):
        assert Path(pdb) == validation
        assert selected is group
        ref = project/'receptor'/'reference_ligand.pdb'
        ref.parent.mkdir(parents=True, exist_ok=True)
        ref.write_text('REFERENCE\n', encoding='utf-8')
        return {'reference_pdb': ref}

    monkeypatch.setattr(executors, 'run_reference_stage', fake_reference_stage)

    ctx = ProductionContext(
        project=tmp_path,
        config={},
        receptor=production,
        validation_structure=validation,
    )
    result = _default_global_callback(ctx, 'reference')()

    assert chosen_from == [validation]
    assert result['selected_ligand'] == group.key
    assert result['source_metadata']['structure_used'] == str(validation)
    assert result['input_hashes']['production_receptor'] != result['input_hashes'][str(validation)]
