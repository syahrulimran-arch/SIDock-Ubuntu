from pathlib import Path

from sidock_ubuntu.dependencies import resolve_meeko_receptor_preparer


def test_meeko_resolver_finds_active_python_script(monkeypatch, tmp_path):
    script=tmp_path/'mk_prepare_receptor.py'; script.write_text('#!/usr/bin/env python\n'); script.chmod(0o755)
    monkeypatch.setattr('sidock_ubuntu.dependencies.shutil.which', lambda name: None)
    monkeypatch.setattr('sidock_ubuntu.dependencies.sysconfig.get_path', lambda name: str(tmp_path))
    result=resolve_meeko_receptor_preparer()
    assert result['available'] and result['path']==str(script)


def test_meeko_resolver_distinguishes_missing_script(monkeypatch, tmp_path):
    monkeypatch.setattr('sidock_ubuntu.dependencies.shutil.which', lambda name: None)
    monkeypatch.setattr('sidock_ubuntu.dependencies.sysconfig.get_path', lambda name: str(tmp_path))
    result=resolve_meeko_receptor_preparer()
    assert not result['available'] and 'mk_prepare_receptor.py' in result['message']


def test_doctor_reports_meeko_module_and_receptor_preparer_separately():
    from sidock_ubuntu.cli import doctor_status
    data=doctor_status(); names={x['name'] for x in data['dependencies']}
    assert 'Meeko Python module' in names
    assert 'Meeko receptor preparer (mk_prepare_receptor.py)' in names
