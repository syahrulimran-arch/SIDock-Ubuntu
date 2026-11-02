"""Canonical production executable resolution."""
from __future__ import annotations
import importlib.util, os, shutil, sys, sysconfig
from pathlib import Path

def resolve_meeko_receptor_preparer() -> dict:
    module = importlib.util.find_spec('meeko')
    version = None
    if module:
        try:
            from importlib.metadata import version as package_version
            version = package_version('meeko')
        except Exception: pass
    candidates=[]
    found=shutil.which('mk_prepare_receptor.py')
    if found: candidates.append(Path(found))
    candidates.append(Path(sysconfig.get_path('scripts'))/'mk_prepare_receptor.py')
    for path in candidates:
        if path.is_file() and os.access(path, os.X_OK):
            return {'available':True,'path':str(path.resolve()),'source':'PATH' if found and Path(found)==path else 'active-python-scripts','version':version,'module_available':bool(module),'message':'Meeko receptor preparer available.'}
    return {'available':False,'path':None,'source':None,'version':version,'module_available':bool(module),'message':('Meeko Python package is installed, but mk_prepare_receptor.py required for receptor preparation was not found.' if module else 'Meeko Python package and mk_prepare_receptor.py were not found.')}

def resolve_meeko_script(name: str) -> str | None:
    found=shutil.which(name)
    candidates=[Path(found)] if found else []
    candidates.append(Path(sysconfig.get_path('scripts'))/name)
    for path in candidates:
        if path.is_file() and os.access(path,os.X_OK): return str(path.resolve())
    return None

def resolve_plip_executable() -> dict:
    try:
        from importlib.metadata import version as package_version
        version = package_version('plip')
    except Exception:
        version = None
    found = shutil.which('plip')
    candidates = [Path(found)] if found else []
    candidates.append(Path(sysconfig.get_path('scripts')) / 'plip')
    for path in candidates:
        if path.is_file() and os.access(path, os.X_OK):
            return {'available': True, 'path': str(path.resolve()), 'source': 'PATH' if found else 'active-python-scripts', 'version': version, 'message': 'PLIP executable available.'}
    module = bool(importlib.util.find_spec('plip'))
    return {'available': False, 'path': None, 'source': None, 'version': version,
            'message': 'PLIP Python package is installed, but the plip executable was not found.' if module else 'PLIP Python package and executable were not found.'}
