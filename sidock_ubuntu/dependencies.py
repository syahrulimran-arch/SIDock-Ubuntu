"""Canonical production executable resolution."""
from __future__ import annotations
import importlib.util, os, shutil, sys, sysconfig
from pathlib import Path

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
