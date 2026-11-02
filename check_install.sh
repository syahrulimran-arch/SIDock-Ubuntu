#!/usr/bin/env bash
set -u
fail=0
for cmd in python3 autogrid4 autodock4; do
  if command -v "$cmd" >/dev/null; then echo "OK   $cmd -> $(command -v "$cmd")"; else echo "MISS $cmd"; fail=1; fi
done
if command -v pymol >/dev/null; then echo "OK   pymol -> $(command -v pymol)"; else echo "INFO pymol not installed; PML scripts will still be generated"; fi
python3 - <<'PY' || fail=1
mods=['rdkit','meeko','scipy','gemmi','matplotlib']
for m in mods:
    try: __import__(m); print('OK  ',m)
    except Exception as e: print('MISS',m,e); raise
PY
exit $fail
