#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

echo "SIDock Ubuntu installer"
for cmd in python3; do command -v "$cmd" >/dev/null || { echo "Missing $cmd"; exit 1; }; done
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip wheel setuptools
python -m pip install -e .

cat <<EOF

Python environment installed.
External executables still required on PATH:
  autogrid4
  autodock4
Optional for publication pose figures:
  pymol

Activate with:
  source .venv/bin/activate
Check with:
  sidock-ubuntu --help
EOF
