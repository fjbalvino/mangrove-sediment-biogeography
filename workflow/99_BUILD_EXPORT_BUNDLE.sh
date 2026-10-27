#!/usr/bin/env bash
# Export source only. Inputs, databases, results and figures are external assets.
set -euo pipefail
PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec "${PYTHON_BIN:-python3}" "$PACKAGE_ROOT/tools/build_code_release.py" "$@"
