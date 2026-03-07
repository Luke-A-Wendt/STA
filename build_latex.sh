#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v lualatex >/dev/null 2>&1; then
  echo "Error: lualatex not found. Install TeX Live first." >&2
  exit 1
fi

# Use writable TeX cache paths (helps avoid permission/cache issues).
export TEXMFVAR="${TEXMFVAR:-/tmp/texmf-var}"
export TEXMFCONFIG="${TEXMFCONFIG:-/tmp/texmf-config}"
export TEXMFHOME="${TEXMFHOME:-/tmp/texmf-home}"
mkdir -p "$TEXMFVAR" "$TEXMFCONFIG" "$TEXMFHOME"

GENERATOR="${ROOT_DIR}/tools/generate_readme_equations.py"
if [[ -x "$GENERATOR" ]]; then
  echo "=== Regenerating README conversion files ==="
  "$GENERATOR"
fi

VERIFIER="${ROOT_DIR}/tools/verify_sta_results.py"
if [[ -f "$VERIFIER" ]]; then
  echo "=== Running symbolic verification checks ==="
  python3 "$VERIFIER"
fi

DOC_DIR="manuscript"

echo "=== Building ${DOC_DIR} ==="
(
  cd "${ROOT_DIR}/${DOC_DIR}"
  lualatex -interaction=nonstopmode -halt-on-error main.tex
  lualatex -interaction=nonstopmode -halt-on-error main.tex
)

echo "Build complete."
