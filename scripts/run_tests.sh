#!/usr/bin/env bash
# Run the shredtools test suite in the conda env named "shredtools".
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="${SHREDTOOLS_CONDA_ENV:-shredtools}"

cd "$REPO_ROOT"

if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "Conda environment '$ENV_NAME' not found." >&2
  echo "Create it and install shredtools + pytest first." >&2
  exit 1
fi

# Install dev deps if pytest missing (idempotent).
conda run -n "$ENV_NAME" python -c "import pytest" 2>/dev/null || \
  conda run -n "$ENV_NAME" pip install -q pytest

# Reinstall package from source so tests exercise the working tree.
conda run -n "$ENV_NAME" pip install -q -e .

export SHREDTOOLS_TEST_DATA="${SHREDTOOLS_TEST_DATA:-$REPO_ROOT/test_data}"

echo "Using test data: $SHREDTOOLS_TEST_DATA"
echo "Running: pytest $*"
exec conda run -n "$ENV_NAME" pytest "$@"
