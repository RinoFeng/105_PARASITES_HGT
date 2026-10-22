#!/usr/bin/env bash
set -euo pipefail
example_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
python "$example_dir/family_comparison/run.py"
python "$example_dir/functional_contribution/run.py"
