#!/bin/zsh
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "usage: $0 INPUT_CSV TARGETS_CSV OUTDIR [extra pipeline args...]" >&2
  exit 1
fi

ROOT="$(cd "$(dirname "$0")" && pwd)"
INPUT_CSV="$1"
TARGETS_CSV="$2"
OUTDIR="$3"
shift 3

python3 "$ROOT/run_sequence_tm_pipeline.py" run \
  --input "$INPUT_CSV" \
  --targets "$TARGETS_CSV" \
  --outdir "$OUTDIR" \
  "$@"
