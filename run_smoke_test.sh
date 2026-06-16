#!/bin/zsh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
INPUT_DIR="$ROOT/examples/smoke_test_inputs"
OUTDIR="${1:-$ROOT/outputs/smoke_test}"

python3 "$ROOT/run_sequence_tm_pipeline.py" make-smoke-test --outdir "$INPUT_DIR"
python3 "$ROOT/run_sequence_tm_pipeline.py" run \
  --input "$INPUT_DIR/sequences.csv" \
  --targets "$INPUT_DIR/targets.csv" \
  --outdir "$OUTDIR" \
  --limit 2 \
  --max-concurrent-requests 1 \
  --max-concurrent-tmalign 1 \
  --progress-log-every 1
