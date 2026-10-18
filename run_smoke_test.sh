#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${1:-$ROOT/outputs/smoke_test}"
INPUT_DIR="$OUTDIR/inputs"
MODE="${2:---prepare-only}"
if [[ "$MODE" != "--prepare-only" && "$MODE" != "--live" ]]; then
  echo "usage: $0 [OUTDIR] [--prepare-only|--live]" >&2
  exit 1
fi
EXTRA_ARGS=(--progress-log-every 1)
if [[ "$MODE" == "--prepare-only" ]]; then
  EXTRA_ARGS+=(--prepare-only)
fi

python3 "$ROOT/run_sequence_tm_pipeline.py" make-smoke-test --outdir "$INPUT_DIR"
python3 "$ROOT/run_sequence_tm_pipeline.py" run \
  --input "$INPUT_DIR/sequences.csv" \
  --targets "$INPUT_DIR/targets.csv" \
  --outdir "$OUTDIR/results" \
  --limit 2 \
  --max-concurrent-requests 1 \
  --max-concurrent-tmalign 1 \
  "${EXTRA_ARGS[@]}"
