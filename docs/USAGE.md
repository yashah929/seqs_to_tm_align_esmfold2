# Usage and operations

Start with the [README quick start](../README.md#quick-start-offline-first). Run commands from the repository root unless paths are adjusted.

## Tables and validation

`--input` accepts CSV or TSV with `sequence_id,target_name,sequence`. IDs are strings, including numeric-looking IDs with leading zeros. IDs and target names must be nonblank; unknown targets are errors. Repeated IDs require `--allow-duplicate-sequence-ids`. Distinct IDs with the same target, sequence, and model naturally share one task.

The sequence normalizer removes whitespace, uppercases letters, and accepts only `ACDEFGHIKLMNPQRSTVWY`. Stops, gaps, ambiguous residues, missing CSV cells, and lengths outside 3–9999 are rejected. The PDB workflow imposes this length bound; the service may have stricter model limits. Extra input columns are not retained: keep strategy/assay metadata in a separate table keyed by ID.

`--targets` needs `target_name,template_path`, optionally `template_chain,description`. Target names must match `[A-Za-z0-9][A-Za-z0-9_.-]*`. Reference paths resolve **relative to the target table**, not the working directory. Examples use `../../TEM_pdb.ent` and `../../APH_pdb.ent`. A chain is a single PDB character. Leave it blank only for a reference with exactly one ATOM chain (a blank chain ID is supported in that case).

PDB fixed-column files are supported; mmCIF is not read by the preparation code. Multiple structural models are rejected; extract a single model deliberately first. Every reference is prepared into one Cα-only chain. Only standard amino-acid `ATOM` Cα records with positive occupancy are retained; ligands, water, and `HETATM` modified residues are excluded. Alternate locations use highest occupancy, then blank, then A, then file order for ties. Duplicate residue/conformer keys, nonfinite Cα coordinates, unsupported ATOM residue names, and fewer than three usable Cα atoms are errors. Missing residues are not rebuilt; PDB residue numbering and insertion codes are retained in the reference. See [reference coverage](SCIENTIFIC_CONTEXT.md).

## Commands and wrappers

```bash
python3 run_sequence_tm_pipeline.py run --help
python3 run_sequence_tm_pipeline.py make-smoke-test --help
./run_smoke_test.sh outputs/offline_smoke
```

The smoke wrapper generates inputs and validates the first two unique TEM tasks offline by default. `./run_smoke_test.sh outputs/live_smoke_pair --live` performs two live predictions (TEM baseline and teaching mutant) and requires a Biohub token and `TMalign`. It does not test APH unless you run the full four-row table. The general Bash wrapper accepts `INPUT TARGETS OUTDIR [extra CLI options...]` and passes options through unchanged. Both work from other working directories; relative input arguments are relative to your current shell directory.

Regenerate all four example inputs anywhere offline:

```bash
python3 run_sequence_tm_pipeline.py make-smoke-test --outdir outputs/generated_inputs
```

This extracts chain-A `SEQRES` from the bundled references and adds H1A/M1A teaching variants. To use other structures, `--tem-template` and `--aph-template` must have suitable chain-A `SEQRES`; names remain TEM/APH. For general targets, write tables directly.

`--prepare-only` validates inputs/references and writes `pending` output rows. It does not contact the API or execute TM-align. Use a separate preparation directory: preparation is rejected beside an existing checkpoint to avoid overwriting scored aggregate tables. `--limit N` selects the first N unique tasks and all input rows associated with them; it is not a limit on rows or a random sample.

## Live service and retries

The general CLI always requires `ESM_API_KEY` and `TMalign` for `run` without `--prepare-only`, including resumed runs. Export a Biohub token in your shell; never put it in an input table or commit it. The request uses Bearer authentication, JSON `sequence` and `model`, and the endpoint documented in [service verification](SCIENTIFIC_CONTEXT.md#folding-service-and-coordinate-contract). The default model is `esmfold2-fast-2026-05`. `--model` is passed through; no local availability check or backend fallback occurs. A reported response model that differs from the requested identifier is rejected.

`--max-retries` retains its historical name but means **maximum total attempts**, including the first attempt (default 8). HTTP 429, 500, 502, 503, 504 and network/timeouts are retried with `--retry-sleep-seconds × attempt` delay (default 10 seconds). HTTP 403 receives a 60-second cooldown and another attempt; other HTTP errors end that task. Malformed success payloads fail conversion without repeated requests. The implementation does not honor `Retry-After`. Inspect `attempt_log.jsonl` and each result's `error`; a failure is never assigned a score of zero.

`--request-timeout` defaults to 120 seconds. Request and TM-align concurrency each default to 4; use 1 while checking service access. Retry delays occupy a request slot. `--cooldown-seconds` delays completion of each task; it is not a global rate limiter. Positive concurrency, timeout, attempts, and limits are required. Delays and progress intervals must be nonnegative; zero progress intervals disable periodic logging/saving. Completed tasks always append a checkpoint.

## Output field guide

| File | Contents |
| --- | --- |
| `normalized_input.csv` | Input row index (1-based), sequence ID, target, normalized sequence and input length, model, task ID, prepared reference path/chain, source reference path/hash, scoring signature |
| `sequence_results.csv` | One row per normalized input row plus status, metrics, timings, prediction path/hash, error, and completion time |
| `unique_results.csv` | One row per distinct target + sequence + model, deduplicated from the checkpoint for the current scoring configuration |
| `summary.json` | Input-row and unique-task counts; row status counts; per-target row-weighted TM statistics over scored rows |
| `checkpoint.jsonl` | Append-only completion records, including failures and provenance; authoritative resume log |
| `attempt_log.jsonl` | HTTP/network retry details; conversion failures appear in the checkpoint's `error` |
| `folded_pdbs/` | Cα-only prediction PDBs with JSON provenance sidecars; override directory with `--folded-dir` |
| `prepared_templates/` | Single-chain Cα reference files actually passed to TM-align |

Result status is `pending`, `scored`, `fold_failed`, or `tmalign_failed`. Failed/pending metric cells are empty. Metric columns are defined in the [README](../README.md#reading-the-scores). `fold_attempts=0` means verified prediction-cache reuse. `fold_elapsed_s` includes API retries/waiting or cache validation; `tmalign_elapsed_s` includes waiting for an alignment slot. `completed_unix` is a Unix timestamp. API confidence outputs are not scored or retained as confidence tracks; do not interpret PDB B factors as calibrated confidence.

## Cache and checkpoint provenance

Rerun the same command and `--outdir` to resume. Do not run two processes writing the same outdir/cache concurrently. Paths in generated outputs are absolute; after moving a run, prepare/rescore against the new paths rather than relying on those saved paths.

A `fold_task_id` hashes target name, requested model, and normalized sequence. Prediction reuse also requires matching endpoint/conversion policy metadata and the exact PDB SHA-256. PDBs without sidecars (including legacy caches) are unverified and will be replaced only after a successful fold. Preserve copies before upgrading if you need those historical artifacts. Sidecar and PDB writes use atomic file replacement; an interrupted mismatched pair is not reused. No network call is made by offline preparation to verify a cache.

Scoring signatures hash the original reference file contents, prepared coordinates, selected chain, reference-selection policy, TM-align executable, argument order, conversion policy, and endpoint. Changing any of these invalidates checkpointed scores. Changing sequence/model changes the task ID. Valid cached predictions survive reference/chain/TM-align changes and can be re-aligned without another fold. Even header-only reference changes conservatively invalidate scores. Switching a remote deployment under an unchanged model identifier cannot be detected: archive cache metadata and use a fresh cache/outdir if you need a new prediction realization.

Only current matching checkpoint records are included in results. A scored record also requires its prediction file to exist with its recorded hash. Records without new provenance fields are explicitly warned about and **not treated as verified**. Existing log entries are retained; new verified results append to the log. A malformed/truncated JSON line stops with its location so the original checkpoint can be preserved and repaired deliberately.

By default, matching `scored`, `fold_failed`, and `tmalign_failed` records are complete. `--retry-failed` retries failures, reusing valid predictions after alignment failures. Latest matching task completion wins. CSV/summary snapshots refresh every `--progress-save-every` tasks and at the end; after an abrupt stop they may lag behind the checkpoint until rerun.

## Legacy runner

`run_final_tm_scores.py` remains at its original entry point to preserve provenance. It infers TEM/APH targets from fixed lengths and expects external experiment `configs/`, `outputs/`, MSA, and reference directories that are not bundled. Its scoring, cache, and service conventions are historical and do not inherit the general CLI's new validation/provenance guarantees. It is not the portable example runner. No experiment-specific data or predictions were imported into this cleanup.
