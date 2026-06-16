# Usage Guide

## Setup

Create a virtual environment and install Python dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Make sure `TMalign` is installed and visible on `PATH`:

```bash
which TMalign
```

Export the ESM API key:

```bash
export ESM_API_KEY='YOUR_KEY_HERE'
```

## Input Tables

### Sequence table

Required columns:

- `sequence_id`
- `target_name`
- `sequence`

Example:

```csv
sequence_id,target_name,sequence
variant_0001,tem,HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERF...
variant_0002,aph,MIEQDGLHAGSPAAWVERLFGYDWAQQTIGCSDAAVFRLSA...
```

### Target table

Required columns:

- `target_name`
- `template_path`

Optional columns:

- `template_chain`
- `description`

Example:

```csv
target_name,template_path,template_chain,description
tem,/abs/path/to/TEM_pdb.ent,A,TEM reference structure
aph,/abs/path/to/APH_pdb.ent,A,APH reference structure
```

If `template_chain` is supplied, the pipeline extracts that chain into `prepared_templates/` before running TM-align.

## Included Templates

Repository examples:

- `examples/full_run_template/targets.csv`
- `examples/full_run_template/sequences.csv`
- `examples/full_run_template/sequences.example.csv`

The `sequences.csv` file is meant to be edited for your own run. The `sequences.example.csv` file is a valid concrete example.

## Running The Pipeline

### Smoke test

```bash
./run_smoke_test.sh
```

This:

- builds a small input set from local TEM and APH reference PDBs
- runs two unique tasks
- writes outputs under `outputs/smoke_test` unless you override the destination

### Direct CLI usage

```bash
python3 run_sequence_tm_pipeline.py run \
  --input path/to/sequences.csv \
  --targets path/to/targets.csv \
  --outdir path/to/output_dir
```

### Full wrapper

```bash
./run_full_pipeline.sh path/to/sequences.csv path/to/targets.csv path/to/output_dir
```

## Important Flags

### `--prepare-only`

Validate inputs and write placeholder outputs without contacting the ESM API or running TM-align.

Use this first if you want to sanity-check:

- column names
- target paths
- chain extraction
- normalized inputs

### `--limit`

Restrict the run to the first `N` unique tasks after normalization.

Useful for:

- small dry runs
- API testing
- queue debugging

### `--max-concurrent-requests`

Controls concurrent fold API calls.

Increase this for throughput if the API and your quota tolerate it. Keep it low if you are troubleshooting failures or rate limiting.

### `--max-concurrent-tmalign`

Controls concurrent local TM-align processes.

### `--allow-duplicate-sequence-ids`

By default, duplicate `sequence_id` values raise an error. This is intentional, because duplicated row identifiers make downstream interpretation messy.

Use this flag only if duplicate IDs are truly deliberate.

### `--retry-failed`

By default, checkpoint entries with status:

- `scored`
- `fold_failed`
- `tmalign_failed`

are treated as complete.

If you want to retry failed tasks without editing `checkpoint.jsonl`, rerun the same command with:

```bash
--retry-failed
```

In that mode, only previously `scored` tasks are considered complete.

## Output Files

### `normalized_input.csv`

Validated and normalized input rows with:

- resolved template paths
- normalized sequences
- computed `fold_task_id`

### `unique_results.csv`

One row per unique fold-and-align task.

### `sequence_results.csv`

One row per original input row. This is the main analysis table if the same unique task appears in multiple input rows.

### `checkpoint.jsonl`

Append-only progress log used for resuming interrupted runs.

### `summary.json`

Top-level counts plus per-target summary statistics for the primary `tm_score`.

### `attempt_log.jsonl`

Retry and exception log for fold API calls.

### `folded_pdbs/`

Cached predicted structures keyed by `fold_task_id`.

### `prepared_templates/`

Chain-filtered WT/reference PDBs when `template_chain` is provided.

## Resuming After Interruption

The repository is designed so you can rerun the same command safely.

Typical cases:

- terminal closed
- laptop slept
- API timeout
- network interruption

What happens on rerun:

- `checkpoint.jsonl` is loaded
- completed tasks are skipped
- cached PDBs are reused
- aggregate CSV and JSON outputs are refreshed

## Troubleshooting

### `ESM_API_KEY environment variable not set`

Export the key in the current shell:

```bash
export ESM_API_KEY='YOUR_KEY_HERE'
```

### `TMalign not found in PATH`

Install TM-align and confirm:

```bash
which TMalign
```

### DNS or connection failures during folding

Check:

- network access
- API availability
- local firewall or sandbox restrictions

Look at:

- `attempt_log.jsonl`

### Duplicate `sequence_id` error

Either:

- fix the input so each row has a unique `sequence_id`
- or rerun with `--allow-duplicate-sequence-ids` if the duplication is intentional

### Wrong target or wrong chain

Check:

- `target_name` values in the sequence table
- `template_path` values in the target table
- `template_chain` values
- generated files in `prepared_templates/`

## Provenance

This repository was generalized from an experiment-specific runner kept as:

- `run_final_tm_scores.py`

The current pipeline keeps the working fold and TM-align mechanics while removing experiment-specific parsing and scoring logic.
