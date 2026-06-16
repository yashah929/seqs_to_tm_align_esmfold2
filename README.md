# Sequence-to-TM-Score Pipeline

This repository contains a sequence-to-structure comparison pipeline for protein variant sets tied to known reference structures.

Given a table of amino-acid sequences and a mapping from each sequence to its reference structure, the pipeline:

1. folds each unique sequence through the ESM API
2. converts the fold response into a PDB when needed
3. runs `TMalign` against the appropriate WT/reference structure
4. records a WT-normalized TM-score and supporting alignment metadata
5. checkpoints progress so large runs can be resumed safely

The main entry point is `run_sequence_tm_pipeline.py`.

## Repository Summary

Suggested GitHub repository description:

`General-purpose ESM API to TM-align pipeline for scoring protein sequences against WT reference structures.`

Suggested short tagline:

`Fold sequences with ESM, align to WT with TM-align, and checkpoint the results.`

## Scientific Aim

In many design and variant-screening settings, sequence-level scores are not enough. The relevant question is structural: does a candidate sequence still adopt the wild-type fold, or has it moved away from the reference state?

The pipeline addresses that question in batch form. It folds each candidate through the ESM API, aligns the predicted structure to the correct WT/reference structure with TM-align, and records a structural similarity score suitable for downstream analysis.

The scientific rationale and scoring conventions are described in `docs/SCIENTIFIC_CONTEXT.md`.

## What The Pipeline Measures

The primary output is `tm_score`, defined here as the TM-align score normalized by the WT/reference structure length.

This choice is intentional. In the present workflow the reference structure is the fixed object of interest. The quantity of interest is similarity to that WT/reference structure, not self-normalized agreement on the candidate's own length scale.

The pipeline also records:

- `tm_score_query`: TM-align score normalized by the predicted-query length
- `aligned_length`
- `rmsd`
- `seq_identity_aligned`
- `predicted_length`
- `wt_length`

More detail is in `docs/SCIENTIFIC_CONTEXT.md`.

## Repository Layout

- `run_sequence_tm_pipeline.py`: main CLI
- `run_smoke_test.sh`: small end-to-end smoke run
- `run_full_pipeline.sh`: convenience wrapper for full runs
- `examples/full_run_template/`: starter input files
- `TEM_pdb.ent`, `APH_pdb.ent`: example reference structures
- `run_final_tm_scores.py`: imported experiment-specific source runner kept for provenance
- `docs/SCIENTIFIC_CONTEXT.md`: scientific motivation and TM-align interpretation
- `docs/USAGE.md`: setup, inputs, outputs, resumability, and troubleshooting

## Requirements

- Python 3.10+
- `TMalign` installed and available on `PATH`
- an `ESM_API_KEY` exported in the shell environment

Python dependencies are listed in `requirements.txt`.

Install them with:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Quick Start

### 1. Set the ESM API key

```bash
export ESM_API_KEY='YOUR_KEY_HERE'
```

### 2. Run the smoke test

```bash
./run_smoke_test.sh
```

This builds a small test set from the bundled TEM and APH reference structures and runs two unique fold-and-align tasks.

### 3. Run a real batch

Edit:

- `examples/full_run_template/targets.csv`
- `examples/full_run_template/sequences.csv`

Then run:

```bash
./run_full_pipeline.sh \
  examples/full_run_template/sequences.csv \
  examples/full_run_template/targets.csv \
  outputs/full_run
```

Or call the CLI directly:

```bash
python3 run_sequence_tm_pipeline.py run \
  --input examples/full_run_template/sequences.csv \
  --targets examples/full_run_template/targets.csv \
  --outdir outputs/full_run
```

## Input Format

Sequence table:

```csv
sequence_id,target_name,sequence
seq_0001,tem,HPETLVKV...
seq_0002,aph,MIEQDGLH...
```

Target table:

```csv
target_name,template_path,template_chain,description
tem,/abs/path/to/TEM_pdb.ent,A,TEM reference structure
aph,/abs/path/to/APH_pdb.ent,A,APH reference structure
```

Required sequence columns:

- `sequence_id`
- `target_name`
- `sequence`

Required target columns:

- `target_name`
- `template_path`

Optional target columns:

- `template_chain`
- `description`

If `template_chain` is provided, the pipeline writes a chain-filtered reference PDB into `prepared_templates/` and uses that structure for TM-align.

## Output Files

Each run writes:

- `normalized_input.csv`: validated and normalized input rows
- `sequence_results.csv`: one row per original input row
- `checkpoint.jsonl`: append-only per-unique-task checkpoint log
- `summary.json`: run-level and per-target summary statistics

Additional supporting outputs:

- `unique_results.csv`: deduplicated per-unique-task results
- `attempt_log.jsonl`: ESM API retry and exception log
- `folded_pdbs/`: cached predicted structures
- `prepared_templates/`: chain-filtered WT/reference structures when requested

## Resumability

The pipeline is designed for interrupted or long-running jobs.

If a run stops partway through:

- rerun the same command against the same `outdir`
- completed tasks will be skipped
- cached predicted PDBs will be reused
- aggregate outputs will be refreshed at the end of the rerun

By default, checkpoint entries with status `scored`, `fold_failed`, or `tmalign_failed` are treated as complete. If you want to retry failures, rerun with `--retry-failed`.

## Common Options

Useful flags on `run`:

- `--limit N`: process only the first `N` unique tasks after normalization
- `--prepare-only`: validate inputs and write placeholder outputs without calling the API or TM-align
- `--max-concurrent-requests N`: control ESM API concurrency
- `--max-concurrent-tmalign N`: control TM-align subprocess concurrency
- `--allow-duplicate-sequence-ids`: permit repeated `sequence_id` values
- `--retry-failed`: re-attempt checkpointed failures

Run `python3 run_sequence_tm_pipeline.py run --help` for the full CLI.

## Documentation

- Usage and operational details: `docs/USAGE.md`
- Scientific motivation and TM-align interpretation: `docs/SCIENTIFIC_CONTEXT.md`
- Release and publishing checklist: `docs/RELEASE_CHECKLIST.md`

## References

- Zhang Y, Skolnick J. TM-align: a protein structure alignment algorithm based on the TM-score. *Nucleic Acids Research*. 2005;33(7):2302-2309. https://pmc.ncbi.nlm.nih.gov/articles/PMC1084323/
- Lin Z, Akin H, Rao R, Hie B, Zhu Z, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science*. 2023;379(6637):1123-1130. https://www.science.org/doi/10.1126/science.ade2574
- Hu Y, Cheng W, Wang J, Liu Y. EasyNano: rapid epitope-targeted nanobody CDR design via differentiable distogram optimization with ESMFold2. *arXiv*. 2026. https://arxiv.org/abs/2606.12772
- Biohub fold endpoint used by this pipeline: `POST https://biohub.ai/api/v1/fold`

At the time of writing, I have not identified a standalone published methods paper for ESMFold2 itself. The documentation therefore cites the original ESMFold paper for the model family and a recent preprint that explicitly uses the `ESMFold2-Fast` / `ESMFold2` nomenclature.
