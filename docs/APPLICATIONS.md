# Reproducible application workflows

These are potential applications, not demonstrated design or experimental results. No live model predictions or assay measurements are bundled. Offline parser/coordinate fixtures in `tests/` are explicitly labeled and must not be used as scientific results.

## Small worked example

`examples/full_run_template/targets.csv` maps TEM to 1AXB A and APH to 1ND4 A using portable paths. `sequences.csv` contains their deposited chain-A sequences as prediction baselines and synthetic H1A/M1A teaching variants. See [construct provenance](SCIENTIFIC_CONTEXT.md#bundled-experimental-references); “H1A” is the first input residue, not TEM's author-numbered residue 1.

First validate locally:

```bash
python3 run_sequence_tm_pipeline.py run \
  --input examples/full_run_template/sequences.csv \
  --targets examples/full_run_template/targets.csv \
  --outdir outputs/prepared --prepare-only
python3 examples/analyze_results.py outputs/prepared/sequence_results.csv \
  --outdir outputs/prepared/analysis
```

Preparation yields four pending rows, four unique tasks, no scores, and references with 263 (TEM) and 255 (APH) observed Cα atoms. The analysis writes status counts and all pending rows to `unscored_rows.csv`; score summaries contain no scored candidates. These are preparation facts, not model results.

After setting a Biohub token and installing TM-align, run the four predictions explicitly:

```bash
./run_full_pipeline.sh examples/full_run_template/sequences.csv \
  examples/full_run_template/targets.csv outputs/example \
  --model esmfold2-fast-2026-05 --max-concurrent-requests 1
python3 examples/analyze_results.py outputs/example/sequence_results.csv \
  --outdir outputs/example/analysis
```

Inspect `sequence_results.csv`: each row contains `sequence_id,target_name,sequence,sequence_length,model,fold_task_id,status,tm_score,tm_score_query,aligned_length,rmsd,seq_identity_aligned,predicted_length,wt_length` plus provenance, paths, timing, and errors. A successful row has `status=scored`; failures retain empty scores and an `error`. Completion does not guarantee every candidate was scored.

For failed tasks inspect `error` and `attempt_log.jsonl`, resolve the cause, then rerun the live command with `--retry-failed`. No expected TM-score values are supplied: these must come from actual model output. Compare each teaching variant to its own deposited-sequence baseline, without treating a difference as an experimentally established mutation effect.

## Screen a redesigned library

Replace or extend the candidate table with your actual design sequences and unique IDs. Map each target deliberately to the appropriate experimental parent construct and chain. Include a deposited-sequence or independently verified WT-sequence prediction baseline; label which sequence was used. Use the same model for all candidates being compared.

```bash
./run_full_pipeline.sh my_library.csv my_targets.csv outputs/library
python3 examples/analyze_results.py outputs/library/sequence_results.csv \
  --outdir outputs/library/analysis
```

`strategy_summary.csv` reports count/mean/median/min/max for scored unique tasks within each target/set, not experimental replicates. `row_status_counts.csv` counts all input rows; `unique_status_counts.csv` counts target+sequence+model tasks. Repeated sequences with distinct IDs share folding/scoring, but still appear in `sequence_results.csv`; they are not independent predictions. Failed and pending rows are excluded from score summaries and exported separately, never imputed as zero. Compare failure rates as well as scored distributions to avoid selection bias.

## Compare candidate sets or design strategies

Keep candidate IDs unique and supply a separate CSV with `sequence_id,design_strategy`, one row per ID. For the bundled teaching set, a real metadata table (no scores) is included:

```bash
python3 examples/analyze_results.py outputs/example/sequence_results.csv \
  --metadata examples/full_run_template/metadata.csv \
  --outdir outputs/example/strategy_analysis
```

Use `baseline` versus `synthetic_teaching_variant` to exercise the join; two candidates per target cannot establish a useful distributional result. For larger real libraries, compare scored score distributions within a target using fixed reference contents, chain, model, and scoring settings. The script rejects mixed scoring signatures/models within a target. It deduplicates `fold_task_id` **within each strategy**, so an identical sequence deliberately included in two strategies contributes once to each, not as multiple independent predictions. Missing metadata forms a separate missing-label group.

The supplied analysis is tabular and executable without plotting dependencies. `scored_candidates.csv` preserves the per-set score distributions for ECDFs, histograms, or other analysis once genuine outputs exist. No plot or numerical scientific claim is supplied without scored inputs; record the unit of analysis and failure handling before adding statistical comparisons.

## Prioritize outliers for structural inspection

`lowest_scoring_candidates.csv` exports up to five scored candidates per target in ascending `tm_score`, with paths and `reference_coverage = aligned_length / wt_length`. These are relative outliers, not a universal pass/fail classification. Check input/reference coverage, missing predicted positions, domain truncations, construct selection, and model uncertainty before interpreting low scores as design defects.

Once a prediction is scored, render its actual alignment:

```bash
python3 -m pip install -r docs/requirements-figures.txt
python3 docs/assets/render_comparison.py \
  --results outputs/example/sequence_results.csv --sequence-id tem_wt \
  --protein TEM-1 --reference-accession 1AXB \
  --out outputs/example/tem_wt_comparison.svg
```

The script validates prediction and reference provenance, recomputes metrics with the same TM-align executable, and applies its `-m` rotation/translation to the candidate. It does not fit structures independently for display. Blue is prediction, orange is experimental reference; both are drawn as Cα traces. Labels include the user-supplied protein/accession, input candidate ID, selected chain, requested model, and computed reference-normalized TM-score. Trace gaps above 4.5 Å are not connected. Inspect the named PDB accession against the source header; the renderer checks coordinates/provenance, not the truth of a user-supplied display label. A companion JSON retains metrics and the transformation. Side-chain/ligand inspection needs fuller structures from appropriate sources.

## Combine structure with independent measurements

Provide your actual external assay CSV with `sequence_id,activity`, one row per candidate ID. Define units, controls, assay conditions, missing values, and replicate aggregation separately. No measurements are included in this repository.

```bash
python3 examples/analyze_results.py outputs/example/sequence_results.csv \
  --measurements path/to/independently_measured_activity.csv \
  --outdir outputs/example/assay_analysis
```

`structure_activity.csv` keeps scored candidate IDs, structure metrics, measured `activity`, and a merge indicator so unmatched assay data are visible. Repeated row IDs are rejected for this join; resolve them deliberately. Distinct biological candidate IDs can share a prediction if sequences match. Preserve that dependence when analyzing associations. Failed/pending candidates remain in `unscored_rows.csv` and do not acquire a fabricated structural score. The join does not establish causation or imply activity from TM-score.
