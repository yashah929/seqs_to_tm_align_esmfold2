# Cleanup behavior changes and validation

The general CLI entry point and primary output names are retained. The experiment-specific runner is kept at its original path with a provenance note; it has not been modernized.

## Material changes

- TM-align scores are parsed by explicit normalization labels, not positional score order; incomplete metric output is rejected. `tm_score` remains reference-normalized.
- Every reference is reduced to a single observed standard-ATOM Cα chain, with occupancy-based alternate selection. Multimodel files and ambiguous multichain references are rejected. HETATM residues are excluded, not treated as canonical protein positions.
- Coordinate conversion requires `[L,37,3]`, validates sequence length and finite Cα values, preserves missing-position numbering, and labels skipped all-null Cα positions. Coordinate-derived B factors changed from arbitrary 50 to placeholder 0 and are not confidence.
- Cached predictions require matching sequence/model/endpoint/conversion metadata and a PDB hash. Old bare PDB caches are not trusted automatically; a successful new fold replaces the unverified artifact. Scores track reference contents/chain/preparation policy and the TM-align executable. Reference-only changes retain valid predictions for re-alignment.
- Legacy checkpoint records without new provenance are warned about and excluded from verified results. They are retained in the log. Latest matching completion wins; corrupted checkpoint JSON produces a repair-oriented error.
- Required fields and runtime settings now reject missing values, unsafe target filenames, invalid sequences, and zero/negative concurrency, attempt counts, timeouts, or limits. Failed tasks keep empty metric columns.
- `--prepare-only` refuses to overwrite aggregates beside an existing checkpoint. The smoke wrapper now defaults to offline preparation; `--live` explicitly opts into two TEM predictions. Wrappers use Bash on macOS/Linux.
- Runnable teaching sequences replace placeholders; target paths are relative to their table. The duplicate `sequences.example.csv` and publishing checklist were removed; their tracked originals remain recoverable from git. The four-row `sequences.csv` is the single runnable example.

## Reproduce local checks (no credentials or network)

```bash
python3 -m unittest discover -s tests -v
python3 -m py_compile run_sequence_tm_pipeline.py run_final_tm_scores.py \
  examples/analyze_results.py docs/assets/generate_diagrams.py docs/assets/render_comparison.py
python3 run_sequence_tm_pipeline.py --help
python3 run_sequence_tm_pipeline.py run --help
python3 run_sequence_tm_pipeline.py make-smoke-test --help
bash -n run_full_pipeline.sh run_smoke_test.sh
./run_smoke_test.sh outputs/offline_smoke
python3 run_sequence_tm_pipeline.py run \
  --input examples/full_run_template/sequences.csv \
  --targets examples/full_run_template/targets.csv \
  --outdir outputs/prepared --prepare-only
python3 examples/analyze_results.py outputs/prepared/sequence_results.csv \
  --metadata examples/full_run_template/metadata.csv --outdir outputs/prepared/analysis
python3 docs/assets/generate_diagrams.py
git diff --check
```

The tests use standard-library unittest, core pipeline dependencies, small synthetic coordinate fixtures, and a labeled parser fixture obtained by aligning the two bundled experimental references with TM-align 20220412. That fixture is not a candidate prediction or evidence for an application result. Local TM-align integration is optional and skips if the executable is absent. API/cache behavior is mocked; tests never require credentials or network access.

The SVG generator uses only the Python standard library. Molecular rendering dependencies are optional and isolated in `requirements-figures.txt`. The conceptual comparison must not be relabeled as a result: replace it only after a genuine scored prediction has been rendered and its provenance checked. The exact optional one-prediction live smoke command is in the [README](../README.md#quick-start-offline-first); no live API call is part of this cleanup validation.

Validation on 2026-10-18: all 21 offline tests passed, including two local TM-align 20220412 integrations. Python syntax, Bash syntax, all CLI help commands, the offline wrapper, four-row preparation and downstream analysis, relative documentation/heading/image links, SVG parsing, and `git diff --check` passed. All three diagrams were visually inspected at README width. The renderer correctly rejected pending rows; its TM-align matrix reader was tested against a real local alignment. A genuine prediction rendering remains untested because no live predictions are available. No folding API requests, commits, pushes, or PRs were made.
