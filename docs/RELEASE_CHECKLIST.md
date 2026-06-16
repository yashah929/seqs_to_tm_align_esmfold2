# Release Checklist

Use this checklist before publishing the repository or sharing it with collaborators.

## Repository Content

- Confirm that `README.md` matches the current CLI behavior.
- Confirm that `docs/USAGE.md` matches the current input and output files.
- Confirm that `docs/SCIENTIFIC_CONTEXT.md` still reflects the score definition used in code.
- Confirm that example files under `examples/full_run_template/` are still valid.
- Confirm that no local output directories or cached run artifacts are present.

## Environment And Dependencies

- Confirm that `requirements.txt` installs cleanly in a fresh virtual environment.
- Confirm that `python3 -m py_compile run_sequence_tm_pipeline.py` passes.
- Confirm that `TMalign` is named correctly in the documentation for the target environment.

## Functional Checks

- Run `python3 run_sequence_tm_pipeline.py run --help`
- Run `python3 run_sequence_tm_pipeline.py make-smoke-test --help`
- Run a `--prepare-only` validation against the example inputs.
- If possible, run a small live smoke test with a valid `ESM_API_KEY`.

## Security And Sharing

- Confirm that no API keys are hardcoded in source or docs.
- Confirm that no local absolute paths remain in files intended for general users, except where they are clearly example placeholders or bundled example assets.
- Confirm that `.gitignore` excludes local environments and generated outputs.

## GitHub Setup

- Add a repository description.
- Add repository topics if useful, for example:
  - `protein-structure`
  - `bioinformatics`
  - `esmfold`
  - `tmalign`
  - `structural-biology`
- Choose whether issues and pull requests should be enabled.
- Add a citation file if this repository becomes part of a published workflow.

## After Publishing

- Clone the repository into a fresh directory.
- Follow the README from scratch.
- Verify that a new user could get to a `--prepare-only` run without private context.
