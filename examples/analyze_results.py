#!/usr/bin/env python3
"""Analyze genuine pipeline output; pending/failed tasks never become zero scores.

No embedded scores, measurements, plots, or example numerical results.
"""
import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, help="sequence_results.csv")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, help="Optional CSV: sequence_id,design_strategy (one row per ID)")
    parser.add_argument("--measurements", type=Path, help="Optional external CSV: sequence_id,activity (one row per ID; aggregate assay replicates explicitly first)")
    args = parser.parse_args()
    frame = pd.read_csv(args.results, dtype={"sequence_id": str, "fold_task_id": str}, keep_default_na=False)
    required = {"sequence_id", "target_name", "fold_task_id", "status", "tm_score", "aligned_length", "wt_length", "scoring_signature", "model"}
    if not required <= set(frame):
        parser.error(f"Missing pipeline columns: {sorted(required - set(frame))}")
    for column in ("tm_score", "aligned_length", "wt_length"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    args.outdir.mkdir(parents=True, exist_ok=True)
    for target, group in frame.groupby("target_name"):
        if group["scoring_signature"].nunique() > 1 or group["model"].nunique() > 1:
            parser.error(f"{target}: mixed model or reference/scoring configurations; analyze separately")
    # Report failures/pending at both row and unique target+sequence+model task levels.
    frame.groupby(["target_name", "status"], dropna=False).size().rename("input_rows").to_csv(args.outdir / "row_status_counts.csv")
    unique = frame.drop_duplicates("fold_task_id").copy()
    unique.groupby(["target_name", "status"], dropna=False).size().rename("unique_tasks").to_csv(args.outdir / "unique_status_counts.csv")
    failed = frame.loc[frame["status"].ne("scored") | frame["tm_score"].isna()]
    failed.to_csv(args.outdir / "unscored_rows.csv", index=False)
    scored = frame.loc[frame["status"].eq("scored") & frame["tm_score"].notna()].copy()
    scored["reference_coverage"] = scored["aligned_length"] / scored["wt_length"]
    if args.metadata:
        meta = pd.read_csv(args.metadata, dtype={"sequence_id": str}, keep_default_na=False)
        if not {"sequence_id", "design_strategy"} <= set(meta):
            parser.error("Metadata needs sequence_id,design_strategy")
        scored = scored.merge(meta[["sequence_id", "design_strategy"]], on="sequence_id", how="left", validate="many_to_one")
    else:
        scored["design_strategy"] = "all_candidates"
    # One task within each strategy; the same sequence may deliberately occur in two sets.
    candidates = scored.drop_duplicates(["target_name", "design_strategy", "fold_task_id"]).copy()
    candidates.to_csv(args.outdir / "scored_candidates.csv", index=False)
    candidates.groupby(["target_name", "design_strategy"], dropna=False)["tm_score"].agg(["count", "mean", "median", "min", "max"]).to_csv(args.outdir / "strategy_summary.csv")
    candidates.sort_values(["target_name", "tm_score"]).groupby("target_name", group_keys=False).head(5).to_csv(args.outdir / "lowest_scoring_candidates.csv", index=False)
    if args.measurements:
        if scored["sequence_id"].duplicated().any():
            parser.error("Activity joins require unique scored sequence_id values; resolve repeated row IDs explicitly")
        measurements = pd.read_csv(args.measurements, dtype={"sequence_id": str}, keep_default_na=False)
        if not {"sequence_id", "activity"} <= set(measurements):
            parser.error("Measurements need sequence_id,activity")
        joined = scored.merge(measurements[["sequence_id", "activity"]], on="sequence_id", how="left", validate="one_to_one", indicator=True)
        joined.to_csv(args.outdir / "structure_activity.csv", index=False)
        # Independent biological candidate IDs are retained here, even for equal sequences.
        # No correlation or claim of function is produced from the tiny teaching set.
    print(f"Scored input rows: {len(scored)}; scored unique tasks: {scored['fold_task_id'].nunique()}; unscored rows: {len(failed)}")


if __name__ == "__main__":
    main()
