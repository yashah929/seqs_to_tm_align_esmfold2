#!/usr/bin/env python3
"""Render a real scored candidate using the TM-align transformation, never a refit.

Optional dependencies: docs/requirements-figures.txt. No network access.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import run_sequence_tm_pipeline as pipeline


def read_matrix(text):
    rows = []
    for line in text.splitlines():
        if re.match(r"^\s*[012]\s+", line):
            values = line.split()
            if len(values) == 5:
                rows.append([float(v) for v in values[1:]])
    if len(rows) != 3:
        raise ValueError("Could not parse the three TM-align rotation/translation rows")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--sequence-id", required=True)
    parser.add_argument("--protein", required=True)
    parser.add_argument("--reference-accession", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    frame = pd.read_csv(args.results, dtype={"sequence_id": str}, keep_default_na=False)
    selected = frame.loc[frame["sequence_id"].eq(args.sequence_id)]
    if len(selected) != 1 or selected.iloc[0]["status"] != "scored":
        parser.error("Select exactly one scored input row; pending/failed rows cannot be rendered")
    row = selected.iloc[0]
    prediction = Path(row["predicted_pdb_path"])
    reference = Path(row["template_path"])
    source = Path(row["reference_source_path"])
    sidecar = json.loads(prediction.with_suffix(".json").read_text())
    if any(sidecar.get(k) != v for k, v in pipeline.prediction_provenance(row["sequence"], row["model"]).items()):
        parser.error("Prediction metadata does not match the selected sequence/model")
    if row["prediction_sha256"] != pipeline.sha256_file(prediction) or sidecar.get("pdb_sha256") != row["prediction_sha256"]:
        parser.error("Prediction changed since scoring")
    chain = str(row["template_chain"]).strip() or None
    prepared_target = pipeline.TargetDefinition(row["target_name"], reference, chain)
    source_target = pipeline.TargetDefinition(row["target_name"], source, chain)
    executable = pipeline.shutil.which("TMalign")
    if not executable:
        parser.error("TMalign is required")
    signature = pipeline.scoring_signature(prepared_target, source_target, pipeline.sha256_file(Path(executable)))
    if signature != row["scoring_signature"]:
        parser.error("Reference or TM-align changed; rerun scoring before rendering")
    with tempfile.TemporaryDirectory() as tmp:
        matrix = Path(tmp) / "rotation.txt"
        result = subprocess.run([executable, str(prediction), str(reference), "-m", str(matrix)], check=True, capture_output=True, text=True)
        metrics = pipeline.parse_tmalign_output(result.stdout)
        if abs(metrics["tm_score"] - float(row["tm_score"])) > 1e-5:
            parser.error("Recomputed TM-score differs from saved result")
        transform = np.asarray(read_matrix(matrix.read_text()))
    def coordinates(path):
        return np.asarray([[float(line[s:s + 8]) for s in (30, 38, 46)] for line in pipeline.selected_ca_records(path.read_text())])
    pred_xyz = coordinates(prediction) @ transform[:, 1:].T + transform[:, 0]
    ref_xyz = coordinates(reference)
    # Shared principal-axis projection after the actual TM-align superposition.
    combined = np.vstack([pred_xyz, ref_xyz])
    center = combined.mean(axis=0)
    _, _, axes = np.linalg.svd(combined - center, full_matrices=False)
    pred_xy = (pred_xyz - center) @ axes[:2].T
    ref_xy = (ref_xyz - center) @ axes[:2].T
    fig, ax = plt.subplots(figsize=(10, 6), facecolor="white")
    for xyz, xy, color, label in ((ref_xyz, ref_xy, "#b15b16", "Experimental reference"), (pred_xyz, pred_xy, "#2864b4", "Candidate prediction")):
        # Avoid drawing false bonds across missing residues / chain breaks.
        xy = xy.copy()
        breaks = np.where(np.linalg.norm(np.diff(xyz, axis=0), axis=1) > 4.5)[0] + 1
        chunks = np.split(xy, breaks)
        for i, chunk in enumerate(chunks):
            ax.plot(chunk[:, 0], chunk[:, 1], color=color, lw=1.5, alpha=.85, label=label if i == 0 else None)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.legend(loc="upper right", frameon=False)
    fig.suptitle(f"{args.protein} · {args.sequence_id}\nReference {args.reference_accession} chain {chain or '(blank)'} · reference-normalized TM = {metrics['tm_score']:.5f}", fontsize=13)
    fig.text(.5, .035, f"{row['model']} · Cα-only prediction and Cα traces\nActual TM-align superposition; gaps > 4.5 Å are not connected", ha="center", fontsize=10)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=180, bbox_inches="tight", facecolor="white")
    args.out.with_suffix(".json").write_text(json.dumps({"sequence_id": args.sequence_id, "model": row["model"], "reference_accession": args.reference_accession, "chain": chain, "metrics": metrics, "scoring_signature": signature, "representation": "CA-only", "transform": transform.tolist()}, indent=2) + "\n")


if __name__ == "__main__":
    main()
