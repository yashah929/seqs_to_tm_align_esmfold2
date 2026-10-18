#!/usr/bin/env python3
"""Generate editable README diagrams using only the Python standard library.

These diagrams contain no molecular coordinates, predictions, or numerical results.
"""
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parent
INK = "#203348"
BLUE = "#2864b4"  # candidate / predicted structure
ORANGE = "#b15b16"  # experimental reference
TEAL = "#18786b"


def start(title, description, height):
    return [f'''<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="{height}" viewBox="0 0 1000 {height}" role="img" aria-labelledby="title desc">
<title id="title">{escape(title)}</title><desc id="desc">{escape(description)}</desc>
<defs><marker id="arrow" markerWidth="9" markerHeight="9" refX="7" refY="3" orient="auto"><path d="M0,0 L0,6 L8,3 Z" fill="#607387"/></marker></defs>
<rect width="1000" height="{height}" fill="white"/>
<g font-family="Arial, Helvetica, sans-serif" fill="{INK}">''']


def text(parts, x, y, lines, size=18, color=INK, weight="normal"):
    for i, line in enumerate(lines):
        parts.append(f'<text x="{x}" y="{y + i * (size + 8)}" font-size="{size}" fill="{color}" font-weight="{weight}">{escape(line)}</text>')


def card(parts, x, y, w, h, title, lines, color=BLUE, fill="#f3f7fd"):
    parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="12" fill="{fill}" stroke="{color}" stroke-width="1.5"/>')
    text(parts, x + 18, y + 31, [title], 20, color, "bold")
    text(parts, x + 18, y + 64, lines)


def arrow(parts, path, dashed=False):
    dash = ' stroke-dasharray="6 5"' if dashed else ""
    parts.append(f'<path d="{path}" fill="none" stroke="#607387" stroke-width="2"{dash} marker-end="url(#arrow)"/>')


def save(name, parts):
    (ROOT / name).write_text("\n".join([*parts, "</g></svg>"]) + "\n")


def workflow():
    parts = start("From candidate sequence to reference-anchored similarity", "Blue is the candidate prediction branch; orange is the experimental reference branch. Both feed TM-align. Cache and checkpoint provenance support resumption.", 555)
    text(parts, 30, 38, ["Computational assessment of predicted fold preservation"], 25, weight="bold")
    card(parts, 30, 75, 295, 130, "1  Candidate sequence table", ["sequence_id · target_name", "amino-acid sequence"])
    card(parts, 380, 75, 260, 130, "2  ESMFold2 prediction", ["Configured API model", "Retries on transient errors"])
    card(parts, 695, 75, 275, 130, "3  Prediction cache", ["Cα PDB + metadata", "Sequence / model / file hash"])
    arrow(parts, "M325 140 H375")
    arrow(parts, "M640 140 H690")
    card(parts, 30, 260, 370, 145, "Experimental WT reference", ["Target → PDB file + selected chain", "Observed protein Cα coordinates", "Single model; explicit chain selection"], ORANGE, "#fff7ee")
    arrow(parts, "M175 205 V253")
    text(parts, 190, 237, ["target mapping"], 16, ORANGE)
    card(parts, 465, 260, 220, 145, "4  TM-align", ["Prediction = Chain_1", "Reference = Chain_2", "WT-normalized score"], TEAL, "#f0f8f5")
    arrow(parts, "M400 330 H458")
    arrow(parts, "M825 205 V230 H575 V253")
    card(parts, 750, 260, 220, 145, "5  Results", ["Rows + unique tasks", "TM · RMSD · lengths", "Separate failure counts"], TEAL, "#f0f8f5")
    arrow(parts, "M685 330 H743")
    parts.append('<rect x="30" y="452" width="940" height="72" rx="10" fill="#f4f6f8"/>')
    text(parts, 50, 480, ["Checkpoint + resume"], 18, weight="bold")
    text(parts, 285, 480, ["Reuse verified scores; re-align after reference or scorer changes.", "Keep valid cached predictions when only alignment inputs change."], 17)
    arrow(parts, "M860 405 V445", dashed=True)
    save("workflow.svg", parts)


def applications():
    parts = start("Potential applications", "Four potential downstream workflows use scored candidates. Experimental activity is an independent external input, not inferred by this pipeline.", 470)
    text(parts, 30, 38, ["Potential applications · compare within a fixed target and reference"], 24, weight="bold")
    card(parts, 30, 70, 455, 125, "Screen redesigned libraries", ["Rank candidates by predicted fold similarity.", "Include a predicted WT-sequence baseline."])
    card(parts, 515, 70, 455, 125, "Compare candidate sets", ["Compare distributions across design strategies.", "Use the same reference, chain, and model."])
    card(parts, 30, 225, 455, 125, "Inspect low-scoring outliers", ["Check alignment, coverage, and truncations.", "Inspect local active sites and interfaces."])
    card(parts, 515, 225, 455, 125, "Relate structure to measured properties", ["Join scored candidates to independent assays.", "Test associations; do not infer activity from TM."], TEAL, "#f0f8f5")
    parts.append('<rect x="30" y="390" width="940" height="54" rx="10" fill="#fff7ee" stroke="#b15b16"/>')
    text(parts, 50, 423, ["External input: independently measured activity, binding, stability, or other properties"], 19, ORANGE)
    arrow(parts, "M745 390 V357", dashed=True)
    save("applications.svg", parts)


def conceptual():
    parts = start("Conceptual structural comparison — no prediction available", "Diagram only. A real comparison requires an API prediction aligned to the experimental reference; no structure or score is depicted.", 300)
    text(parts, 30, 38, ["Conceptual comparison · no model-generated result shown"], 24, weight="bold")
    card(parts, 30, 80, 280, 125, "Candidate prediction", ["TEM: tem_wt baseline", "esmfold2-fast-2026-05"], BLUE)
    card(parts, 690, 80, 280, 125, "Experimental reference", ["TEM-1 · PDB 1AXB · A", "X-ray · inhibitor-bound"], ORANGE, "#fff7ee")
    card(parts, 370, 80, 260, 125, "Structural superposition", ["Requires a real prediction", "No TM-score reported"], TEAL, "#f0f8f5")
    arrow(parts, "M310 140 H363")
    arrow(parts, "M690 140 H637")
    text(parts, 30, 250, ["Once predictions exist, render_comparison.py uses the actual TM-align rotation and metrics.", "The pipeline retains Cα-only predictions; the rendering shows Cα traces."], 18)
    save("comparison_concept.svg", parts)


if __name__ == "__main__":
    workflow()
    applications()
    conceptual()
