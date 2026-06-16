#!/usr/bin/env python3

from __future__ import annotations

import argparse
import asyncio
import hashlib
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import aiohttp
import numpy as np
import pandas as pd
from tqdm.asyncio import tqdm_asyncio


ROOT = Path(__file__).resolve().parent
OUTPUTS_ROOT = ROOT / "outputs"
CONFIGS_ROOT = ROOT / "configs"
REFERENCE_ROOT = ROOT.parent.parent / "5_12_2026_dms_then_refinement_second_implementation"
PDB_DIR = REFERENCE_ROOT / "pdbs"
MSA_DIR = REFERENCE_ROOT / "MSAs"

LEGACY_ESM_API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
OFFICIAL_FOLD_API_URL = "https://biohub.ai/api/v1/fold"
TMALIGN_RE = re.compile(r"TM-score\s*=\s*([\d.]+)")
RETRY_HTTP_STATUSES = {429, 500, 502, 503, 504}
COOLDOWN_ON_403 = 60
AA3 = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS",
    "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN",
    "R": "ARG", "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}
TARGET_METADATA = {
    "tem": {
        "label": "TEM",
        "length": 263,
        "msa_path": MSA_DIR / "TEM.a3m",
        "template_path": PDB_DIR / "TEM_pdb.ent",
    },
    "aph": {
        "label": "APH",
        "length": 264,
        "msa_path": MSA_DIR / "APH.a3m",
        "template_path": PDB_DIR / "APH_pdb.ent",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute TM-scores for final TEM/APH production survivors.")
    parser.add_argument("--outdir", type=Path, default=ROOT / "tm_score_results_final1000", help="Output directory.")
    parser.add_argument("--folded-dir", type=Path, default=None, help="Optional fold cache directory.")
    parser.add_argument("--max-concurrent-requests", type=int, default=4, help="Max concurrent ESM API requests.")
    parser.add_argument("--max-concurrent-tmalign", type=int, default=4, help="Max concurrent TM-align processes.")
    parser.add_argument("--max-retries", type=int, default=8, help="Max ESM API retries per unique sequence.")
    parser.add_argument("--request-timeout", type=int, default=120, help="Per-request timeout in seconds.")
    parser.add_argument("--cooldown-seconds", type=float, default=2.0, help="Sleep between completed sequences.")
    parser.add_argument("--retry-sleep-seconds", type=float, default=10.0, help="Base sleep after retryable failures.")
    parser.add_argument("--limit-unique", type=int, default=None, help="Optional limit on unique fold tasks.")
    parser.add_argument("--progress-log-every", type=int, default=1, help="Emit a progress heartbeat every N completed unique tasks.")
    parser.add_argument("--progress-save-every", type=int, default=10, help="Refresh CSV checkpoints every N completed unique tasks.")
    parser.add_argument("--task-order", choices=["dms_desc", "wt_identity_desc", "input"], default="dms_desc")
    parser.add_argument("--prepare-only", action="store_true", help="Only build input tables.")
    return parser.parse_args()


def parse_first_a3m_sequence(path: Path) -> str:
    header_seen = False
    seq_chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                if header_seen:
                    break
                header_seen = True
                continue
            if header_seen:
                seq_chunks.append(line)
    if not seq_chunks:
        raise ValueError(f"No query sequence found in {path}")
    sequence = "".join(seq_chunks)
    return "".join(char for char in sequence if not char.islower()).replace(".", "-")


def load_wt_sequences() -> dict[str, str]:
    wt_sequences: dict[str, str] = {}
    for target_key, meta in TARGET_METADATA.items():
        wt_sequences[target_key] = parse_first_a3m_sequence(meta["msa_path"])
    return wt_sequences


def compute_identity(seq_a: str, seq_b: str) -> float:
    if len(seq_a) != len(seq_b):
        raise ValueError(f"Length mismatch: {len(seq_a)} vs {len(seq_b)}")
    return float(sum(left == right for left, right in zip(seq_a, seq_b, strict=True)) / len(seq_a))


def infer_target_key(sequence: str) -> str:
    seq_len = len(sequence)
    for target_key, meta in TARGET_METADATA.items():
        if seq_len == meta["length"]:
            return target_key
    raise ValueError(f"Could not infer target from sequence length {seq_len}")


def build_fold_task_id(target_key: str, sequence: str) -> str:
    digest = hashlib.sha1(f"{target_key}\n{sequence}".encode("utf-8")).hexdigest()[:16]
    return f"{target_key}__sha1_{digest}"


def load_design_state_table() -> pd.DataFrame:
    wt_sequences = load_wt_sequences()
    rows: list[dict[str, object]] = []
    for config_path in sorted(CONFIGS_ROOT.glob("*.json")):
        config = json.loads(config_path.read_text(encoding="utf-8"))
        variant = config_path.stem
        designs_path = OUTPUTS_ROOT / variant / "designs.json"
        payload = json.loads(designs_path.read_text(encoding="utf-8"))
        for phase2_row in payload["phase2"]:
            protein_a_seq = str(phase2_row["protein_a_final"])
            protein_b_seq = str(phase2_row["protein_b_final"])
            protein_a_target = infer_target_key(protein_a_seq)
            protein_b_target = infer_target_key(protein_b_seq)
            row = {
                "variant": variant,
                "config_path": str(config_path),
                "designs_path": str(designs_path),
                "arrangement": config["arrangement"],
                "offset": int(config["offset"]),
                "survivor_rank": int(phase2_row["rank"]),
                "seed": int(phase2_row["seed"]),
                "protein_a_sequence": protein_a_seq,
                "protein_b_sequence": protein_b_seq,
                "protein_a_target": protein_a_target,
                "protein_b_target": protein_b_target,
                "protein_a_label": TARGET_METADATA[protein_a_target]["label"],
                "protein_b_label": TARGET_METADATA[protein_b_target]["label"],
                "dms_score_a_final": float(phase2_row["dms_score_a_final"]),
                "dms_score_b_final": float(phase2_row["dms_score_b_final"]),
                "dms_score_total_final": float(phase2_row["dms_score_a_final"]) + float(phase2_row["dms_score_b_final"]),
                "nuc_seq_final": phase2_row["nuc_seq_final"],
                "wt_identity_non_overlap_reported": phase2_row.get("wt_identity_non_overlap"),
                "wt_identity_overlap_a_reported": phase2_row.get("wt_identity_overlap_a"),
                "wt_identity_overlap_b_reported": phase2_row.get("wt_identity_overlap_b"),
                "unique_aa_overlap_a": phase2_row.get("unique_aa_overlap_a"),
                "unique_aa_overlap_b": phase2_row.get("unique_aa_overlap_b"),
                "sd_mode": phase2_row.get("sd_mode"),
                "sd_realized_dG": phase2_row.get("sd_realized_dG"),
                "sd_realized_motif": phase2_row.get("sd_realized_motif"),
                "sd_spacing": phase2_row.get("sd_spacing"),
                "sd_start_codon": phase2_row.get("sd_start_codon"),
            }
            for slot in ["a", "b"]:
                seq = row[f"protein_{slot}_sequence"]
                target_key = row[f"protein_{slot}_target"]
                wt_seq = wt_sequences[target_key]
                row[f"protein_{slot}_wt_identity_recomputed"] = compute_identity(seq, wt_seq)
            rows.append(row)
    df = pd.DataFrame(rows)
    if df.empty:
        raise FileNotFoundError(f"No phase2 survivor rows found under {OUTPUTS_ROOT}")
    return df.sort_values(["variant", "survivor_rank"]).reset_index(drop=True)


def build_sequence_long_table(design_state_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for record in design_state_df.to_dict("records"):
        for slot in ["a", "b"]:
            target_key = str(record[f"protein_{slot}_target"])
            sequence = str(record[f"protein_{slot}_sequence"])
            rows.append(
                {
                    **record,
                    "sequence_slot": f"protein_{slot}",
                    "target": target_key,
                    "target_label": TARGET_METADATA[target_key]["label"],
                    "sequence": sequence,
                    "sequence_length": len(sequence),
                    "dms_score_target_final": float(record[f"dms_score_{slot}_final"]),
                    "wt_identity": float(record[f"protein_{slot}_wt_identity_recomputed"]),
                    "fold_task_id": build_fold_task_id(target_key, sequence),
                }
            )
    return pd.DataFrame(rows).sort_values(["variant", "survivor_rank", "sequence_slot"]).reset_index(drop=True)


def run_tmalign(template_pdb: Path, predicted_pdb: Path) -> float | None:
    result = subprocess.run(
        ["TMalign", str(predicted_pdb), str(template_pdb)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        return None
    tm_scores = TMALIGN_RE.findall(result.stdout)
    return float(tm_scores[1]) if len(tm_scores) >= 2 else None


def build_ca_pdb_from_fold_json(sequence: str, payload: dict[str, object]) -> str:
    coordinates = payload.get("coordinates")
    if not isinstance(coordinates, list):
        raise ValueError("Official fold response missing coordinates")
    lines: list[str] = []
    atom_serial = 1
    for residue_index, (aa, residue_atoms) in enumerate(zip(sequence, coordinates, strict=True), start=1):
        if not isinstance(residue_atoms, list) or len(residue_atoms) < 2:
            continue
        ca_atom = residue_atoms[1]
        if not isinstance(ca_atom, list) or len(ca_atom) != 3 or any(coord is None for coord in ca_atom):
            continue
        x, y, z = [float(coord) for coord in ca_atom]
        resname = AA3.get(aa, "UNK")
        lines.append(
            f"ATOM  {atom_serial:5d}  CA  {resname:>3s} A{residue_index:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 50.00           C"
        )
        atom_serial += 1
    lines.append("END")
    if atom_serial == 1:
        raise ValueError("Official fold response did not contain usable CA coordinates")
    return "\n".join(lines) + "\n"


def persist_outputs(
    outdir: Path,
    existing_df: pd.DataFrame,
    new_results: list[dict[str, object]],
    sequence_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, int]]:
    frames = [frame for frame in [existing_df, pd.DataFrame(new_results)] if not frame.empty]
    unique_results_df = deduplicate_results(pd.concat(frames, ignore_index=True)) if frames else pd.DataFrame()
    if "tm_score" not in unique_results_df.columns:
        unique_results_df["tm_score"] = np.nan
        unique_results_df["folded_pdb"] = None
    unique_results_df.to_csv(outdir / "unique_tm_scores.csv", index=False)

    mapped_cols = ["fold_task_id", "folded_pdb", "tm_score", "fold_attempts", "fold_elapsed_s", "tmalign_elapsed_s", "fold_error"]
    results_df = sequence_df.merge(unique_results_df[mapped_cols], on="fold_task_id", how="left", validate="many_to_one")
    results_df.to_csv(outdir / "tm_scores_by_sequence.csv", index=False)

    design_tm_df = summarize_by_design(results_df)
    design_tm_df.to_csv(outdir / "tm_scores_by_design.csv", index=False)

    target_summary_df, design_summary_df = summarize_by_variant(results_df, design_tm_df)
    target_summary_df.to_csv(outdir / "tm_scores_summary_by_variant_target.csv", index=False)
    design_summary_df.to_csv(outdir / "tm_scores_summary_by_variant_design.csv", index=False)

    summary = {
        "n_unique_scored": int(unique_results_df["tm_score"].notna().sum()),
        "n_sequence_records_scored": int(results_df["tm_score"].notna().sum()),
        "n_designs_scored_both": int(design_tm_df["tm_score_min_both"].notna().sum()) if "tm_score_min_both" in design_tm_df.columns else 0,
    }
    (outdir / "tm_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return unique_results_df, results_df, design_tm_df, summary


async def fold_sequence(
    session: aiohttp.ClientSession,
    semaphore: asyncio.Semaphore,
    api_key: str,
    fold_task_id: str,
    sequence: str,
    folded_dir: Path,
    max_retries: int,
    request_timeout: int,
    retry_sleep_seconds: float,
    attempt_log_path: Path,
) -> tuple[Path | None, int, str | None]:
    out_pdb = folded_dir / f"{fold_task_id}.pdb"
    if out_pdb.exists() and out_pdb.stat().st_size > 0:
        return out_pdb, 0, None
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}
    last_error: str | None = None
    for attempt in range(1, max_retries + 1):
        async with semaphore:
            try:
                async with session.post(
                    OFFICIAL_FOLD_API_URL,
                    headers=headers,
                    json={"sequence": sequence, "model": "esmfold2-fast-2026-05"},
                    timeout=request_timeout,
                ) as official_response:
                    if official_response.status == 200:
                        official_payload = await official_response.json()
                        pdb_text = build_ca_pdb_from_fold_json(sequence, official_payload)
                        out_pdb.write_text(pdb_text, encoding="utf-8")
                        return out_pdb, attempt, None
                    official_error_text = (await official_response.text())[:240].replace("\n", " ")
                    last_error = f"official HTTP {official_response.status}: {official_error_text}"
                    with attempt_log_path.open("a", encoding="utf-8") as handle:
                        handle.write(
                            json.dumps(
                                {
                                    "fold_task_id": fold_task_id,
                                    "attempt": attempt,
                                    "endpoint": "official",
                                    "status": official_response.status,
                                    "error": last_error,
                                    "timestamp_unix": time.time(),
                                }
                            )
                            + "\n"
                        )
                    if official_response.status == 403:
                        await asyncio.sleep(COOLDOWN_ON_403)
                        continue
                    if official_response.status in RETRY_HTTP_STATUSES:
                        await asyncio.sleep(retry_sleep_seconds * attempt)
                        continue
                    return None, attempt, last_error
            except (aiohttp.ClientError, asyncio.TimeoutError) as exc:
                last_error = f"{type(exc).__name__}: {exc}"
                with attempt_log_path.open("a", encoding="utf-8") as handle:
                    handle.write(
                        json.dumps(
                            {
                                "fold_task_id": fold_task_id,
                                "attempt": attempt,
                                "status": "exception",
                                "error": last_error,
                                "timestamp_unix": time.time(),
                            }
                        )
                        + "\n"
                    )
                await asyncio.sleep(retry_sleep_seconds * attempt)
            except Exception as exc:
                return None, attempt, f"{type(exc).__name__}: {exc}"
    return None, max_retries, last_error


def load_existing_results(outdir: Path) -> pd.DataFrame:
    csv_path = outdir / "unique_tm_scores.csv"
    jsonl_path = outdir / "unique_tm_progress.jsonl"
    frames: list[pd.DataFrame] = []
    if csv_path.exists():
        frames.append(pd.read_csv(csv_path))
    if jsonl_path.exists():
        rows = [json.loads(line) for line in jsonl_path.read_text(encoding="utf-8").splitlines() if line.strip()]
        if rows:
            frames.append(pd.DataFrame(rows))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def deduplicate_results(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    clean = df.copy()
    clean["tm_score_present"] = clean["tm_score"].notna()
    clean["folded_pdb_present"] = clean["folded_pdb"].notna()
    if "completed_unix" not in clean.columns:
        clean["completed_unix"] = 0.0
    clean = clean.sort_values(
        ["fold_task_id", "tm_score_present", "folded_pdb_present", "completed_unix"],
        ascending=[True, False, False, False],
    )
    return clean.drop_duplicates(subset=["fold_task_id"], keep="first").drop(columns=["tm_score_present", "folded_pdb_present"])


def order_unique_tasks(unique_df: pd.DataFrame, task_order: str) -> pd.DataFrame:
    if task_order == "input":
        return unique_df
    if task_order == "dms_desc":
        sort_cols = ["dms_score_total_final", "wt_identity", "variant", "survivor_rank", "sequence_slot"]
        ascending = [False, False, True, True, True]
    elif task_order == "wt_identity_desc":
        sort_cols = ["wt_identity", "dms_score_total_final", "variant", "survivor_rank", "sequence_slot"]
        ascending = [False, False, True, True, True]
    else:
        raise ValueError(f"Unsupported task order: {task_order}")
    return unique_df.sort_values(sort_cols, ascending=ascending, kind="mergesort").reset_index(drop=True)


async def fold_and_score(
    unique_df: pd.DataFrame,
    folded_dir: Path,
    api_key: str,
    max_concurrent_requests: int,
    max_concurrent_tmalign: int,
    max_retries: int,
    request_timeout: int,
    progress_path: Path,
    attempt_log_path: Path,
    outdir: Path,
    existing_df: pd.DataFrame,
    sequence_df: pd.DataFrame,
    progress_log_every: int,
    progress_save_every: int,
    cooldown_seconds: float,
    retry_sleep_seconds: float,
) -> pd.DataFrame:
    records = unique_df.to_dict("records")
    results: list[dict[str, object]] = []
    write_lock = asyncio.Lock()
    api_semaphore = asyncio.Semaphore(max_concurrent_requests)
    tmalign_semaphore = asyncio.Semaphore(max_concurrent_tmalign)
    connector = aiohttp.TCPConnector(limit=max(64, max_concurrent_requests * 4), ttl_dns_cache=300)
    timeout = aiohttp.ClientTimeout(total=None, connect=30, sock_connect=30, sock_read=request_timeout)
    async with aiohttp.ClientSession(connector=connector, timeout=timeout) as session:
        async def handle_record(record: dict[str, object]) -> None:
            start_s = time.time()
            pdb_path, attempts, fold_error = await fold_sequence(
                session=session,
                semaphore=api_semaphore,
                api_key=api_key,
                fold_task_id=str(record["fold_task_id"]),
                sequence=str(record["sequence"]),
                folded_dir=folded_dir,
                max_retries=max_retries,
                request_timeout=request_timeout,
                retry_sleep_seconds=retry_sleep_seconds,
                attempt_log_path=attempt_log_path,
            )
            fold_elapsed_s = time.time() - start_s
            tm_score = None
            tmalign_elapsed_s = 0.0
            if pdb_path is not None:
                align_start_s = time.time()
                async with tmalign_semaphore:
                    tm_score = await asyncio.to_thread(
                        run_tmalign,
                        TARGET_METADATA[str(record["target"])]["template_path"],
                        pdb_path,
                    )
                tmalign_elapsed_s = time.time() - align_start_s
            result = {
                **record,
                "folded_pdb": str(pdb_path) if pdb_path is not None else None,
                "tm_score": tm_score,
                "fold_attempts": attempts,
                "fold_elapsed_s": fold_elapsed_s,
                "tmalign_elapsed_s": tmalign_elapsed_s,
                "fold_error": fold_error,
                "completed_unix": time.time(),
            }
            results.append(result)
            async with write_lock:
                with progress_path.open("a", encoding="utf-8") as handle:
                    handle.write(json.dumps(result) + "\n")
                if progress_save_every > 0 and len(results) % progress_save_every == 0:
                    persist_outputs(outdir, existing_df, results, sequence_df)
                if progress_log_every > 0 and len(results) % progress_log_every == 0:
                    print(
                        f"[PROGRESS] completed={len(results)}/{len(records)} "
                        f"fold_task_id={record['fold_task_id']} tm_score={tm_score} fold_error={fold_error}",
                        flush=True,
                    )
            if cooldown_seconds > 0:
                await asyncio.sleep(cooldown_seconds)

        await tqdm_asyncio.gather(*[handle_record(record) for record in records])
    return pd.DataFrame(results)


def summarize_by_design(results_df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        results_df.pivot_table(
            index=["variant", "arrangement", "offset", "survivor_rank", "seed", "dms_score_total_final"],
            columns="target",
            values="tm_score",
            aggfunc="first",
        )
        .reset_index()
        .rename(columns={"tem": "tm_score_tem", "aph": "tm_score_aph"})
    )
    if "tm_score_tem" not in summary.columns:
        summary["tm_score_tem"] = np.nan
    if "tm_score_aph" not in summary.columns:
        summary["tm_score_aph"] = np.nan
    summary["tm_score_mean_both"] = summary[["tm_score_tem", "tm_score_aph"]].mean(axis=1)
    summary["tm_score_min_both"] = summary[["tm_score_tem", "tm_score_aph"]].min(axis=1)
    return summary.sort_values(["variant", "survivor_rank"]).reset_index(drop=True)


def summarize_by_variant(results_df: pd.DataFrame, design_df: pd.DataFrame) -> pd.DataFrame:
    target_summary = (
        results_df.groupby(["variant", "target"], dropna=False)
        .agg(
            n_sequences=("tm_score", "size"),
            n_scored=("tm_score", lambda s: int(s.notna().sum())),
            tm_score_mean=("tm_score", "mean"),
            tm_score_median=("tm_score", "median"),
            tm_score_min=("tm_score", "min"),
            tm_score_max=("tm_score", "max"),
            wt_identity_mean=("wt_identity", "mean"),
            dms_score_target_mean=("dms_score_target_final", "mean"),
        )
        .reset_index()
    )
    design_summary = (
        design_df.groupby("variant", dropna=False)
        .agg(
            n_designs=("survivor_rank", "size"),
            n_scored_both=("tm_score_min_both", lambda s: int(s.notna().sum())),
            tm_score_mean_both_mean=("tm_score_mean_both", "mean"),
            tm_score_mean_both_median=("tm_score_mean_both", "median"),
            tm_score_min_both_mean=("tm_score_min_both", "mean"),
            tm_score_min_both_median=("tm_score_min_both", "median"),
        )
        .reset_index()
    )
    return target_summary, design_summary


def main() -> None:
    args = parse_args()
    outdir = args.outdir.resolve()
    folded_dir = args.folded_dir.resolve() if args.folded_dir is not None else (outdir / "folded_pdbs")
    outdir.mkdir(parents=True, exist_ok=True)
    folded_dir.mkdir(parents=True, exist_ok=True)

    for meta in TARGET_METADATA.values():
        if not meta["template_path"].exists():
            raise FileNotFoundError(f"Missing template PDB: {meta['template_path']}")
        if not meta["msa_path"].exists():
            raise FileNotFoundError(f"Missing MSA: {meta['msa_path']}")

    design_state_df = load_design_state_table()
    sequence_df = build_sequence_long_table(design_state_df)
    unique_df = sequence_df.drop_duplicates(subset=["fold_task_id"], keep="first").copy()
    unique_df = order_unique_tasks(unique_df, args.task_order)
    if args.limit_unique is not None:
        unique_ids = set(unique_df.head(args.limit_unique)["fold_task_id"])
        unique_df = unique_df.loc[unique_df["fold_task_id"].isin(unique_ids)].copy()
        sequence_df = sequence_df.loc[sequence_df["fold_task_id"].isin(unique_ids)].copy()
        design_state_df = design_state_df.loc[design_state_df["survivor_rank"].isin(sequence_df["survivor_rank"])].copy()

    design_state_df.to_csv(outdir / "design_states.csv", index=False)
    sequence_df.to_csv(outdir / "tm_input_sequences.csv", index=False)
    unique_df.to_csv(outdir / "unique_tm_input_sequences.csv", index=False)

    print(
        f"[INFO] designs={len(design_state_df)}; "
        f"sequence records={len(sequence_df)}; unique fold tasks={len(unique_df)}"
    )

    if args.prepare_only:
        placeholder_unique = unique_df.copy()
        placeholder_unique["folded_pdb"] = None
        placeholder_unique["tm_score"] = np.nan
        placeholder_unique["fold_attempts"] = np.nan
        placeholder_unique["fold_elapsed_s"] = np.nan
        placeholder_unique["tmalign_elapsed_s"] = np.nan
        placeholder_unique["fold_error"] = None
        placeholder_unique.to_csv(outdir / "unique_tm_scores.csv", index=False)
        results_df = sequence_df.merge(
            placeholder_unique[
                ["fold_task_id", "folded_pdb", "tm_score", "fold_attempts", "fold_elapsed_s", "tmalign_elapsed_s", "fold_error"]
            ],
            on="fold_task_id",
            how="left",
            validate="many_to_one",
        )
        results_df.to_csv(outdir / "tm_scores_by_sequence.csv", index=False)
        design_tm_df = summarize_by_design(results_df)
        design_tm_df.to_csv(outdir / "tm_scores_by_design.csv", index=False)
        target_summary_df, design_summary_df = summarize_by_variant(results_df, design_tm_df)
        target_summary_df.to_csv(outdir / "tm_scores_summary_by_variant_target.csv", index=False)
        design_summary_df.to_csv(outdir / "tm_scores_summary_by_variant_design.csv", index=False)
        print("[INFO] prepare-only complete; no ESM API or TM-align calls were made.")
        return

    api_key = os.environ.get("ESM_API_KEY")
    if api_key is None:
        raise RuntimeError("ESM_API_KEY environment variable not set")
    if shutil.which("TMalign") is None:
        raise RuntimeError("TMalign not found in PATH")

    existing_df = deduplicate_results(load_existing_results(outdir))
    done_ids = set()
    if not existing_df.empty and "tm_score" in existing_df.columns:
        done_ids = set(existing_df.loc[existing_df["tm_score"].notna(), "fold_task_id"].astype(str))
    pending_df = unique_df.loc[~unique_df["fold_task_id"].isin(done_ids)].copy()
    print(
        f"[INFO] already scored unique tasks={len(done_ids)}; "
        f"pending unique tasks={len(pending_df)}; "
        f"api concurrency={args.max_concurrent_requests}; "
        f"tmalign concurrency={args.max_concurrent_tmalign}"
    )

    new_df = pd.DataFrame()
    if not pending_df.empty:
        new_df = asyncio.run(
            fold_and_score(
                unique_df=pending_df,
                folded_dir=folded_dir,
                api_key=api_key,
                max_concurrent_requests=args.max_concurrent_requests,
                max_concurrent_tmalign=args.max_concurrent_tmalign,
                max_retries=args.max_retries,
                request_timeout=args.request_timeout,
                progress_path=outdir / "unique_tm_progress.jsonl",
                attempt_log_path=outdir / "attempt_log.jsonl",
                outdir=outdir,
                existing_df=existing_df,
                sequence_df=sequence_df,
                progress_log_every=args.progress_log_every,
                progress_save_every=args.progress_save_every,
                cooldown_seconds=args.cooldown_seconds,
                retry_sleep_seconds=args.retry_sleep_seconds,
            )
        )

    unique_results_df, results_df, design_tm_df, progress_summary = persist_outputs(
        outdir=outdir,
        existing_df=existing_df,
        new_results=new_df.to_dict("records"),
        sequence_df=sequence_df,
    )
    summary = {
        "n_designs": int(len(design_state_df)),
        "n_sequence_records": int(len(sequence_df)),
        "n_unique_fold_tasks": int(len(unique_df)),
        **progress_summary,
    }
    (outdir / "tm_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(
        f"[INFO] scored unique tasks={summary['n_unique_scored']}/{summary['n_unique_fold_tasks']}; "
        f"designs with both TM-scores={summary['n_designs_scored_both']}/{summary['n_designs']}"
    )


if __name__ == "__main__":
    main()
