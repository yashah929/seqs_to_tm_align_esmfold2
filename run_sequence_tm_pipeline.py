#!/usr/bin/env python3

from __future__ import annotations

import argparse
import asyncio
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, median

import aiohttp
import pandas as pd
from tqdm.asyncio import tqdm_asyncio


ROOT = Path(__file__).resolve().parent
OFFICIAL_FOLD_API_URL = "https://biohub.ai/api/v1/fold"
DEFAULT_MODEL = "esmfold2-fast-2026-05"
RETRY_HTTP_STATUSES = {429, 500, 502, 503, 504}
COOLDOWN_ON_403 = 60.0
VALID_SEQUENCE_RE = re.compile(r"^[A-Z]+$")
TMALIGN_SCORE_RE = re.compile(r"TM-score\s*=\s*([0-9.]+)")
TMALIGN_LEN_RE = re.compile(r"Length of Chain_([12]):\s*(\d+)")
TMALIGN_ALIGNED_RE = re.compile(r"Aligned length=\s*(\d+),\s*RMSD=\s*([0-9.]+),\s*Seq_ID=n_identical/n_aligned=\s*([0-9.]+)")

AA1_TO_AA3 = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}
AA3_TO_AA1 = {value: key for key, value in AA1_TO_AA3.items()}


@dataclass(frozen=True)
class TargetDefinition:
    target_name: str
    template_path: Path
    template_chain: str | None = None
    description: str | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="General-purpose ESM API to TM-align pipeline for amino-acid sequences."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    smoke = subparsers.add_parser(
        "make-smoke-test",
        help="Write a self-contained smoke-test target table and sequence table from local template PDBs.",
    )
    smoke.add_argument("--outdir", type=Path, default=ROOT / "examples" / "smoke_test_inputs")
    smoke.add_argument("--tem-template", type=Path, default=ROOT / "TEM_pdb.ent")
    smoke.add_argument("--aph-template", type=Path, default=ROOT / "APH_pdb.ent")

    run = subparsers.add_parser("run", help="Fold sequences through the ESM API and score them with TM-align.")
    run.add_argument("--input", type=Path, required=True, help="CSV/TSV with sequence_id,target_name,sequence.")
    run.add_argument("--targets", type=Path, required=True, help="CSV/TSV with target_name,template_path[,template_chain,description].")
    run.add_argument("--outdir", type=Path, required=True, help="Output directory for normalized inputs, results, and checkpoints.")
    run.add_argument("--folded-dir", type=Path, default=None, help="Optional directory for cached predicted PDBs.")
    run.add_argument("--model", default=DEFAULT_MODEL, help="ESM folding model name to send to the API.")
    run.add_argument("--max-concurrent-requests", type=int, default=4, help="Max in-flight ESM API requests.")
    run.add_argument("--max-concurrent-tmalign", type=int, default=4, help="Max in-flight TM-align subprocesses.")
    run.add_argument("--max-retries", type=int, default=8, help="Max API retries per unique fold task.")
    run.add_argument("--request-timeout", type=int, default=120, help="Per-request timeout in seconds.")
    run.add_argument("--retry-sleep-seconds", type=float, default=10.0, help="Base retry sleep; multiplied by attempt.")
    run.add_argument("--cooldown-seconds", type=float, default=0.0, help="Optional sleep after each completed unique task.")
    run.add_argument("--progress-log-every", type=int, default=10, help="Print a progress heartbeat every N completed unique tasks.")
    run.add_argument("--progress-save-every", type=int, default=10, help="Refresh output CSVs every N completed unique tasks.")
    run.add_argument("--limit", type=int, default=None, help="Optional limit on unique tasks after normalization.")
    run.add_argument("--allow-duplicate-sequence-ids", action="store_true", help="Permit repeated sequence_id values in the input table.")
    run.add_argument("--retry-failed", action="store_true", help="Retry checkpointed tasks whose last recorded status was fold_failed or tmalign_failed.")
    run.add_argument("--prepare-only", action="store_true", help="Validate and write normalized inputs without calling the API or TM-align.")
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def normalize_sequence(raw: object) -> str:
    sequence = str(raw).upper().replace(" ", "").replace("\n", "")
    if not sequence:
        raise ValueError("Empty sequence")
    if not VALID_SEQUENCE_RE.fullmatch(sequence):
        raise ValueError(f"Sequence contains invalid characters: {sequence!r}")
    unknown = sorted(set(sequence) - set(AA1_TO_AA3))
    if unknown:
        raise ValueError(f"Sequence contains unsupported residues for CA-PDB conversion: {''.join(unknown)}")
    return sequence


def stable_fold_task_id(target_name: str, sequence: str, model: str) -> str:
    digest = hashlib.sha1(f"{target_name}\n{model}\n{sequence}".encode("utf-8")).hexdigest()[:16]
    return f"{target_name}__sha1_{digest}"


def prepare_template_pdb(template: TargetDefinition, prepared_dir: Path) -> TargetDefinition:
    if template.template_chain is None:
        return template
    prepared_dir.mkdir(parents=True, exist_ok=True)
    out_path = prepared_dir / f"{template.target_name}__chain_{template.template_chain}.pdb"
    kept_lines: list[str] = []
    with template.template_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            if raw_line.startswith(("ATOM", "HETATM", "TER")) and len(raw_line) > 21 and raw_line[21].strip() == template.template_chain:
                kept_lines.append(raw_line.rstrip("\n"))
    if not kept_lines:
        raise ValueError(
            f"No ATOM/HETATM records found for chain {template.template_chain!r} in {template.template_path}"
        )
    kept_lines.append("END")
    out_path.write_text("\n".join(kept_lines) + "\n", encoding="utf-8")
    return TargetDefinition(
        target_name=template.target_name,
        template_path=out_path,
        template_chain=template.template_chain,
        description=template.description,
    )


def resolve_existing_path(raw_path: object, base_dir: Path) -> Path:
    candidate = Path(str(raw_path))
    if not candidate.is_absolute():
        candidate = (base_dir / candidate).resolve()
    if not candidate.exists():
        raise FileNotFoundError(f"Missing file: {candidate}")
    return candidate


def load_targets(targets_path: Path) -> dict[str, TargetDefinition]:
    target_df = read_table(targets_path)
    required = {"target_name", "template_path"}
    missing = required - set(target_df.columns)
    if missing:
        raise ValueError(f"Targets table missing required columns: {sorted(missing)}")
    definitions: dict[str, TargetDefinition] = {}
    for row in target_df.to_dict("records"):
        target_name = str(row["target_name"]).strip()
        if not target_name:
            raise ValueError("Encountered blank target_name in targets table")
        if target_name in definitions:
            raise ValueError(f"Duplicate target_name in targets table: {target_name}")
        chain = row.get("template_chain")
        chain_value = str(chain).strip() if pd.notna(chain) and str(chain).strip() else None
        definitions[target_name] = TargetDefinition(
            target_name=target_name,
            template_path=resolve_existing_path(row["template_path"], targets_path.parent),
            template_chain=chain_value,
            description=str(row["description"]).strip() if "description" in row and pd.notna(row["description"]) else None,
        )
    return definitions


def build_normalized_input(input_path: Path, targets: dict[str, TargetDefinition], model: str) -> pd.DataFrame:
    input_df = read_table(input_path)
    required = {"sequence_id", "target_name", "sequence"}
    missing = required - set(input_df.columns)
    if missing:
        raise ValueError(f"Input table missing required columns: {sorted(missing)}")
    records: list[dict[str, object]] = []
    for input_row_index, row in enumerate(input_df.to_dict("records"), start=1):
        sequence_id = str(row["sequence_id"]).strip()
        target_name = str(row["target_name"]).strip()
        if not sequence_id:
            raise ValueError(f"Blank sequence_id on input row {input_row_index}")
        if target_name not in targets:
            raise ValueError(f"Unknown target_name {target_name!r} on input row {input_row_index}")
        sequence = normalize_sequence(row["sequence"])
        records.append(
            {
                "input_row_index": input_row_index,
                "sequence_id": sequence_id,
                "target_name": target_name,
                "template_path": str(targets[target_name].template_path),
                "template_chain": targets[target_name].template_chain,
                "sequence": sequence,
                "sequence_length": len(sequence),
                "model": model,
                "fold_task_id": stable_fold_task_id(target_name, sequence, model),
            }
        )
    normalized_df = pd.DataFrame(records)
    if normalized_df.empty:
        raise ValueError("Input table contains no sequence rows")
    return normalized_df


def validate_normalized_input(normalized_df: pd.DataFrame, allow_duplicate_sequence_ids: bool) -> None:
    duplicate_sequence_ids = normalized_df.loc[
        normalized_df["sequence_id"].duplicated(keep=False),
        "sequence_id",
    ].astype(str)
    if not allow_duplicate_sequence_ids and not duplicate_sequence_ids.empty:
        duplicates = sorted(duplicate_sequence_ids.unique())
        preview = ", ".join(duplicates[:10])
        raise ValueError(
            "Duplicate sequence_id values detected. "
            "Each input row should usually have a unique sequence_id. "
            f"Duplicates: {preview}"
            + (" ..." if len(duplicates) > 10 else "")
            + ". Use --allow-duplicate-sequence-ids to override."
        )


def load_existing_checkpoint(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def completed_statuses_for_run(retry_failed: bool) -> set[str]:
    if retry_failed:
        return {"scored"}
    return {"scored", "fold_failed", "tmalign_failed"}


def deduplicate_results(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    clean = df.copy()
    if "completed_unix" not in clean.columns:
        clean["completed_unix"] = 0.0
    clean["has_tm_score"] = clean["tm_score"].notna() if "tm_score" in clean.columns else False
    clean = clean.sort_values(["fold_task_id", "has_tm_score", "completed_unix"], ascending=[True, False, False])
    return clean.drop_duplicates(subset=["fold_task_id"], keep="first").drop(columns=["has_tm_score"])


def merge_results(normalized_df: pd.DataFrame, unique_results_df: pd.DataFrame) -> pd.DataFrame:
    result_columns = [
        "fold_task_id",
        "status",
        "tm_score",
        "tm_score_query",
        "aligned_length",
        "rmsd",
        "seq_identity_aligned",
        "predicted_length",
        "wt_length",
        "fold_attempts",
        "fold_elapsed_s",
        "tmalign_elapsed_s",
        "predicted_pdb_path",
        "error",
        "completed_unix",
    ]
    available = [column for column in result_columns if column in unique_results_df.columns]
    if not available:
        merged = normalized_df.copy()
        merged["status"] = "pending"
        return merged
    return normalized_df.merge(unique_results_df[available], on="fold_task_id", how="left", validate="many_to_one")


def summarize_results(results_df: pd.DataFrame) -> dict[str, object]:
    scored = results_df.loc[results_df["tm_score"].notna()].copy() if "tm_score" in results_df.columns else pd.DataFrame()
    summary: dict[str, object] = {
        "n_rows": int(len(results_df)),
        "n_unique_tasks": int(results_df["fold_task_id"].nunique()),
        "n_scored_rows": int(len(scored)),
        "n_scored_unique_tasks": int(scored["fold_task_id"].nunique()) if not scored.empty else 0,
        "status_counts": {},
        "targets": {},
    }
    if "status" in results_df.columns:
        summary["status_counts"] = {str(key): int(value) for key, value in results_df["status"].fillna("pending").value_counts().items()}
    for target_name, target_df in results_df.groupby("target_name", dropna=False):
        target_scored = target_df.loc[target_df["tm_score"].notna()].copy() if "tm_score" in target_df.columns else pd.DataFrame()
        target_summary: dict[str, object] = {
            "n_rows": int(len(target_df)),
            "n_unique_tasks": int(target_df["fold_task_id"].nunique()),
            "n_scored_rows": int(len(target_scored)),
        }
        if not target_scored.empty:
            scores = [float(value) for value in target_scored["tm_score"]]
            target_summary.update(
                {
                    "tm_score_mean": mean(scores),
                    "tm_score_median": median(scores),
                    "tm_score_min": min(scores),
                    "tm_score_max": max(scores),
                }
            )
        summary["targets"][str(target_name)] = target_summary
    return summary


def build_ca_pdb_from_coordinates(sequence: str, coordinates: list[object]) -> str:
    lines: list[str] = []
    atom_serial = 1
    for residue_index, (aa, residue_atoms) in enumerate(zip(sequence, coordinates, strict=True), start=1):
        if not isinstance(residue_atoms, list) or len(residue_atoms) < 2:
            continue
        ca_atom = residue_atoms[1]
        if not isinstance(ca_atom, list) or len(ca_atom) != 3 or any(coord is None for coord in ca_atom):
            continue
        x, y, z = [float(coord) for coord in ca_atom]
        lines.append(
            f"ATOM  {atom_serial:5d}  CA  {AA1_TO_AA3[aa]:>3s} A{residue_index:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 50.00           C"
        )
        atom_serial += 1
    if atom_serial == 1:
        raise ValueError("Fold response did not include usable CA coordinates")
    lines.append("END")
    return "\n".join(lines) + "\n"


def fold_payload_to_pdb(sequence: str, payload: dict[str, object]) -> str:
    pdb_text = payload.get("pdb")
    if isinstance(pdb_text, str) and pdb_text.strip():
        return pdb_text if pdb_text.endswith("\n") else pdb_text + "\n"
    coordinates = payload.get("coordinates")
    if isinstance(coordinates, list):
        return build_ca_pdb_from_coordinates(sequence, coordinates)
    raise ValueError("Fold response missing both pdb text and coordinates")


def parse_tmalign_output(stdout: str) -> dict[str, float | int | None]:
    tm_scores = [float(match) for match in TMALIGN_SCORE_RE.findall(stdout)]
    if len(tm_scores) < 2:
        raise ValueError("Could not parse TM-score values from TM-align output")
    lengths = {int(chain): int(length) for chain, length in TMALIGN_LEN_RE.findall(stdout)}
    aligned_match = TMALIGN_ALIGNED_RE.search(stdout)
    aligned_length = int(aligned_match.group(1)) if aligned_match else None
    rmsd = float(aligned_match.group(2)) if aligned_match else None
    seq_identity = float(aligned_match.group(3)) if aligned_match else None
    return {
        "tm_score_query": tm_scores[0],
        "tm_score": tm_scores[1],
        "predicted_length": lengths.get(1),
        "wt_length": lengths.get(2),
        "aligned_length": aligned_length,
        "rmsd": rmsd,
        "seq_identity_aligned": seq_identity,
    }


def run_tmalign(predicted_pdb: Path, template_pdb: Path) -> dict[str, float | int | None]:
    result = subprocess.run(
        ["TMalign", str(predicted_pdb), str(template_pdb)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        stderr = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(f"TMalign failed with exit code {result.returncode}: {stderr[:240]}")
    return parse_tmalign_output(result.stdout)


async def fold_sequence(
    session: aiohttp.ClientSession,
    semaphore: asyncio.Semaphore,
    api_key: str,
    sequence: str,
    model: str,
    out_pdb: Path,
    max_retries: int,
    request_timeout: int,
    retry_sleep_seconds: float,
    attempt_log_path: Path,
) -> tuple[Path | None, int, str | None]:
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
                    json={"sequence": sequence, "model": model},
                    timeout=request_timeout,
                ) as response:
                    if response.status == 200:
                        payload = await response.json()
                        pdb_text = fold_payload_to_pdb(sequence, payload)
                        out_pdb.write_text(pdb_text, encoding="utf-8")
                        return out_pdb, attempt, None
                    error_text = (await response.text())[:400].replace("\n", " ")
                    last_error = f"HTTP {response.status}: {error_text}"
                    with attempt_log_path.open("a", encoding="utf-8") as handle:
                        handle.write(
                            json.dumps(
                                {
                                    "fold_task_id": out_pdb.stem,
                                    "attempt": attempt,
                                    "status": response.status,
                                    "error": last_error,
                                    "timestamp_unix": time.time(),
                                }
                            )
                            + "\n"
                        )
                    if response.status == 403:
                        await asyncio.sleep(COOLDOWN_ON_403)
                        continue
                    if response.status in RETRY_HTTP_STATUSES:
                        await asyncio.sleep(retry_sleep_seconds * attempt)
                        continue
                    return None, attempt, last_error
            except (aiohttp.ClientError, asyncio.TimeoutError) as exc:
                last_error = f"{type(exc).__name__}: {exc}"
                with attempt_log_path.open("a", encoding="utf-8") as handle:
                    handle.write(
                        json.dumps(
                            {
                                "fold_task_id": out_pdb.stem,
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


def persist_outputs(
    outdir: Path,
    normalized_df: pd.DataFrame,
    existing_results_df: pd.DataFrame,
    new_results: list[dict[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    frames = [frame for frame in [existing_results_df, pd.DataFrame(new_results)] if not frame.empty]
    unique_results_df = deduplicate_results(pd.concat(frames, ignore_index=True)) if frames else pd.DataFrame()
    normalized_df.to_csv(outdir / "normalized_input.csv", index=False)
    unique_results_df.to_csv(outdir / "unique_results.csv", index=False)
    sequence_results_df = merge_results(normalized_df, unique_results_df)
    sequence_results_df.to_csv(outdir / "sequence_results.csv", index=False)
    summary = summarize_results(sequence_results_df)
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return unique_results_df, sequence_results_df, summary


async def fold_and_score(
    pending_df: pd.DataFrame,
    targets: dict[str, TargetDefinition],
    folded_dir: Path,
    api_key: str,
    model: str,
    max_concurrent_requests: int,
    max_concurrent_tmalign: int,
    max_retries: int,
    request_timeout: int,
    retry_sleep_seconds: float,
    cooldown_seconds: float,
    progress_log_every: int,
    progress_save_every: int,
    checkpoint_path: Path,
    attempt_log_path: Path,
    outdir: Path,
    normalized_df: pd.DataFrame,
    existing_results_df: pd.DataFrame,
) -> pd.DataFrame:
    records = pending_df.to_dict("records")
    results: list[dict[str, object]] = []
    write_lock = asyncio.Lock()
    api_semaphore = asyncio.Semaphore(max_concurrent_requests)
    tmalign_semaphore = asyncio.Semaphore(max_concurrent_tmalign)
    connector = aiohttp.TCPConnector(limit=max(64, max_concurrent_requests * 4), ttl_dns_cache=300)
    timeout = aiohttp.ClientTimeout(total=None, connect=30, sock_connect=30, sock_read=request_timeout)

    async with aiohttp.ClientSession(connector=connector, timeout=timeout) as session:
        async def handle_record(record: dict[str, object]) -> None:
            fold_task_id = str(record["fold_task_id"])
            predicted_pdb_path = folded_dir / f"{fold_task_id}.pdb"
            fold_start = time.time()
            pdb_path, attempts, fold_error = await fold_sequence(
                session=session,
                semaphore=api_semaphore,
                api_key=api_key,
                sequence=str(record["sequence"]),
                model=model,
                out_pdb=predicted_pdb_path,
                max_retries=max_retries,
                request_timeout=request_timeout,
                retry_sleep_seconds=retry_sleep_seconds,
                attempt_log_path=attempt_log_path,
            )
            fold_elapsed_s = time.time() - fold_start
            result: dict[str, object] = {
                **record,
                "status": "fold_failed" if pdb_path is None else "folded",
                "fold_attempts": attempts,
                "fold_elapsed_s": fold_elapsed_s,
                "tmalign_elapsed_s": 0.0,
                "predicted_pdb_path": str(pdb_path) if pdb_path is not None else None,
                "error": fold_error,
                "completed_unix": time.time(),
            }
            if pdb_path is not None:
                align_start = time.time()
                try:
                    async with tmalign_semaphore:
                        tmalign_result = await asyncio.to_thread(
                            run_tmalign,
                            pdb_path,
                            targets[str(record["target_name"])].template_path,
                        )
                    result.update(tmalign_result)
                    result["status"] = "scored"
                except Exception as exc:
                    result["status"] = "tmalign_failed"
                    result["error"] = f"{type(exc).__name__}: {exc}"
                result["tmalign_elapsed_s"] = time.time() - align_start
            results.append(result)
            async with write_lock:
                with checkpoint_path.open("a", encoding="utf-8") as handle:
                    handle.write(json.dumps(result) + "\n")
                completed = len(results)
                if progress_save_every > 0 and completed % progress_save_every == 0:
                    persist_outputs(outdir, normalized_df, existing_results_df, results)
                if progress_log_every > 0 and completed % progress_log_every == 0:
                    print(
                        f"[PROGRESS] completed={completed}/{len(records)} "
                        f"fold_task_id={fold_task_id} status={result['status']} "
                        f"tm_score={result.get('tm_score')}",
                        flush=True,
                    )
            if cooldown_seconds > 0:
                await asyncio.sleep(cooldown_seconds)

        await tqdm_asyncio.gather(*[handle_record(record) for record in records])
    return pd.DataFrame(results)


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_seqres_sequences(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            if not raw_line.startswith("SEQRES"):
                continue
            chain = raw_line[11].strip() or "A"
            residues = raw_line[19:].split()
            sequences.setdefault(chain, [])
            sequences[chain].extend(AA3_TO_AA1[residue] for residue in residues if residue in AA3_TO_AA1)
    return {chain: "".join(residues) for chain, residues in sequences.items() if residues}


def mutate_for_smoke(sequence: str) -> str:
    if not sequence:
        raise ValueError("Cannot mutate an empty sequence")
    replacement = "A" if sequence[0] != "A" else "G"
    return replacement + sequence[1:]


def make_smoke_test(outdir: Path, tem_template: Path, aph_template: Path) -> None:
    tem_path = tem_template.resolve()
    aph_path = aph_template.resolve()
    tem_sequences = parse_seqres_sequences(tem_path)
    aph_sequences = parse_seqres_sequences(aph_path)
    if "A" not in tem_sequences or "A" not in aph_sequences:
        raise ValueError("Expected chain A in both template PDBs for smoke test generation")
    targets_rows = [
        {
            "target_name": "tem",
            "template_path": str(tem_path),
            "template_chain": "A",
            "description": "TEM template from local PDB",
        },
        {
            "target_name": "aph",
            "template_path": str(aph_path),
            "template_chain": "A",
            "description": "APH template from local PDB",
        },
    ]
    sequence_rows = [
        {"sequence_id": "tem_wt", "target_name": "tem", "sequence": tem_sequences["A"]},
        {"sequence_id": "tem_smoke_mutant", "target_name": "tem", "sequence": mutate_for_smoke(tem_sequences["A"])},
        {"sequence_id": "aph_wt", "target_name": "aph", "sequence": aph_sequences["A"]},
        {"sequence_id": "aph_smoke_mutant", "target_name": "aph", "sequence": mutate_for_smoke(aph_sequences["A"])},
    ]
    write_csv(outdir / "targets.csv", targets_rows, ["target_name", "template_path", "template_chain", "description"])
    write_csv(outdir / "sequences.csv", sequence_rows, ["sequence_id", "target_name", "sequence"])
    print(f"[INFO] wrote smoke test targets to {outdir / 'targets.csv'}")
    print(f"[INFO] wrote smoke test sequences to {outdir / 'sequences.csv'}")


def run_pipeline(args: argparse.Namespace) -> None:
    outdir = args.outdir.resolve()
    folded_dir = args.folded_dir.resolve() if args.folded_dir else (outdir / "folded_pdbs")
    outdir.mkdir(parents=True, exist_ok=True)
    folded_dir.mkdir(parents=True, exist_ok=True)

    targets = {
        target_name: prepare_template_pdb(target, outdir / "prepared_templates")
        for target_name, target in load_targets(args.targets.resolve()).items()
    }
    normalized_df = build_normalized_input(args.input.resolve(), targets, args.model)
    validate_normalized_input(normalized_df, allow_duplicate_sequence_ids=args.allow_duplicate_sequence_ids)
    unique_df = normalized_df.drop_duplicates(subset=["fold_task_id"], keep="first").reset_index(drop=True)
    if args.limit is not None:
        unique_ids = set(unique_df.head(args.limit)["fold_task_id"])
        unique_df = unique_df.loc[unique_df["fold_task_id"].isin(unique_ids)].reset_index(drop=True)
        normalized_df = normalized_df.loc[normalized_df["fold_task_id"].isin(unique_ids)].reset_index(drop=True)
    normalized_df.to_csv(outdir / "normalized_input.csv", index=False)

    placeholder_results_df = pd.DataFrame(
        {
            "fold_task_id": unique_df["fold_task_id"],
            "status": "pending",
            "tm_score": None,
            "tm_score_query": None,
            "aligned_length": None,
            "rmsd": None,
            "seq_identity_aligned": None,
            "predicted_length": None,
            "wt_length": None,
            "fold_attempts": None,
            "fold_elapsed_s": None,
            "tmalign_elapsed_s": None,
            "predicted_pdb_path": None,
            "error": None,
            "completed_unix": None,
        }
    )
    if args.prepare_only:
        persist_outputs(outdir, normalized_df, placeholder_results_df, [])
        print("[INFO] prepare-only complete; normalized inputs and placeholder outputs were written.")
        return

    api_key = os.environ.get("ESM_API_KEY")
    if not api_key:
        raise RuntimeError("ESM_API_KEY environment variable not set")
    if shutil.which("TMalign") is None:
        raise RuntimeError("TMalign not found in PATH")

    checkpoint_path = outdir / "checkpoint.jsonl"
    attempt_log_path = outdir / "attempt_log.jsonl"
    existing_results_df = deduplicate_results(load_existing_checkpoint(checkpoint_path))
    completed_ids = set()
    if not existing_results_df.empty and "status" in existing_results_df.columns:
        completed_statuses = completed_statuses_for_run(retry_failed=args.retry_failed)
        completed_ids = set(
            existing_results_df.loc[
                existing_results_df["status"].isin(completed_statuses),
                "fold_task_id",
            ].astype(str)
        )
    pending_df = unique_df.loc[~unique_df["fold_task_id"].isin(completed_ids)].reset_index(drop=True)

    print(
        f"[INFO] rows={len(normalized_df)} unique_tasks={len(unique_df)} "
        f"completed={len(completed_ids)} pending={len(pending_df)} model={args.model}",
        flush=True,
    )

    new_results_df = pd.DataFrame()
    if not pending_df.empty:
        new_results_df = asyncio.run(
            fold_and_score(
                pending_df=pending_df,
                targets=targets,
                folded_dir=folded_dir,
                api_key=api_key,
                model=args.model,
                max_concurrent_requests=args.max_concurrent_requests,
                max_concurrent_tmalign=args.max_concurrent_tmalign,
                max_retries=args.max_retries,
                request_timeout=args.request_timeout,
                retry_sleep_seconds=args.retry_sleep_seconds,
                cooldown_seconds=args.cooldown_seconds,
                progress_log_every=args.progress_log_every,
                progress_save_every=args.progress_save_every,
                checkpoint_path=checkpoint_path,
                attempt_log_path=attempt_log_path,
                outdir=outdir,
                normalized_df=normalized_df,
                existing_results_df=existing_results_df,
            )
        )

    _, _, summary = persist_outputs(outdir, normalized_df, existing_results_df, new_results_df.to_dict("records"))
    print(
        f"[INFO] finished unique_scored={summary['n_scored_unique_tasks']}/{summary['n_unique_tasks']} "
        f"status_counts={summary['status_counts']}",
        flush=True,
    )


def main() -> None:
    args = parse_args()
    if args.command == "make-smoke-test":
        make_smoke_test(args.outdir.resolve(), args.tem_template, args.aph_template)
        return
    if args.command == "run":
        run_pipeline(args)
        return
    raise ValueError(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    main()
