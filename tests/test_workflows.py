"""End-to-end OFFLINE FIXTURES: no API predictions or experimental assay results."""
import asyncio
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

import run_sequence_tm_pipeline as p


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def args(self, *extra):
        with patch("sys.argv", ["pipeline", "run", "--input", str(p.ROOT / "examples/full_run_template/sequences.csv"), "--targets", str(p.ROOT / "examples/full_run_template/targets.csv"), "--outdir", str(self.root), *extra]):
            return p.parse_args()

    def test_malformed_checkpoint_and_prepare_guard(self):
        checkpoint = self.root / "checkpoint.jsonl"
        checkpoint.write_text('{"fold_task_id": "fixture"}\n{"truncated":')
        with self.assertRaisesRegex(ValueError, "Invalid checkpoint JSON"):
            p.load_existing_checkpoint(checkpoint)
        with self.assertRaisesRegex(ValueError, "existing checkpoint"):
            p.run_pipeline(self.args("--prepare-only"))

    def test_prepare_only_no_api_or_tmalign_execution(self):
        with patch.object(p, "fold_sequence", side_effect=AssertionError("No folding allowed")), patch.object(p, "run_tmalign", side_effect=AssertionError("No alignment allowed")), patch.dict(os.environ, {}, clear=True):
            p.run_pipeline(self.args("--prepare-only"))
        frame = pd.read_csv(self.root / "sequence_results.csv")
        self.assertEqual(len(frame), 4)
        self.assertTrue(frame["status"].eq("pending").all())
        self.assertTrue(frame["tm_score"].isna().all())
        self.assertFalse((self.root / "checkpoint.jsonl").exists())

    def test_mocked_fold_failure_and_repeated_rows_in_analysis(self):
        """Synthetic metrics exercise file contracts only, never scientific results."""
        args = self.args("--prepare-only")
        p.run_pipeline(args)
        rows = pd.read_csv(self.root / "normalized_input.csv")
        normalized = pd.concat([rows.iloc[[0]], rows.iloc[[0]], rows.iloc[[1]]], ignore_index=True)
        normalized["sequence_id"] = ["NA", "fixture_ok_2", "fixture_failed"]
        normalized["input_row_index"] = [1, 2, 3]
        pending = normalized.drop_duplicates("fold_task_id")
        source = p.TargetDefinition("tem", p.ROOT / "TEM_pdb.ent", "A")
        target = p.prepare_template_pdb(source, self.root / "prepared_templates")
        async def fake_fold(**kwargs):
            if kwargs["sequence"].startswith("A"):
                return None, 1, "OFFLINE FIXTURE: simulated prediction failure"
            path = kwargs["out_pdb"]
            path.write_text(target.template_path.read_text())
            return path, 1, None
        metrics = {"tm_score": .75, "tm_score_query": .75, "aligned_length": 200, "rmsd": 2., "seq_identity_aligned": 1., "predicted_length": 263, "wt_length": 263}
        with patch.object(p, "fold_sequence", side_effect=fake_fold), patch.object(p, "run_tmalign", return_value=metrics):
            completed = asyncio.run(p.fold_and_score(
                pending, {"tem": target}, self.root / "folded_pdbs", "OFFLINE-FAKE-KEY", "OFFLINE-FIXTURE-NOT-A-MODEL",
                1, 1, 1, 1, 0., 0., 0, 1, self.root / "checkpoint.jsonl", self.root / "attempt_log.jsonl",
                self.root, normalized, pd.DataFrame()))
        _, results, summary = p.persist_outputs(self.root, normalized, pd.DataFrame(), completed.to_dict("records"))
        self.assertEqual(len(results), 3)
        self.assertEqual((summary["n_scored_rows"], summary["n_scored_unique_tasks"]), (2, 1))
        self.assertEqual(results.iloc[2]["status"], "fold_failed")
        self.assertTrue(pd.isna(results.iloc[2]["tm_score"]))
        analysis = self.root / "analysis"
        subprocess.run([sys.executable, str(p.ROOT / "examples/analyze_results.py"), str(self.root / "sequence_results.csv"), "--outdir", str(analysis)], check=True, capture_output=True)
        self.assertEqual(len(pd.read_csv(analysis / "scored_candidates.csv")), 1)
        self.assertEqual(len(pd.read_csv(analysis / "unscored_rows.csv")), 1)
        self.assertEqual(pd.read_csv(analysis / "strategy_summary.csv")["count"].iloc[0], 1)
        meta = self.root / "metadata.csv"
        meta.write_text("sequence_id,design_strategy\nNA,baseline\nfixture_ok_2,redesign\nfixture_failed,redesign\n")
        assays = self.root / "assays.csv"
        assays.write_text("sequence_id,activity\nNA,1.0\nfixture_ok_2,2.0\n")
        subprocess.run([sys.executable, str(p.ROOT / "examples/analyze_results.py"), str(self.root / "sequence_results.csv"), "--metadata", str(meta), "--measurements", str(assays), "--outdir", str(analysis)], check=True, capture_output=True)
        joined = pd.read_csv(analysis / "structure_activity.csv", keep_default_na=False)
        self.assertEqual(joined["sequence_id"].tolist(), ["NA", "fixture_ok_2"])
        self.assertTrue(joined["_merge"].eq("both").all())

    def test_reference_change_realigns_but_preserves_cache(self):
        """Mock folding delegates to real cache validation; network is forbidden."""
        targets_file = self.root / "targets.csv"
        reference_file = self.root / "reference.ent"
        reference_file.write_text((p.ROOT / "TEM_pdb.ent").read_text())
        targets_file.write_text("target_name,template_path,template_chain\ntem,reference.ent,A\n")
        args = self.args("--limit", "1")
        args.targets = targets_file
        args.input = self.root / "sequences.csv"
        pd.read_csv(p.ROOT / "examples/full_run_template/sequences.csv").iloc[:1].to_csv(args.input, index=False)
        prediction_dir = self.root / "folded_pdbs"
        prediction_dir.mkdir()
        sequence = p.parse_seqres_sequences(reference_file)["A"]
        task_id = p.stable_fold_task_id("tem", sequence, args.model)
        pred = prediction_dir / f"{task_id}.pdb"
        # Re-number the experimental reference solely as an OFFLINE cache fixture.
        records = p.selected_ca_records(reference_file.read_text(), "A")
        pred.write_text("REMARK OFFLINE FIXTURE NOT A MODEL PREDICTION\n" + "\n".join(line[:22] + f"{i:4d} " + line[27:] for i, line in enumerate(records, 1)) + "\nEND\n")
        pred.with_suffix(".json").write_text(json.dumps({**p.prediction_provenance(sequence, args.model), "pdb_sha256": p.sha256_file(pred)}))
        class ForbiddenNetwork:
            def post(self, *a, **kw):
                raise AssertionError("No network allowed")
        original_fold = p.fold_sequence
        async def cache_only(**kwargs):
            kwargs["session"] = ForbiddenNetwork()
            return await original_fold(**kwargs)
        metrics = {"tm_score": 1., "tm_score_query": 1., "aligned_length": 263, "rmsd": 0., "seq_identity_aligned": 1., "predicted_length": 263, "wt_length": 263}
        with patch.dict(os.environ, {"ESM_API_KEY": "OFFLINE-FAKE-KEY"}), patch.object(p.shutil, "which", return_value=sys.executable), patch.object(p, "fold_sequence", side_effect=cache_only) as fold, patch.object(p, "run_tmalign", return_value=metrics) as align:
            p.run_pipeline(args)
            self.assertEqual((fold.call_count, align.call_count), (1, 1))
            p.run_pipeline(args)
            self.assertEqual((fold.call_count, align.call_count), (1, 1))
            original_hash = p.sha256_file(pred)
            reference_file.write_text(reference_file.read_text() + "REMARK OFFLINE provenance change\n")
            with self.assertWarns(UserWarning):
                p.run_pipeline(args)
            self.assertEqual((fold.call_count, align.call_count), (2, 2))
            self.assertEqual(p.sha256_file(pred), original_hash)
            result = pd.read_csv(self.root / "sequence_results.csv").iloc[0]
            self.assertEqual(result["fold_attempts"], 0)

    @unittest.skipUnless(shutil.which("TMalign"), "optional local rotation-matrix integration")
    def test_renderer_reads_actual_tmalign_matrix(self):
        spec = importlib.util.spec_from_file_location("render_comparison", p.ROOT / "docs/assets/render_comparison.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        target = p.prepare_template_pdb(p.TargetDefinition("tem", p.ROOT / "TEM_pdb.ent", "A"), self.root)
        matrix = self.root / "matrix.txt"
        subprocess.run(["TMalign", str(target.template_path), str(target.template_path), "-m", str(matrix)], check=True, capture_output=True)
        values = module.read_matrix(matrix.read_text())
        self.assertEqual(len(values), 3)
        for i in range(3):
            self.assertAlmostEqual(values[i][0], 0., places=4)
            for j in range(3):
                self.assertAlmostEqual(values[i][j + 1], float(i == j), places=4)


if __name__ == "__main__":
    unittest.main()
