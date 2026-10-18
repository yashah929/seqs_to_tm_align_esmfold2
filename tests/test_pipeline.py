"""Offline reliability tests. Synthetic coordinates are fixtures, never model results."""
import asyncio
import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

import run_sequence_tm_pipeline as p


def atom37(n=4):
    return [[[float(i * 3), float(j), 0.] for j in range(37)] for i in range(n)]


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_reference_normalization_labels_not_line_order(self):
        fixture = Path(__file__).parent / "fixtures/tmalign_reference_pair.txt"
        text = fixture.read_text()
        result = p.parse_tmalign_output(text)
        self.assertEqual(result["tm_score"], .30917)
        self.assertEqual(result["tm_score_query"], .303)
        self.assertEqual((result["predicted_length"], result["wt_length"]), (263, 255))
        self.assertEqual(result["aligned_length"], 148)
        self.assertEqual(result["rmsd"], 6.78)
        self.assertEqual(result["seq_identity_aligned"], .054)
        lines = text.splitlines()
        lines[-1], lines[-2] = lines[-2], lines[-1]
        self.assertEqual(p.parse_tmalign_output("\n".join(lines)), result)
        with self.assertRaises(ValueError):
            p.parse_tmalign_output("TM-score=0.9\nTM-score=0.8")

    def test_coordinate_conversion_ca_index_and_bfactor(self):
        text = p.build_ca_pdb_from_coordinates("ACDE", atom37())
        ca = p.selected_ca_records(text)
        self.assertEqual(float(ca[0][38:46]), 1.)
        self.assertEqual(float(ca[0][60:66]), 0.)
        p.validate_prediction_pdb("ACDE", text)

    def test_bad_coordinates(self):
        for coords in (atom37(3), [[[0., 0., 0.]] * 3] * 4):
            with self.subTest(coords=coords), self.assertRaises(ValueError):
                p.build_ca_pdb_from_coordinates("ACDE", coords)
        for value in (float("nan"), float("inf"), None, 10000.):
            coords = atom37()
            coords[0][1][0] = value
            with self.subTest(value=value), self.assertRaises((ValueError, TypeError)):
                p.build_ca_pdb_from_coordinates("ACDE", coords)

    def test_missing_ca_preserves_sequence_number(self):
        coords = atom37()
        coords[1][1] = [None] * 3
        with self.assertWarns(UserWarning):
            text = p.build_ca_pdb_from_coordinates("ACDE", coords)
        self.assertEqual([int(x[22:26]) for x in p.selected_ca_records(text)], [1, 3, 4])

    def test_reject_models_chains_sequence_mismatch(self):
        text = p.build_ca_pdb_from_coordinates("ACDE", atom37())
        with self.assertRaises(ValueError):
            p.selected_ca_records("MODEL        1\n" + text + "MODEL        2\n" + text)
        with self.assertRaises(ValueError):
            p.selected_ca_records(text + text.replace(" A   ", " B   "))
        with self.assertRaises(ValueError):
            p.validate_prediction_pdb("AAAA", text)
        with self.assertRaises(ValueError):
            p.selected_ca_records(text + text)

    def test_alternate_locations_highest_occupancy(self):
        text = p.build_ca_pdb_from_coordinates("ACDE", atom37())
        first, *rest = p.selected_ca_records(text)
        alt = first[:16] + "B" + first[17:54] + "  0.75" + first[60:]
        other = first[:16] + "A" + first[17:54] + "  0.25" + first[60:]
        records = p.selected_ca_records("\n".join([other, alt, *rest]))
        self.assertEqual(float(records[0][54:60]), .75)
        self.assertEqual(records[0][16], " ")

    def test_required_values_and_portable_examples(self):
        for value in (None, float("nan"), "", " ", "A_X", "NA", "ACX"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                p.normalize_sequence(value)
        targets = p.load_targets(p.ROOT / "examples/full_run_template/targets.csv")
        df = p.build_normalized_input(p.ROOT / "examples/full_run_template/sequences.csv", targets, p.DEFAULT_MODEL)
        self.assertEqual(len(df), 4)
        self.assertEqual(df["sequence_length"].tolist(), [263, 263, 264, 264])
        for field in ("sequence_id", "target_name", "sequence"):
            path = self.root / "bad.csv"
            row = {"sequence_id": "test", "target_name": "tem", "sequence": "ACDE"}
            row[field] = ""
            pd.DataFrame([row]).to_csv(path, index=False)
            with self.subTest(field=field), self.assertRaises(ValueError):
                p.build_normalized_input(path, targets, p.DEFAULT_MODEL)
        with self.assertRaises(ValueError):
            p.build_normalized_input(self.write_table("input.csv", "sequence_id,target_name,sequence\nx,unknown,ACDE\n"), targets, p.DEFAULT_MODEL)
        for name, path in (("", "ref.pdb"), ("../bad", "ref.pdb"), ("tem", "")):
            with self.assertRaises((ValueError, FileNotFoundError)):
                p.load_targets(self.write_table("targets.csv", f"target_name,template_path\n{name},{path}\n"))

    def write_table(self, name, text):
        path = self.root / name
        path.write_text(text)
        return path

    def test_runtime_options(self):
        with patch("sys.argv", ["pipeline", "run", "--input", "i", "--targets", "t", "--outdir", "o"]):
            args = p.parse_args()
        p.validate_runtime_options(args)
        for name, value in (("limit", 0), ("max_retries", 0), ("max_concurrent_requests", -1), ("cooldown_seconds", float("nan"))):
            with self.subTest(name=name), patch.object(args, name, value), self.assertRaises(ValueError):
                p.validate_runtime_options(args)

    def test_reference_policy_and_signature(self):
        source = p.TargetDefinition("aph", p.ROOT / "APH_pdb.ent", "A")
        prepared = p.prepare_template_pdb(source, self.root)
        self.assertEqual(len(p.selected_ca_records(prepared.template_path.read_text())), 255)
        sig = p.scoring_signature(prepared, source, "binary1")
        self.assertNotEqual(sig, p.scoring_signature(prepared, source, "binary2"))
        copied = self.root / "source.ent"
        copied.write_text(source.template_path.read_text() + "REMARK changed\n")
        self.assertNotEqual(sig, p.scoring_signature(prepared, p.TargetDefinition("aph", copied, "A"), "binary1"))
        self.assertNotEqual(sig, p.scoring_signature(prepared, p.TargetDefinition("aph", source.template_path, "B"), "binary1"))

    def test_checkpoint_requires_current_provenance_and_prediction_hash(self):
        pdb = self.write_table("pred.pdb", p.build_ca_pdb_from_coordinates("ACDE", atom37()))
        row = {"fold_task_id": "id", "scoring_signature": "current", "status": "scored", "predicted_pdb_path": str(pdb), "prediction_sha256": p.sha256_file(pdb)}
        current = pd.DataFrame([row])
        self.assertEqual(len(p.verified_checkpoint(current, current)), 1)
        for changed in ({"scoring_signature": "old"}, {"prediction_sha256": "bad"}, {"fold_task_id": "different-sequence-or-model"}):
            with self.subTest(changed=changed), self.assertWarns(UserWarning):
                self.assertTrue(p.verified_checkpoint(pd.DataFrame([{**row, **changed}]), current).empty)
        with self.assertWarns(UserWarning):
            self.assertTrue(p.verified_checkpoint(current.drop(columns="scoring_signature"), current).empty)
        self.assertNotEqual(p.stable_fold_task_id("t", "ACDE", "model1"), p.stable_fold_task_id("t", "ACDE", "model2"))

    def test_cache_reuse_and_invalidation_without_network(self):
        pdb = self.write_table("cache.pdb", p.build_ca_pdb_from_coordinates("ACDE", atom37()))
        metadata = {**p.prediction_provenance("ACDE", "model"), "pdb_sha256": p.sha256_file(pdb)}
        pdb.with_suffix(".json").write_text(json.dumps(metadata))
        class NoNetwork:
            def post(self, *a, **kw):
                raise AssertionError("Network forbidden in offline test")
        async def fold(sequence="ACDE", model="model"):
            return await p.fold_sequence(NoNetwork(), asyncio.Semaphore(1), "fake-key", sequence, model, pdb, 1, 1, 0, self.root / "attempts.jsonl")
        self.assertEqual(asyncio.run(fold())[1], 0)
        for sequence, model in (("ACDF", "model"), ("ACDE", "changed-model")):
            with self.assertWarns(UserWarning):
                self.assertIsNone(asyncio.run(fold(sequence, model))[0])
        pdb.write_text(pdb.read_text() + "REMARK changed\n")
        with self.assertWarns(UserWarning):
            self.assertIsNone(asyncio.run(fold())[0])

    def test_mocked_api_conversion_and_response_model(self):
        """Synthetic API payload tests wire behavior, not prediction quality."""
        class Response:
            status = 200
            async def __aenter__(self):
                return self
            async def __aexit__(self, *args):
                pass
            async def json(self):
                return {"model": self.model, "coordinates": atom37()}
        class Session:
            def post(self, url, **kwargs):
                self.request = (url, kwargs)
                response = Response()
                response.model = self.response_model
                return response
        session = Session()
        session.response_model = "OFFLINE-FIXTURE"
        pred = self.root / "fixture.pdb"
        async def fold():
            return await p.fold_sequence(session, asyncio.Semaphore(1), "OFFLINE-FAKE-KEY", "ACDE", "OFFLINE-FIXTURE", pred, 1, 1, 0, self.root / "attempts.jsonl")
        self.assertEqual(asyncio.run(fold())[:2], (pred, 1))
        self.assertEqual(session.request[0], p.OFFICIAL_FOLD_API_URL)
        self.assertEqual(session.request[1]["json"], {"sequence": "ACDE", "model": "OFFLINE-FIXTURE"})
        self.assertEqual(json.loads(pred.with_suffix(".json").read_text())["pdb_sha256"], p.sha256_file(pred))
        # Corrupt sidecars fail closed and do not crash before a request can be made.
        pred.with_suffix(".json").write_text("[]")
        session.response_model = "WRONG-MODEL"
        with self.assertWarns(UserWarning):
            path, _, error = asyncio.run(fold())
        self.assertIsNone(path)
        self.assertIn("differs from requested", error)

    def test_summary_duplicates_and_failures(self):
        frame = pd.DataFrame([
            {"target_name": "tem", "fold_task_id": "a", "status": "scored", "tm_score": .8},
            {"target_name": "tem", "fold_task_id": "a", "status": "scored", "tm_score": .8},
            {"target_name": "tem", "fold_task_id": "b", "status": "fold_failed", "tm_score": None}])
        summary = p.summarize_results(frame)
        self.assertEqual((summary["n_rows"], summary["n_unique_tasks"], summary["n_scored_rows"], summary["n_scored_unique_tasks"]), (3, 2, 2, 1))
        self.assertEqual(summary["status_counts"]["fold_failed"], 1)

    @unittest.skipUnless(shutil.which("TMalign"), "optional local TMalign integration")
    def test_actual_tmalign_self_alignment(self):
        target = p.prepare_template_pdb(p.TargetDefinition("tem", p.ROOT / "TEM_pdb.ent", "A"), self.root)
        result = p.run_tmalign(target.template_path, target.template_path)
        self.assertEqual(result["tm_score"], 1.)
        self.assertEqual(result["wt_length"], 263)


if __name__ == "__main__":
    unittest.main()
