import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import batch_evaluate


class CommandLogRetentionTest(unittest.TestCase):
    def run_with_result(
        self, run_dir: Path, returncode: int, stdout: str, stderr: str
    ) -> Path:
        command_log = run_dir / "commands.txt"
        completed = batch_evaluate.subprocess.CompletedProcess(
            ["example-tool"], returncode, stdout, stderr
        )
        with mock.patch.object(batch_evaluate.subprocess, "run", return_value=completed):
            batch_evaluate.run_command(
                ["example-tool", "--flag"], run_dir, "example", command_log
            )
        return command_log

    def test_success_without_stderr_keeps_only_commands(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run_dir = Path(temporary)
            (run_dir / "example.stdout.txt").write_text("stale", encoding="utf-8")
            (run_dir / "example.stderr.txt").write_text("stale", encoding="utf-8")

            command_log = self.run_with_result(run_dir, 0, "routine output\n", "")

            self.assertTrue(command_log.is_file())
            self.assertIn("example-tool --flag", command_log.read_text(encoding="utf-8"))
            self.assertIn("exit_code=0", command_log.read_text(encoding="utf-8"))
            self.assertFalse((run_dir / "example.stdout.txt").exists())
            self.assertFalse((run_dir / "example.stderr.txt").exists())

    def test_success_keeps_nonempty_stderr_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run_dir = Path(temporary)

            self.run_with_result(run_dir, 0, "routine output\n", "warning\n")

            self.assertFalse((run_dir / "example.stdout.txt").exists())
            self.assertEqual(
                (run_dir / "example.stderr.txt").read_text(encoding="utf-8"),
                "warning\n",
            )

    def test_failure_keeps_both_streams_and_commands(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run_dir = Path(temporary)

            command_log = self.run_with_result(run_dir, 7, "details\n", "failure\n")

            self.assertTrue(command_log.is_file())
            self.assertEqual(
                (run_dir / "example.stdout.txt").read_text(encoding="utf-8"),
                "details\n",
            )
            self.assertEqual(
                (run_dir / "example.stderr.txt").read_text(encoding="utf-8"),
                "failure\n",
            )
            self.assertIn("exit_code=7", command_log.read_text(encoding="utf-8"))


class BatchManifestTest(unittest.TestCase):
    def test_one_manifest_is_shared_by_all_predictions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            evaluation_dir = Path(temporary) / "evaluation"
            annotation_path = evaluation_dir / "annotations" / "samples.json"
            annotations = {"schema_version": 1, "samples": {}}
            discovered_predictions = []
            audit_predictions = []

            for index in range(2):
                sample_id = f"sample_{index}"
                prediction_path = f"predicted/method/{sample_id}/model.cif"
                reference_path = f"references/{sample_id}.cif"
                prediction_file = evaluation_dir / prediction_path
                reference_file = evaluation_dir / reference_path
                prediction_file.parent.mkdir(parents=True, exist_ok=True)
                reference_file.parent.mkdir(parents=True, exist_ok=True)
                prediction_file.write_text(f"prediction {index}\n", encoding="utf-8")
                reference_file.write_text(f"reference {index}\n", encoding="utf-8")
                annotations["samples"][sample_id] = {
                    "reference": reference_path,
                    "target_chains": ["T"],
                    "binders": [{"id": "binder", "chains": ["B"]}],
                    "predictions": {
                        prediction_path: {
                            "model_chain_id_type": "label",
                            "chain_mapping": {"T": "A", "B": "B"},
                        }
                    },
                }
                discovered_predictions.append(prediction_path)
                audit_predictions.append(
                    {"prediction": prediction_path, "short_binder_chains": []}
                )

            annotation_path.parent.mkdir(parents=True, exist_ok=True)
            annotation_path.write_text(json.dumps(annotations), encoding="utf-8")
            report = {
                "discovered_predictions": discovered_predictions,
                "predictions": audit_predictions,
            }
            records = batch_evaluate.prediction_records(annotations)
            versions = {"openstructure": "OpenStructure 2.12.0"}

            manifest_path, identities = batch_evaluate.write_batch_manifest(
                evaluation_dir,
                annotation_path,
                annotations,
                report,
                records,
                versions,
            )

            self.assertEqual(set(identities), set(discovered_predictions))
            self.assertEqual(
                manifest_path.parent.parent,
                evaluation_dir / "results" / "batches",
            )
            self.assertEqual(list(evaluation_dir.rglob("run_manifest.json")), [manifest_path])
            manifest = batch_evaluate.load_json(manifest_path)
            self.assertEqual(len(manifest["inputs"]["predictions"]), 2)

            reused_path, reused_identities = batch_evaluate.write_batch_manifest(
                evaluation_dir,
                annotation_path,
                annotations,
                report,
                records,
                versions,
            )
            self.assertEqual(reused_path, manifest_path)
            self.assertEqual(set(reused_identities), set(discovered_predictions))
            self.assertEqual(list(evaluation_dir.rglob("run_manifest.json")), [manifest_path])


if __name__ == "__main__":
    unittest.main()
