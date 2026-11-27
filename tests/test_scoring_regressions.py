import csv
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import batch_evaluate as batch
import chirality_angles as chirality
import summarize_batch_metrics as reporting
from target_aligned_binder_rmsd import compute


def write_structure(path, chains):
    path.parent.mkdir(parents=True, exist_ok=True)
    structure = gemmi.Structure()
    model = gemmi.Model("1")
    for auth, label, codes in chains:
        chain = gemmi.Chain(auth)
        for i, code in enumerate(codes, 1):
            residue = gemmi.Residue()
            residue.name = code
            residue.seqid = gemmi.SeqId(i, " ")
            residue.label_seq = i
            residue.subchain = label
            for name, (x, y, z) in {
                "N": (1, 0, 0), "CA": (0, 0, 0), "C": (0, 1, 0),
                "O": (0, 2, 0), "CB": (0, 0, 1),
            }.items():
                atom = gemmi.Atom()
                atom.name = name
                atom.element = gemmi.Element(name[0])
                atom.occ = 1
                atom.pos = gemmi.Position(x + 4 * i, y, z)
                residue.add_atom(atom)
            chain.add_residue(residue)
        model.add_chain(chain)
    structure.add_model(model)
    structure.make_mmcif_document().write_file(str(path))


class ScoringRegressions(unittest.TestCase):
    def test_chirality_uses_the_declared_namespace_and_inclusive_threshold(self):
        with tempfile.TemporaryDirectory() as temp:
            model = Path(temp) / "model.cif"
            write_structure(model, [("A", "B", ["ALA"]), ("B", "A", ["DAL"])])
            rows, report = chirality.evaluate(model, include_laa=True, chains=["A"], chain_id_type="label")
            self.assertEqual([(row["chain"], row["target_dl"]) for row in rows], [("A", "D")])
            self.assertEqual(report["by_chain"]["A"]["DAA"]["n_assessable"], 1)
            rows, _ = chirality.evaluate(model, include_laa=True, chains=["A"], chain_id_type="auth")
            self.assertEqual([row["target_dl"] for row in rows], ["L"])
        self.assertTrue(chirality.angle_metrics(math.radians(15), True, 15)["passed"])

    def test_binder_folding_survives_missing_target_correspondence(self):
        with tempfile.TemporaryDirectory() as temp:
            reference, model = Path(temp) / "ref.cif", Path(temp) / "mdl.cif"
            write_structure(reference, [("T", "T", ["ALA"]), ("B", "B", ["ALA"] * 3)])
            write_structure(model, [("Y", "Y", ["ALA"]), ("X", "X", ["ALA"] * 2)])
            result = compute({
                "reference": str(reference), "model": str(model),
                "reference_chain_id_type": "label", "model_chain_id_type": "label",
                "chain_mapping": {"T": "Y", "B": "X"}, "target_chains": ["T"],
                "binders": [{"id": "binder", "chains": ["B"]}],
            }, {"status": "SUCCESS", "aln": [">reference:B\nAAA\n>model:X\nAA-"]})
            self.assertIsNone(result["target_aligned"]["alignment"]["target_backbone_rmsd_angstrom"])
            placed = result["target_aligned"]["binders"][0]["target_aligned_binder_backbone_rmsd_angstrom"]
            self.assertEqual(placed["status"], "UNAVAILABLE")
            binder = result["self_aligned"]["binders"][0]["binder_self_aligned_backbone_rmsd_angstrom"]
            self.assertAlmostEqual(binder["rmsd_angstrom"], 0)
            self.assertAlmostEqual(binder["reference_residue_coverage"], 2 / 3)
            self.assertEqual(binder["atom_pair_coverage"], 1)
            self.assertEqual(result["residue_correspondence"]["chains"]["B"]["unmatched_reference_residues"], ["3"])

    def test_retry_ignores_old_outputs_and_context_mappings(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            pred = "predicted/method/model.cif"
            sample = {
                "reference": "ref.cif", "reference_chain_id_type": "label", "reference_assembly": {},
                "target_chains": ["T"], "binders": [{"id": "binder", "chains": ["B"]}],
                "context_chains": ["C"], "chirality_chains": ["B"],
            }
            prediction = {"model_chain_id_type": "label", "chain_mapping": {"T": "T", "B": "B", "C": "C"}}
            chains = [(chain, chain, ["ALA"]) for chain in "TBC"]
            write_structure(root / pred, chains)
            write_structure(root / "ref.cif", chains)
            run_root = root / "results/method/model.cif/runs/current"
            fingerprint = "current-full-fingerprint"
            stale = run_root / "chirality/all_residue_chirality_summary.json"
            batch.write_json(stale, {"pooled": {"all_selected": {"n_sites": 1, "n_assessable": 1}}})
            batch.write_json(run_root / "summary.json", {"status": "PARTIAL", "run_fingerprint": fingerprint})
            batch.write_skipped_summary(root, "sample", pred, {"skip_reason": "old skip"}, {})
            identity = {
                "run_dir": run_root, "fingerprint": fingerprint,
                "settings": {"ost": {"min_pep_length": 6}, "chirality": {"scope": "whole-example"}},
            }
            with mock.patch.object(batch, "run_command", return_value=subprocess.CompletedProcess([], 1)) as run:
                summary = batch.run_one(root, root / "annotations.json", "sample", sample, pred,
                                        prediction, {"warnings": [], "short_binder_chains": []},
                                        Path("chirality"), Path("rmsd"), Path("role"), {}, True, identity)
            commands = {call.args[2]: call.args[0] for call in run.call_args_list}
            self.assertNotIn("C:C", commands["openstructure"])
            self.assertEqual(summary["chain_mapping_used"], {"T": "T", "B": "B"})
            self.assertIn("--chain-id-type", commands["chirality"])
            selected = [value for flag, value in zip(commands["chirality"], commands["chirality"][1:])
                        if flag == "--chain"]
            self.assertEqual(set(selected), {"T", "B"})
            actual = Path(summary["run_dir"])
            self.assertEqual(actual.parent, (run_root / "attempts").resolve())
            self.assertEqual(batch.chirality_summary(stale, 1)["status"], "UNAVAILABLE")
            self.assertEqual(summary["chirality"]["status"], "UNAVAILABLE")
            self.assertTrue(stale.exists())
            csv_path = batch.regenerate_batch_csv(root, [actual / "summary.json"])
            with csv_path.open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual([row["status"] for row in rows], ["FAILURE"])
            self.assertEqual(rows[0]["summary_path"], str(actual / "summary.json"))

    def test_batch_counts_pool_assessable_residues_once_per_prediction(self):
        first = {
            "prediction_path": "predicted/method/one.cif", "status": "SUCCESS",
            "chirality_all_correct_count": "2", "chirality_all_assessable_count": "2",
            "chirality_all_skipped_count": "1", "chirality_daa_correct_count": "1",
            "chirality_daa_assessable_count": "1", "chirality_laa_correct_count": "1",
            "chirality_laa_assessable_count": "1",
        }
        second = {
            "prediction_path": "predicted/method/two.cif", "status": "SUCCESS",
            "chirality_all_correct_count": "1", "chirality_all_eligible_count": "10",
            "chirality_all_skipped_count": "2", "chirality_daa_correct_count": "0",
            "chirality_daa_eligible_count": "6", "chirality_daa_skipped_count": "1",
            "chirality_laa_correct_count": "1", "chirality_laa_eligible_count": "4",
            "chirality_laa_skipped_count": "1",
        }
        report = reporting.summarize([first, first, second], Path("batch_summary.csv"))
        accuracies = report["methods"][0]["folding"]["chirality_accuracy"]
        self.assertEqual(accuracies["all"]["accuracy_pct"], 30)
        self.assertEqual(accuracies["daa"]["n_assessable"], 6)
        self.assertEqual(accuracies["laa"]["accuracy_pct"], 50)
        self.assertIn("3/10 (30.0%)", reporting.markdown(report))


if __name__ == "__main__":
    unittest.main()
