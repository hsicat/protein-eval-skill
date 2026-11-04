from pathlib import Path
import sys
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import batch_evaluate
import summarize_batch_metrics


class PerPredictionThresholdTest(unittest.TestCase):
    def classify(self, dockq, rmsd):
        return batch_evaluate.evaluate_docking_thresholds(
            {"dockq_ave": dockq},
            {
                "binders": [
                    {
                        "id": "binder",
                        "target_aligned_binder_backbone_rmsd_angstrom": {
                            "rmsd_angstrom": rmsd
                        },
                    }
                ]
            },
            2.0,
            0.23,
        )

    def test_cutoffs_are_strict(self) -> None:
        at_cutoff = self.classify(0.23, 2.0)
        self.assertFalse(at_cutoff["dockq_ave_above_threshold"])
        self.assertFalse(
            at_cutoff["binders"][0][
                "target_aligned_binder_backbone_rmsd_below_threshold"
            ]
        )

        passing = self.classify(0.230001, 1.999999)
        self.assertTrue(passing["dockq_ave_above_threshold"])
        self.assertTrue(
            passing["binders"][0][
                "target_aligned_binder_backbone_rmsd_below_threshold"
            ]
        )

    def test_unavailable_scores_are_not_failures(self) -> None:
        result = self.classify(None, None)
        self.assertIsNone(result["dockq_ave_above_threshold"])
        self.assertIsNone(
            result["binders"][0][
                "target_aligned_binder_backbone_rmsd_below_threshold"
            ]
        )


class MethodSummaryThresholdTest(unittest.TestCase):
    def test_rates_use_prediction_and_binder_group_denominators(self) -> None:
        rows = [
            {
                "status": "SUCCESS",
                "prediction_path": "predicted/method/sample-a/model.cif",
                "ost_dockq_ave": "0.24",
                "target_aligned_binder_backbone_rmsd_angstrom": "1.9",
            },
            {
                "status": "SUCCESS",
                "prediction_path": "predicted/method/sample-a/model.cif",
                "ost_dockq_ave": "0.24",
                "target_aligned_binder_backbone_rmsd_angstrom": "2.0",
            },
            {
                "status": "PARTIAL",
                "prediction_path": "predicted/method/sample-b/model.cif",
                "ost_dockq_ave": "0.23",
                "target_aligned_binder_backbone_rmsd_angstrom": "",
            },
        ]

        report = summarize_batch_metrics.summarize(rows, Path("batch.csv"), 2.0, 0.23)
        docking = report["methods"][0]["docking"]

        self.assertEqual(
            docking["dockq_above_threshold"],
            {
                "threshold": 0.23,
                "operator": ">",
                "n_success": 1,
                "n_available": 2,
                "success_rate_pct": 50.0,
            },
        )
        self.assertEqual(
            docking["target_aligned_binder_backbone_rmsd_below_threshold"],
            {
                "threshold_angstrom": 2.0,
                "operator": "<",
                "n_success": 1,
                "n_available": 2,
                "success_rate_pct": 50.0,
            },
        )


if __name__ == "__main__":
    unittest.main()
