import json
from pathlib import Path
import sys
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import batch_evaluate
import summarize_batch_metrics


class RoleScopedInterfaceTest(unittest.TestCase):
    def test_excludes_target_target_and_binder_binder_interfaces(self) -> None:
        identities = [
            ["T1", "B1", "A", "C"],
            ["T2", "B2", "B", "D"],
            ["T1", "T2", "A", "B"],
            ["B1", "B2", "C", "D"],
        ]
        raw = {
            "dockq_reference_interfaces": [pair[:2] for pair in identities],
            "per_interface": [
                {"identity": pair, "dockq": score, "irmsd": score, "nnat": contacts}
                for pair, score, contacts in zip(
                    identities, [0.4, 0.1, 0.9, 0.8], [20, 80, 200, 200]
                )
            ],
        }
        role = batch_evaluate.target_binder_dockq_metrics(
            raw, ["T1", "T2"],
            [{"id": "one", "chains": ["B1"]}, {"id": "two", "chains": ["B2"]}],
        )
        self.assertEqual(len(role["per_interface"]), 2)
        self.assertAlmostEqual(role["mean"], 0.25)
        self.assertAlmostEqual(role["wave"], 0.16)
        threshold = batch_evaluate.evaluate_docking_thresholds(
            {"dockq_ave": 0.55}, {"binders": []}, 2.0, 0.23,
            {"target_binder_dockq": role},
        )
        self.assertEqual(threshold["target_binder_dockq_success_fraction"], 0.5)
        self.assertTrue(threshold["dockq_ave_above_threshold"])


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
            {"target_binder_dockq": {"per_interface": [
                {"identity": ["T", "B", "M", "N"], "dockq": dockq}
            ]}},
        )

    def test_cutoffs_are_strict(self) -> None:
        at_cutoff = self.classify(0.23, 2.0)
        self.assertFalse(at_cutoff["dockq_ave_above_threshold"])
        self.assertEqual(at_cutoff["target_binder_dockq_success_fraction"], 0.0)
        self.assertFalse(
            at_cutoff["binders"][0][
                "target_aligned_binder_backbone_rmsd_below_threshold"
            ]
        )

        passing = self.classify(0.230001, 1.999999)
        self.assertTrue(passing["dockq_ave_above_threshold"])
        self.assertEqual(passing["target_binder_dockq_success_fraction"], 1.0)
        self.assertTrue(
            passing["binders"][0][
                "target_aligned_binder_backbone_rmsd_below_threshold"
            ]
        )

    def test_unavailable_scores_are_not_failures(self) -> None:
        result = self.classify(None, None)
        self.assertIsNone(result["dockq_ave_above_threshold"])
        self.assertIsNone(result["target_binder_dockq_success_fraction"])
        self.assertIsNone(
            result["binders"][0][
                "target_aligned_binder_backbone_rmsd_below_threshold"
            ]
        )


class FoldingThresholdTest(unittest.TestCase):
    def classify(self, target_rmsd, binder_rmsd):
        return batch_evaluate.evaluate_folding_thresholds(
            {"alignment": {"target_backbone_rmsd_angstrom": target_rmsd}},
            {
                "binders": [
                    {
                        "id": "binder",
                        "binder_self_aligned_backbone_rmsd_angstrom": {
                            "rmsd_angstrom": binder_rmsd
                        },
                    }
                ]
            },
            2.0,
            2.0,
        )

    def test_cutoffs_are_strict(self) -> None:
        at_cutoff = self.classify(2.0, 2.0)
        self.assertFalse(
            at_cutoff["target_self_aligned_backbone_rmsd_below_threshold"]
        )
        self.assertFalse(
            at_cutoff["binders"][0][
                "binder_self_aligned_backbone_rmsd_below_threshold"
            ]
        )

        passing = self.classify(1.999999, 1.999999)
        self.assertTrue(
            passing["target_self_aligned_backbone_rmsd_below_threshold"]
        )
        self.assertTrue(
            passing["binders"][0][
                "binder_self_aligned_backbone_rmsd_below_threshold"
            ]
        )

    def test_unavailable_scores_are_not_failures(self) -> None:
        result = self.classify(None, None)
        self.assertIsNone(
            result["target_self_aligned_backbone_rmsd_below_threshold"]
        )
        self.assertIsNone(
            result["binders"][0][
                "binder_self_aligned_backbone_rmsd_below_threshold"
            ]
        )


class MethodSummaryThresholdTest(unittest.TestCase):
    def test_complex_balanced_interface_success_and_role_lddt(self) -> None:
        rows = [
            {
                "status": "SUCCESS",
                "prediction_path": "predicted/method/one.cif",
                "target_only_lddt": "0.8",
                "binder_only_lddt_complex_mean": "0.6",
                "ost_lddt": "0.9",
                "target_binder_dockq_mean": "0.25",
                "target_binder_dockq_wave": "0.16",
                "target_binder_dockq_success_fraction": "0.5",
                "target_binder_dockq_interface_success_count": "1",
                "target_binder_dockq_interface_available_count": "2",
                "target_binder_dockq_per_interface_json": json.dumps([
                    {"identity": ["T1", "B1", "A", "C"], "irmsd": 1.0},
                    {"identity": ["T2", "B2", "B", "D"], "irmsd": 3.0},
                ]),
                "target_binder_ilddt_macro": "0.7",
                "target_binder_ilddt_micro": "0.65",
            },
            {
                "status": "SUCCESS",
                "prediction_path": "predicted/method/two.cif",
                "target_only_lddt": "1.0",
                "binder_only_lddt_complex_mean": "0.8",
                "ost_lddt": "0.95",
                "target_binder_dockq_mean": "0.9",
                "target_binder_dockq_wave": "0.9",
                "target_binder_dockq_success_fraction": "1.0",
                "target_binder_dockq_interface_success_count": "1",
                "target_binder_dockq_interface_available_count": "1",
                "target_binder_dockq_per_interface_json": json.dumps([
                    {"identity": ["T", "B", "E", "F"], "irmsd": 5.0}
                ]),
                "target_binder_ilddt_macro": "0.9",
                "target_binder_ilddt_micro": "0.85",
            },
        ]
        report = summarize_batch_metrics.summarize(rows, Path("batch.csv"))
        method = report["methods"][0]
        self.assertEqual(method["docking"]["dockq_above_threshold"]["success_rate_pct"], 75.0)
        self.assertAlmostEqual(method["docking"]["dockq_mean"]["value"], 0.575)
        self.assertEqual(method["docking"]["irmsd_angstrom_median"]["value"], 3.5)
        self.assertEqual(method["supplemental_pooled"]["irmsd_angstrom_median"]["value"], 3.0)
        self.assertAlmostEqual(method["docking"]["ilddt_macro_mean"]["value"], 0.8)
        self.assertAlmostEqual(method["folding"]["target_only_lddt_mean"]["value"], 0.9)
        self.assertAlmostEqual(method["folding"]["binder_only_lddt_mean"]["value"], 0.7)
        rendered = summarize_batch_metrics.markdown(report)
        self.assertIn("Mean target-only lDDT", rendered)
        self.assertIn("75.000% (2 complexes; 2/3 interfaces)", rendered)
        self.assertIn("Supplemental metrics", rendered)

    def test_rates_use_prediction_denominators_with_group_fractions(self) -> None:
        rows = [
            {
                "status": "SUCCESS",
                "prediction_path": "predicted/method/sample-a/model.cif",
                "ost_dockq_ave": "0.24",
                "target_binder_dockq_mean": "0.24",
                "target_binder_dockq_success_fraction": "1.0",
                "target_binder_dockq_interface_success_count": "1",
                "target_binder_dockq_interface_available_count": "1",
                "target_backbone_fit_rmsd_angstrom": "1.9",
                "target_aligned_binder_backbone_rmsd_angstrom": "1.9",
                "binder_self_aligned_backbone_rmsd_angstrom": "1.9",
            },
            {
                "status": "SUCCESS",
                "prediction_path": "predicted/method/sample-a/model.cif",
                "ost_dockq_ave": "0.24",
                "target_binder_dockq_mean": "0.24",
                "target_binder_dockq_success_fraction": "1.0",
                "target_binder_dockq_interface_success_count": "1",
                "target_binder_dockq_interface_available_count": "1",
                "target_backbone_fit_rmsd_angstrom": "1.9",
                "target_aligned_binder_backbone_rmsd_angstrom": "2.0",
                "binder_self_aligned_backbone_rmsd_angstrom": "2.0",
            },
            {
                "status": "PARTIAL",
                "prediction_path": "predicted/method/sample-b/model.cif",
                "ost_dockq_ave": "0.23",
                "target_binder_dockq_mean": "0.23",
                "target_binder_dockq_success_fraction": "0.0",
                "target_binder_dockq_interface_success_count": "0",
                "target_binder_dockq_interface_available_count": "1",
                "target_backbone_fit_rmsd_angstrom": "2.0",
                "target_aligned_binder_backbone_rmsd_angstrom": "",
                "binder_self_aligned_backbone_rmsd_angstrom": "",
            },
        ]

        report = summarize_batch_metrics.summarize(
            rows, Path("batch.csv"), 2.0, 0.23, 2.0, 2.0
        )
        folding = report["methods"][0]["folding"]
        docking = report["methods"][0]["docking"]

        self.assertEqual(
            folding["target_self_aligned_backbone_rmsd_below_threshold"],
            {
                "threshold_angstrom": 2.0,
                "operator": "<",
                "n_success": 1,
                "n_available": 2,
                "success_rate_pct": 50.0,
            },
        )
        self.assertEqual(
            folding["binder_self_aligned_backbone_rmsd_below_threshold"],
            {
                "threshold_angstrom": 2.0,
                "operator": "<",
                "unit": "mean of per-complex available binder-group success fractions",
                "n_available": 1,
                "n_success_groups": 1,
                "n_available_groups": 2,
                "success_rate_pct": 50.0,
            },
        )

        self.assertEqual(
            docking["dockq_above_threshold"],
            {
                "threshold": 0.23,
                "operator": ">",
                "n_available": 2,
                "unit": "mean of per-complex target-binder interface success fractions",
                "n_success_interfaces": 1,
                "n_available_interfaces": 2,
                "success_rate_pct": 50.0,
            },
        )
        self.assertEqual(
            docking["target_aligned_binder_backbone_rmsd_below_threshold"],
            {
                "threshold_angstrom": 2.0,
                "operator": "<",
                "unit": "mean of per-complex available binder-group success fractions",
                "n_available": 1,
                "n_success_groups": 1,
                "n_available_groups": 2,
                "success_rate_pct": 50.0,
            },
        )

    def test_complex_weighting_changes_rmsds_but_not_chirality_pooling(self) -> None:
        rows = []
        for path, group_values, interface_values, chirality in (
            ("predicted/method/a.cif", [1.0, 1.0, 9.0], [1.0, 1.0, 9.0], (1, 2)),
            ("predicted/method/b.cif", [5.0], [5.0], (99, 100)),
        ):
            for index, value in enumerate(group_values):
                rows.append({
                    "status": "SUCCESS",
                    "prediction_path": path,
                    "binder_id": f"binder_{index}",
                    "binder_self_aligned_backbone_rmsd_angstrom": str(value),
                    "target_aligned_binder_backbone_rmsd_angstrom": str(value),
                    "target_binder_dockq_per_interface_json": json.dumps(
                        [{"irmsd": score} for score in interface_values]
                    ),
                    "chirality_all_correct_count": str(chirality[0]),
                    "chirality_all_eligible_count": str(chirality[1]),
                    "chirality_all_skipped_count": "0",
                })

        report = summarize_batch_metrics.summarize(rows, Path("batch.csv"))
        method = report["methods"][0]
        folding = method["folding"]
        docking = method["docking"]
        pooled = method["supplemental_pooled"]

        for section, key in (
            (folding, "binder_self_aligned_backbone_rmsd_angstrom_median"),
            (docking, "target_aligned_binder_backbone_rmsd_angstrom_median"),
            (docking, "irmsd_angstrom_median"),
        ):
            self.assertAlmostEqual(section[key]["value"], 13 / 3)
            self.assertEqual(section[key]["n"], 2)
            self.assertEqual(pooled[key]["value"], 3.0)
            self.assertEqual(pooled[key]["n"], 4)

        for section, key in (
            (folding, "binder_self_aligned_backbone_rmsd_below_threshold"),
            (docking, "target_aligned_binder_backbone_rmsd_below_threshold"),
        ):
            self.assertAlmostEqual(section[key]["success_rate_pct"], 100 / 3)
            self.assertEqual(section[key]["n_available"], 2)
            self.assertEqual(section[key]["n_success_groups"], 2)
            self.assertEqual(section[key]["n_available_groups"], 4)

        self.assertEqual(folding["chirality_violation"]["n_violations"], 2)
        self.assertEqual(folding["chirality_violation"]["n_assessable"], 102)
        self.assertAlmostEqual(folding["chirality_violation"]["violation_rate_pct"], 200 / 102)
        rendered = summarize_batch_metrics.markdown(report)
        self.assertIn("median of complex means", rendered)
        self.assertIn("33.333% (2 complexes; 2/4 groups)", rendered)


if __name__ == "__main__":
    unittest.main()
