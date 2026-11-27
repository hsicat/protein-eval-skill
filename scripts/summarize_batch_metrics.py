#!/usr/bin/env python3
"""Create deterministic method-level Folding and Docking summary tables."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path, PurePosixPath
import statistics
from typing import Any, Callable


DEFAULT_BINDER_RMSD_THRESHOLD_ANGSTROM = 2.0
DEFAULT_DOCKQ_THRESHOLD = 0.23
DEFAULT_TARGET_SELF_ALIGNED_RMSD_THRESHOLD_ANGSTROM = 2.0
DEFAULT_BINDER_SELF_ALIGNED_RMSD_THRESHOLD_ANGSTROM = 2.0


def positive_finite_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be a finite number greater than zero")
    return number


def dockq_threshold_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or not 0 <= number <= 1:
        raise argparse.ArgumentTypeError("must be a finite number between zero and one")
    return number


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch_csv", type=Path)
    parser.add_argument("--method-name", help="Use one method label for all rows, e.g. paths organized by seed")
    parser.add_argument("--json-output", type=Path, required=True)
    parser.add_argument("--markdown-output", type=Path, required=True)
    parser.add_argument(
        "--binder-rmsd-threshold-angstrom",
        type=positive_finite_float,
        help=(
            "Strict upper cutoff; otherwise infer it from the CSV or use "
            f"{DEFAULT_BINDER_RMSD_THRESHOLD_ANGSTROM:g} Angstrom"
        ),
    )
    parser.add_argument(
        "--dockq-threshold",
        type=dockq_threshold_float,
        help=(
            "Strict lower cutoff applied to saved target-binder interface scores; "
            "otherwise infer it from the CSV or use "
            f"{DEFAULT_DOCKQ_THRESHOLD:g}"
        ),
    )
    parser.add_argument(
        "--target-self-aligned-rmsd-threshold-angstrom",
        type=positive_finite_float,
        help=(
            "Strict upper cutoff; otherwise infer it from the CSV or use "
            f"{DEFAULT_TARGET_SELF_ALIGNED_RMSD_THRESHOLD_ANGSTROM:g} Angstrom"
        ),
    )
    parser.add_argument(
        "--binder-self-aligned-rmsd-threshold-angstrom",
        type=positive_finite_float,
        help=(
            "Strict upper cutoff; otherwise infer it from the CSV or use "
            f"{DEFAULT_BINDER_SELF_ALIGNED_RMSD_THRESHOLD_ANGSTROM:g} Angstrom"
        ),
    )
    return parser.parse_args()


def finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def integer(value: Any) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def pooled_chirality(rows: list[dict[str, Any]], group: str) -> dict[str, Any]:
    prefix = f"chirality_{group}_"
    passed = sum(integer(row.get(prefix + "correct_count")) for row in rows)
    unassessable = sum(integer(row.get(prefix + "skipped_count")) for row in rows)
    assessable = sum(
        integer(row[prefix + "assessable_count"])
        if row.get(prefix + "assessable_count") not in (None, "")
        else integer(row.get(prefix + "eligible_count")) - integer(row.get(prefix + "skipped_count"))
        for row in rows
    )
    return {
        "n_pass": passed, "n_assessable": assessable, "n_unassessable": unassessable,
        "accuracy_pct": 100 * passed / assessable if assessable else None,
    }


def method_from_prediction(path: str) -> str:
    parts = PurePosixPath(path).parts
    if len(parts) >= 2 and parts[0] == "predicted":
        return parts[1]
    return parts[0] if parts else "unknown"


def aggregate(values: list[float], operation: Callable[[list[float]], float]) -> dict[str, Any]:
    return {"value": operation(values) if values else None, "n": len(values)}


def lddt_mode_means(
    rows: list[dict[str, str]], columns: dict[str, str]
) -> dict[str, Any]:
    return {
        f"{metric}{suffix}": aggregate(
            [
                score for row in rows
                if (score := finite_float(row.get(f"{column}{suffix}"))) is not None
            ],
            statistics.mean,
        )
        for metric, column in columns.items()
        for suffix in ("", "_no_stereochecks")
    }


def success_rate(values: list[float], predicate: Callable[[float], bool]) -> dict[str, Any]:
    n_success = sum(predicate(value) for value in values)
    return {
        "n_success": n_success,
        "n_available": len(values),
        "success_rate_pct": 100.0 * n_success / len(values) if values else None,
    }


def complex_group_values(
    by_prediction: dict[str, list[dict[str, str]]], column: str, threshold: float
) -> tuple[list[float], list[float], list[float]]:
    """Return prediction means, prediction pass fractions, and original pooled values."""
    complex_means: list[float] = []
    complex_pass_fractions: list[float] = []
    pooled_values: list[float] = []
    for rows in by_prediction.values():
        values = [
            value for row in rows
            if (value := finite_float(row.get(column))) is not None
        ]
        pooled_values.extend(values)
        if values:
            complex_means.append(statistics.mean(values))
            complex_pass_fractions.append(
                sum(value < threshold for value in values) / len(values)
            )
    return complex_means, complex_pass_fractions, pooled_values


def complex_fraction_success(
    fractions: list[float], pooled_values: list[float], threshold: float
) -> dict[str, Any]:
    return {
        "threshold_angstrom": threshold,
        "operator": "<",
        "unit": "mean of per-complex available binder-group success fractions",
        "n_available": len(fractions),
        "n_success_groups": sum(value < threshold for value in pooled_values),
        "n_available_groups": len(pooled_values),
        "success_rate_pct": (
            100.0 * statistics.mean(fractions) if fractions else None
        ),
    }


def resolve_threshold(
    rows: list[dict[str, str]],
    column: str,
    explicit: float | None,
    default: float,
) -> float:
    if explicit is not None:
        return explicit
    values = {
        value
        for row in rows
        if (value := finite_float(row.get(column))) is not None
    }
    if len(values) > 1:
        raise ValueError(f"Mixed values in {column}: {sorted(values)}")
    return next(iter(values), default)


def load_interfaces(raw: str) -> list[dict[str, Any]]:
    if not raw:
        return []
    try:
        value = json.loads(raw)
    except json.JSONDecodeError:
        return []
    return value if isinstance(value, list) else []


def dockq_interface_pass_counts(
    scores: list[float], row: dict[str, str], threshold: float
) -> tuple[int, int]:
    if scores:
        return sum(score > threshold for score in scores), len(scores)

    # Older CSVs can have counts without individual scores.
    recorded_cutoff = row.get("target_binder_dockq_threshold")
    stored_threshold = (
        DEFAULT_DOCKQ_THRESHOLD if recorded_cutoff in (None, "")
        else finite_float(recorded_cutoff)
    )
    if threshold != stored_threshold:
        return 0, 0
    return (
        integer(row.get("target_binder_dockq_interface_success_count")),
        integer(row.get("target_binder_dockq_interface_available_count")),
    )


def summarize(
    rows: list[dict[str, str]],
    source: Path,
    binder_rmsd_threshold_angstrom: float | None = None,
    dockq_threshold: float | None = None,
    target_self_aligned_rmsd_threshold_angstrom: float | None = None,
    binder_self_aligned_rmsd_threshold_angstrom: float | None = None,
    method_name: str | None = None,
) -> dict[str, Any]:
    eligible = [row for row in rows if row.get("status") != "SKIPPED"]
    binder_rmsd_threshold_angstrom = resolve_threshold(
        eligible,
        "target_aligned_binder_backbone_rmsd_threshold_angstrom",
        binder_rmsd_threshold_angstrom,
        DEFAULT_BINDER_RMSD_THRESHOLD_ANGSTROM,
    )
    dockq_threshold = resolve_threshold(
        eligible,
        "target_binder_dockq_threshold",
        dockq_threshold,
        DEFAULT_DOCKQ_THRESHOLD,
    )
    target_self_aligned_rmsd_threshold_angstrom = resolve_threshold(
        eligible,
        "target_self_aligned_backbone_rmsd_threshold_angstrom",
        target_self_aligned_rmsd_threshold_angstrom,
        DEFAULT_TARGET_SELF_ALIGNED_RMSD_THRESHOLD_ANGSTROM,
    )
    binder_self_aligned_rmsd_threshold_angstrom = resolve_threshold(
        eligible,
        "binder_self_aligned_backbone_rmsd_threshold_angstrom",
        binder_self_aligned_rmsd_threshold_angstrom,
        DEFAULT_BINDER_SELF_ALIGNED_RMSD_THRESHOLD_ANGSTROM,
    )
    by_method: dict[str, list[dict[str, str]]] = {}
    for row in eligible:
        by_method.setdefault(method_name or method_from_prediction(row.get("prediction_path", "")), []).append(row)

    methods: list[dict[str, Any]] = []
    for method, method_rows in sorted(by_method.items()):
        by_prediction: dict[str, list[dict[str, str]]] = {}
        for row in method_rows:
            by_prediction.setdefault(row["prediction_path"], []).append(row)
        prediction_rows = [group[0] for group in by_prediction.values()]

        target_fold = [
            value
            for row in prediction_rows
            if (value := finite_float(row.get("target_backbone_fit_rmsd_angstrom"))) is not None
        ]
        binder_fold, binder_fold_success, binder_fold_pooled = complex_group_values(
            by_prediction,
            "binder_self_aligned_backbone_rmsd_angstrom",
            binder_self_aligned_rmsd_threshold_angstrom,
        )
        chirality = {group: pooled_chirality(prediction_rows, group)
                     for group in ("all", "daa", "laa")}
        chirality_assessable = chirality["all"]["n_assessable"]
        chirality_violations = chirality_assessable - chirality["all"]["n_pass"]

        dockq = [
            value
            for row in prediction_rows
            if (value := finite_float(row.get("target_binder_dockq_mean"))) is not None
        ]
        dockq_wave = [
            value for row in prediction_rows
            if (value := finite_float(row.get("target_binder_dockq_wave"))) is not None
        ]
        dockq_success_fractions: list[float] = []
        dockq_interface_successes = 0
        dockq_interface_available = 0
        target_aligned_binder, target_aligned_success, target_aligned_pooled = (
            complex_group_values(
                by_prediction,
                "target_aligned_binder_backbone_rmsd_angstrom",
                binder_rmsd_threshold_angstrom,
            )
        )
        irmsd: list[float] = []
        irmsd_pooled: list[float] = []
        for row in prediction_rows:
            interfaces = load_interfaces(row.get("target_binder_dockq_per_interface_json", ""))
            dockq_scores = [
                score for interface in interfaces
                if (score := finite_float(interface.get("dockq"))) is not None
            ]
            passed, available = dockq_interface_pass_counts(dockq_scores, row, dockq_threshold)
            if available:
                dockq_success_fractions.append(passed / available)
                dockq_interface_successes += passed
                dockq_interface_available += available
            irmsd_scores = [
                score for interface in interfaces
                if (score := finite_float(interface.get("irmsd"))) is not None
            ]
            irmsd_pooled.extend(irmsd_scores)
            if irmsd_scores:
                irmsd.append(statistics.mean(irmsd_scores))

        methods.append(
            {
                "method": method,
                "n_predictions": len(by_prediction),
                "n_binder_groups": len(method_rows),
                "folding": {
                    "target_self_aligned_backbone_rmsd_angstrom_median": aggregate(
                        target_fold, statistics.median
                    ),
                    "target_self_aligned_backbone_rmsd_below_threshold": {
                        "threshold_angstrom": (
                            target_self_aligned_rmsd_threshold_angstrom
                        ),
                        "operator": "<",
                        **success_rate(
                            target_fold,
                            lambda value: (
                                value < target_self_aligned_rmsd_threshold_angstrom
                            ),
                        ),
                    },
                    "binder_self_aligned_backbone_rmsd_angstrom_median": aggregate(
                        binder_fold, statistics.median
                    ),
                    "binder_self_aligned_backbone_rmsd_below_threshold": complex_fraction_success(
                        binder_fold_success,
                        binder_fold_pooled,
                        binder_self_aligned_rmsd_threshold_angstrom,
                    ),
                    **lddt_mode_means(prediction_rows, {
                        "target_only_lddt_mean": "target_only_lddt",
                        "binder_only_lddt_mean": "binder_only_lddt_complex_mean",
                        "global_scored_subset_lddt_mean": "ost_lddt",
                    }),
                    "chirality_accuracy": chirality,
                    "chirality_violation": {
                        "n_violations": chirality_violations,
                        "n_assessable": chirality_assessable,
                        "n_unassessable": chirality["all"]["n_unassessable"],
                        "violation_rate_pct": (
                            100.0 * chirality_violations / chirality_assessable
                            if chirality_assessable
                            else None
                        ),
                    },
                },
                "docking": {
                    "dockq_mean": aggregate(dockq, statistics.mean),
                    "dockq_wave_mean": aggregate(dockq_wave, statistics.mean),
                    "dockq_above_threshold": {
                        "threshold": dockq_threshold,
                        "operator": ">",
                        "unit": "mean of per-complex target-binder interface success fractions",
                        "n_available": len(dockq_success_fractions),
                        "n_success_interfaces": dockq_interface_successes,
                        "n_available_interfaces": dockq_interface_available,
                        "success_rate_pct": (
                            100.0 * statistics.mean(dockq_success_fractions)
                            if dockq_success_fractions else None
                        ),
                    },
                    "target_aligned_binder_backbone_rmsd_angstrom_median": aggregate(
                        target_aligned_binder, statistics.median
                    ),
                    "target_aligned_binder_backbone_rmsd_below_threshold": complex_fraction_success(
                        target_aligned_success,
                        target_aligned_pooled,
                        binder_rmsd_threshold_angstrom,
                    ),
                    "irmsd_angstrom_median": aggregate(irmsd, statistics.median),
                    **lddt_mode_means(prediction_rows, {
                        "ilddt_macro_mean": "target_binder_ilddt_macro",
                        "ilddt_micro_mean": "target_binder_ilddt_micro",
                    }),
                },
                "supplemental_pooled": {
                    "binder_self_aligned_backbone_rmsd_angstrom_median": aggregate(
                        binder_fold_pooled, statistics.median
                    ),
                    "binder_self_aligned_backbone_rmsd_below_threshold": {
                        "threshold_angstrom": binder_self_aligned_rmsd_threshold_angstrom,
                        "operator": "<",
                        **success_rate(
                            binder_fold_pooled,
                            lambda value: value < binder_self_aligned_rmsd_threshold_angstrom,
                        ),
                    },
                    "target_aligned_binder_backbone_rmsd_angstrom_median": aggregate(
                        target_aligned_pooled, statistics.median
                    ),
                    "target_aligned_binder_backbone_rmsd_below_threshold": {
                        "threshold_angstrom": binder_rmsd_threshold_angstrom,
                        "operator": "<",
                        **success_rate(
                            target_aligned_pooled,
                            lambda value: value < binder_rmsd_threshold_angstrom,
                        ),
                    },
                    "irmsd_angstrom_median": aggregate(irmsd_pooled, statistics.median),
                },
            }
        )

    return {
        "schema_version": 6,
        "source_batch_csv": str(source.resolve()),
        "folding_success_thresholds": {
            "target_self_aligned_backbone_rmsd_angstrom": {
                "operator": "<",
                "value": target_self_aligned_rmsd_threshold_angstrom,
            },
            "binder_self_aligned_backbone_rmsd_angstrom": {
                "operator": "<",
                "value": binder_self_aligned_rmsd_threshold_angstrom,
            },
        },
        "docking_success_thresholds": {
            "target_aligned_binder_backbone_rmsd_angstrom": {
                "operator": "<",
                "value": binder_rmsd_threshold_angstrom,
            },
            "target_binder_interface_dockq": {"operator": ">", "value": dockq_threshold},
        },
        "aggregation": {
            "prediction_level_scores": "Deduplicated by prediction path before aggregation",
            "rmsd": (
                "Median across available prediction values; binder-group and "
                "target-binder-interface RMSDs are first averaged within each prediction"
            ),
            "bounded_scores": "Arithmetic mean across available prediction values",
            "folding_success": (
                "Independent strict-threshold success proportions; target RMSD is "
                "counted per prediction, binder-group RMSD is classified per group "
                "then averaged within and across predictions, and unavailable values "
                "are excluded from each denominator"
            ),
            "docking_success": (
                "DockQ threshold is applied per target-binder interface, then the "
                "success fraction is averaged within each prediction and across "
                "available predictions; binder RMSD threshold success uses the "
                "same group-then-prediction hierarchy independently"
            ),
            "within_complex": (
                "Binder-only lDDT and target-binder iLDDT macro are equal means "
                "of eligible binder groups and chain-pair interfaces respectively; "
                "iLDDT micro pools distance tests only within each prediction"
            ),
            "chirality_violation": (
                "Residue-pooled (inverted + below-threshold) / assessable residues; "
                "unassessable residues reported separately"
            ),
            "chirality_accuracy": (
                "Passed / assessable residues, pooled across included complex–seed predictions "
                "for all (D+L), D-only and L-only; the pass cutoff is at least 15 degrees by default"
            ),
        },
        "methods": methods,
    }


def format_number(value: Any, digits: int = 3) -> str:
    return "NA" if value is None else f"{value:.{digits}f}"


def format_chirality(value: dict[str, Any]) -> str:
    if value["accuracy_pct"] is None:
        return "NA"
    return f"{value['n_pass']}/{value['n_assessable']} ({value['accuracy_pct']:.1f}%)"


def format_success(value: dict[str, Any]) -> str:
    if not value["n_available"]:
        return "NA"
    return (
        f"{value['n_success']}/{value['n_available']} "
        f"({format_number(value['success_rate_pct'])}%)"
    )


def format_dockq_success(value: dict[str, Any]) -> str:
    if not value["n_available"]:
        return "NA"
    return (
        f"{format_number(value['success_rate_pct'])}% "
        f"({value['n_available']} predictions; "
        f"{value['n_success_interfaces']}/{value['n_available_interfaces']} interfaces)"
    )


def format_binder_success(value: dict[str, Any]) -> str:
    if not value["n_available"]:
        return "NA"
    return (
        f"{format_number(value['success_rate_pct'])}% "
        f"({value['n_available']} predictions; "
        f"{value['n_success_groups']}/{value['n_available_groups']} groups)"
    )


def markdown(report: dict[str, Any]) -> str:
    folding_thresholds = report["folding_success_thresholds"]
    target_folding_threshold = folding_thresholds[
        "target_self_aligned_backbone_rmsd_angstrom"
    ]["value"]
    binder_folding_threshold = folding_thresholds[
        "binder_self_aligned_backbone_rmsd_angstrom"
    ]["value"]
    thresholds = report["docking_success_thresholds"]
    rmsd_threshold = thresholds["target_aligned_binder_backbone_rmsd_angstrom"]["value"]
    dockq_threshold = thresholds["target_binder_interface_dockq"]["value"]
    lines = [
        "# Batch metric summary",
        "",
        "lDDT and iLDDT columns show OpenStructure scores with stereochemistry checks "
        "on and off. NA indicates that a score is unavailable.",
        "",
        "## Folding",
        "",
        f"Strict success criteria: target self-aligned backbone RMSD < "
        f"{target_folding_threshold:g} Å; binder self-aligned backbone RMSD < "
        f"{binder_folding_threshold:g} Å. Unavailable values are excluded from "
        "each denominator. Chirality accuracy pools passed and assessable residue counts "
        "across complex–seed predictions.",
        "",
        "| Method | Predictions | Target self-aligned backbone RMSD, median (Å) | "
        "Target RMSD success | Binder self-aligned backbone RMSD, median of prediction means (Å) | "
        "Binder RMSD success | Mean target-only lDDT (checks on) | Mean target-only lDDT (checks off) | "
        "Mean binder-only lDDT (checks on) | Mean binder-only lDDT (checks off) | "
        "All (D+L) chirality accuracy | D-only chirality accuracy | L-only chirality accuracy | "
        "Unassessable chirality sites |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for item in report["methods"]:
        folding = item["folding"]
        chirality = folding["chirality_accuracy"]
        lines.append(
            "| {method} | {n} | {target} | {target_success} | {binder} | "
            "{binder_success} | {target_lddt} | {target_lddt_off} | "
            "{binder_lddt} | {binder_lddt_off} | {all_chirality} | {daa_chirality} | "
            "{laa_chirality} | {unassessable} |".format(
                method=item["method"],
                n=item["n_predictions"],
                target=format_number(
                    folding["target_self_aligned_backbone_rmsd_angstrom_median"]["value"]
                ),
                target_success=format_success(
                    folding["target_self_aligned_backbone_rmsd_below_threshold"]
                ),
                binder=format_number(
                    folding["binder_self_aligned_backbone_rmsd_angstrom_median"]["value"]
                ),
                binder_success=format_binder_success(
                    folding["binder_self_aligned_backbone_rmsd_below_threshold"]
                ),
                target_lddt=format_number(folding["target_only_lddt_mean"]["value"]),
                target_lddt_off=format_number(folding["target_only_lddt_mean_no_stereochecks"]["value"]),
                binder_lddt=format_number(folding["binder_only_lddt_mean"]["value"]),
                binder_lddt_off=format_number(folding["binder_only_lddt_mean_no_stereochecks"]["value"]),
                all_chirality=format_chirality(chirality["all"]),
                daa_chirality=format_chirality(chirality["daa"]),
                laa_chirality=format_chirality(chirality["laa"]),
                unassessable=chirality["all"]["n_unassessable"],
            )
        )

    lines.extend(
        [
            "",
            "## Docking",
            "",
            f"Strict success criteria: DockQ > {dockq_threshold:g}; "
            f"L-RMSD < {rmsd_threshold:g} Å. Unavailable values are excluded from "
            "each denominator.",
            "",
            "L-RMSD is binder backbone (N/CA/C/O) RMSD after jointly fitting all "
            "annotated target chains on their backbone atoms, without refitting the binder.",
            "",
            "| Method | Predictions | Mean target–binder DockQ | Interface DockQ success | "
            "L-RMSD, median of prediction means (Å) | L-RMSD success | "
            "Target–binder iRMSD, median of prediction means (Å) | "
            "Mean target–binder iLDDT (macro, checks on) | Mean target–binder iLDDT (macro, checks off) |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for item in report["methods"]:
        docking = item["docking"]
        lines.append(
            "| {method} | {n} | {dockq} | {dockq_success} | {binder} | "
            "{binder_success} | {irmsd} | {ilddt} | {ilddt_off} |".format(
                method=item["method"],
                n=item["n_predictions"],
                dockq=format_number(docking["dockq_mean"]["value"]),
                dockq_success=format_dockq_success(docking["dockq_above_threshold"]),
                binder=format_number(
                    docking["target_aligned_binder_backbone_rmsd_angstrom_median"]["value"]
                ),
                binder_success=format_binder_success(
                    docking[
                        "target_aligned_binder_backbone_rmsd_below_threshold"
                    ]
                ),
                irmsd=format_number(docking["irmsd_angstrom_median"]["value"]),
                ilddt=format_number(docking["ilddt_macro_mean"]["value"]),
                ilddt_off=format_number(docking["ilddt_macro_mean_no_stereochecks"]["value"]),
            )
        )

    lines.extend(
        [
            "",
            "RMSDs are medians across predictions; binder and interface RMSDs "
            "first average available groups or interfaces within each prediction. "
            "lDDT, DockQ, and iLDDT are arithmetic means across predictions. DockQ success "
            "first applies its cutoff to each target–binder interface, then "
            "averages interface success within each prediction and across predictions. "
            "Binder self-aligned RMSD success and L-RMSD success apply their strict "
            "cutoffs to each available group, "
            "then average group pass fractions within and across predictions. "
            "The parenthetical pooled group counts need not divide to the "
            "reported prediction-balanced percentage. "
            "All threshold classifications are independent. "
            "Skipped predictions are excluded.",
            "",
            "## Supplemental metrics",
            "",
            "| Method | Mean global scored-subset lDDT (checks on) | "
            "Mean global scored-subset lDDT (checks off) | Mean target–binder DockQ-wave | "
            "Mean target–binder iLDDT (micro, checks on) | Mean target–binder iLDDT (micro, checks off) |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for item in report["methods"]:
        folding = item["folding"]
        docking = item["docking"]
        lines.append(
            "| {method} | {lddt} | {lddt_off} | {wave} | {micro} | {micro_off} |".format(
                method=item["method"],
                lddt=format_number(folding["global_scored_subset_lddt_mean"]["value"]),
                lddt_off=format_number(folding["global_scored_subset_lddt_mean_no_stereochecks"]["value"]),
                wave=format_number(docking["dockq_wave_mean"]["value"]),
                micro=format_number(docking["ilddt_micro_mean"]["value"]),
                micro_off=format_number(docking["ilddt_micro_mean_no_stereochecks"]["value"]),
            )
        )
    lines.append("")
    return "\n".join(lines)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def main() -> int:
    args = parse_args()
    with args.batch_csv.open(encoding="utf-8", newline="") as handle:
        report = summarize(
            list(csv.DictReader(handle)),
            args.batch_csv,
            args.binder_rmsd_threshold_angstrom,
            args.dockq_threshold,
            args.target_self_aligned_rmsd_threshold_angstrom,
            args.binder_self_aligned_rmsd_threshold_angstrom,
            args.method_name,
        )
    write_text(args.json_output, json.dumps(report, indent=2, allow_nan=False) + "\n")
    write_text(args.markdown_output, markdown(report))
    print(args.markdown_output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
