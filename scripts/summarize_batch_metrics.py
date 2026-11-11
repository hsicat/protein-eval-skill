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
            "Strict lower cutoff; otherwise infer it from the CSV or use "
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


def method_from_prediction(path: str) -> str:
    parts = PurePosixPath(path).parts
    if len(parts) >= 2 and parts[0] == "predicted":
        return parts[1]
    return parts[0] if parts else "unknown"


def aggregate(values: list[float], operation: Callable[[list[float]], float]) -> dict[str, Any]:
    return {"value": operation(values) if values else None, "n": len(values)}


def success_rate(values: list[float], predicate: Callable[[float], bool]) -> dict[str, Any]:
    n_success = sum(predicate(value) for value in values)
    return {
        "n_success": n_success,
        "n_available": len(values),
        "success_rate_pct": 100.0 * n_success / len(values) if values else None,
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
        binder_fold = [
            value
            for row in method_rows
            if (value := finite_float(row.get("binder_self_aligned_backbone_rmsd_angstrom")))
            is not None
        ]
        global_lddt = [
            value
            for row in prediction_rows
            if (value := finite_float(row.get("ost_lddt"))) is not None
        ]
        target_lddt = [
            value for row in prediction_rows
            if (value := finite_float(row.get("target_only_lddt"))) is not None
        ]
        binder_lddt = [
            value for row in prediction_rows
            if (value := finite_float(row.get("binder_only_lddt_complex_mean"))) is not None
        ]

        chirality_sites = sum(integer(row.get("chirality_all_eligible_count")) for row in prediction_rows)
        chirality_pass = sum(integer(row.get("chirality_all_correct_count")) for row in prediction_rows)
        chirality_unassessable = sum(
            integer(row.get("chirality_all_skipped_count")) for row in prediction_rows
        )
        chirality_assessable = chirality_sites - chirality_unassessable
        chirality_violations = max(0, chirality_assessable - chirality_pass)

        dockq = [
            value
            for row in prediction_rows
            if (value := finite_float(row.get("target_binder_dockq_mean"))) is not None
        ]
        dockq_wave = [
            value for row in prediction_rows
            if (value := finite_float(row.get("target_binder_dockq_wave"))) is not None
        ]
        dockq_success_fractions = [
            value for row in prediction_rows
            if (value := finite_float(row.get("target_binder_dockq_success_fraction"))) is not None
        ]
        dockq_interface_successes = sum(
            integer(row.get("target_binder_dockq_interface_success_count"))
            for row in prediction_rows
        )
        dockq_interface_available = sum(
            integer(row.get("target_binder_dockq_interface_available_count"))
            for row in prediction_rows
        )
        target_aligned_binder = [
            value
            for row in method_rows
            if (
                value := finite_float(
                    row.get("target_aligned_binder_backbone_rmsd_angstrom")
                )
            )
            is not None
        ]
        ilddt_macro = [
            value
            for row in prediction_rows
            if (value := finite_float(row.get("target_binder_ilddt_macro"))) is not None
        ]
        ilddt_micro = [
            value for row in prediction_rows
            if (value := finite_float(row.get("target_binder_ilddt_micro"))) is not None
        ]
        irmsd: list[float] = []
        for row in prediction_rows:
            for interface in load_interfaces(row.get("target_binder_dockq_per_interface_json", "")):
                value = finite_float(interface.get("irmsd"))
                if value is not None:
                    irmsd.append(value)

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
                    "binder_self_aligned_backbone_rmsd_below_threshold": {
                        "threshold_angstrom": (
                            binder_self_aligned_rmsd_threshold_angstrom
                        ),
                        "operator": "<",
                        **success_rate(
                            binder_fold,
                            lambda value: (
                                value < binder_self_aligned_rmsd_threshold_angstrom
                            ),
                        ),
                    },
                    "target_only_lddt_mean": aggregate(target_lddt, statistics.mean),
                    "binder_only_lddt_mean": aggregate(binder_lddt, statistics.mean),
                    "global_scored_subset_lddt_mean": aggregate(global_lddt, statistics.mean),
                    "chirality_violation": {
                        "n_violations": chirality_violations,
                        "n_assessable": chirality_assessable,
                        "n_unassessable": chirality_unassessable,
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
                    "target_aligned_binder_backbone_rmsd_below_threshold": {
                        "threshold_angstrom": binder_rmsd_threshold_angstrom,
                        "operator": "<",
                        **success_rate(
                            target_aligned_binder,
                            lambda value: value < binder_rmsd_threshold_angstrom,
                        ),
                    },
                    "irmsd_angstrom_median": aggregate(irmsd, statistics.median),
                    "ilddt_macro_mean": aggregate(ilddt_macro, statistics.mean),
                    "ilddt_micro_mean": aggregate(ilddt_micro, statistics.mean),
                },
            }
        )

    return {
        "schema_version": 4,
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
            "rmsd": "Median across available prediction, binder-group, or interface values",
            "bounded_scores": "Arithmetic mean across available prediction values",
            "folding_success": (
                "Independent strict-threshold success proportions; target RMSD is "
                "counted per prediction, binder RMSD per binder group, and unavailable "
                "values are excluded from each denominator"
            ),
            "docking_success": (
                "DockQ threshold is applied per target-binder interface, then the "
                "success fraction is averaged within each complex and across "
                "available predictions; RMSD thresholds remain independent"
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
        },
        "methods": methods,
    }


def format_number(value: Any, digits: int = 3) -> str:
    return "NA" if value is None else f"{value:.{digits}f}"


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
        f"({value['n_available']} complexes; "
        f"{value['n_success_interfaces']}/{value['n_available_interfaces']} interfaces)"
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
        "## Folding",
        "",
        f"Strict success criteria: target self-aligned backbone RMSD < "
        f"{target_folding_threshold:g} Å; binder self-aligned backbone RMSD < "
        f"{binder_folding_threshold:g} Å. Unavailable values are excluded from "
        "each denominator.",
        "",
        "| Method | Predictions | Target self-aligned backbone RMSD, median (Å) | "
        "Target RMSD success | Binder self-aligned backbone RMSD, median (Å) | "
        "Binder RMSD success | Mean target-only lDDT | Mean binder-only lDDT | Chirality violations |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for item in report["methods"]:
        folding = item["folding"]
        violation = folding["chirality_violation"]
        violation_text = (
            f"{violation['n_violations']}/{violation['n_assessable']} "
            f"({format_number(violation['violation_rate_pct'])}%)"
        )
        if violation["n_unassessable"]:
            violation_text += f"; {violation['n_unassessable']} unassessable"
        lines.append(
            "| {method} | {n} | {target} | {target_success} | {binder} | "
            "{binder_success} | {target_lddt} | {binder_lddt} | {violation} |".format(
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
                binder_success=format_success(
                    folding["binder_self_aligned_backbone_rmsd_below_threshold"]
                ),
                target_lddt=format_number(folding["target_only_lddt_mean"]["value"]),
                binder_lddt=format_number(folding["binder_only_lddt_mean"]["value"]),
                violation=violation_text,
            )
        )

    lines.extend(
        [
            "",
            "## Docking",
            "",
            f"Strict success criteria: DockQ > {dockq_threshold:g}; target-aligned binder "
            f"backbone RMSD < {rmsd_threshold:g} Å. Unavailable values are excluded from "
            "each denominator.",
            "",
            "| Method | Predictions | Mean target–binder DockQ | Interface DockQ success | Target-aligned binder "
            "backbone RMSD, median (Å) | RMSD success | Target–binder iRMSD, median (Å) | Mean target–binder iLDDT (macro) |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for item in report["methods"]:
        docking = item["docking"]
        lines.append(
            "| {method} | {n} | {dockq} | {dockq_success} | {binder} | "
            "{binder_success} | {irmsd} | {ilddt} |".format(
                method=item["method"],
                n=item["n_predictions"],
                dockq=format_number(docking["dockq_mean"]["value"]),
                dockq_success=format_dockq_success(docking["dockq_above_threshold"]),
                binder=format_number(
                    docking["target_aligned_binder_backbone_rmsd_angstrom_median"]["value"]
                ),
                binder_success=format_success(
                    docking[
                        "target_aligned_binder_backbone_rmsd_below_threshold"
                    ]
                ),
                irmsd=format_number(docking["irmsd_angstrom_median"]["value"]),
                ilddt=format_number(docking["ilddt_macro_mean"]["value"]),
            )
        )

    lines.extend(
        [
            "",
            "RMSDs are medians; lDDT, DockQ, and iLDDT are arithmetic means across "
            "the same prediction or binder-group units as before. DockQ success "
            "first applies its cutoff to each target–binder interface, then "
            "averages interface success within each complex and across predictions. "
            "Target and binder self-aligned backbone RMSD successes are evaluated "
            "independently with strict inequalities. "
            "DockQ and target-aligned binder backbone RMSD successes are evaluated "
            "independently with strict inequalities. "
            "Chirality violation is residue-pooled over assessable residues, with "
            "unassessable residues shown separately. Skipped predictions are excluded.",
            "",
            "## Supplemental metrics",
            "",
            "| Method | Mean global scored-subset lDDT | Mean target–binder DockQ-wave | Mean target–binder iLDDT (micro) |",
            "|---|---:|---:|---:|",
        ]
    )
    for item in report["methods"]:
        folding = item["folding"]
        docking = item["docking"]
        lines.append(
            "| {method} | {lddt} | {wave} | {micro} |".format(
                method=item["method"],
                lddt=format_number(folding["global_scored_subset_lddt_mean"]["value"]),
                wave=format_number(docking["dockq_wave_mean"]["value"]),
                micro=format_number(docking["ilddt_micro_mean"]["value"]),
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
