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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch_csv", type=Path)
    parser.add_argument("--json-output", type=Path, required=True)
    parser.add_argument("--markdown-output", type=Path, required=True)
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


def load_interfaces(raw: str) -> list[dict[str, Any]]:
    if not raw:
        return []
    try:
        value = json.loads(raw)
    except json.JSONDecodeError:
        return []
    return value if isinstance(value, list) else []


def summarize(rows: list[dict[str, str]], source: Path) -> dict[str, Any]:
    eligible = [row for row in rows if row.get("status") != "SKIPPED"]
    by_method: dict[str, list[dict[str, str]]] = {}
    for row in eligible:
        by_method.setdefault(method_from_prediction(row.get("prediction_path", "")), []).append(row)

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
        lddt = [
            value
            for row in prediction_rows
            if (value := finite_float(row.get("ost_lddt"))) is not None
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
            if (value := finite_float(row.get("ost_dockq_ave"))) is not None
        ]
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
        ilddt = [
            value
            for row in prediction_rows
            if (value := finite_float(row.get("ost_ilddt"))) is not None
        ]
        irmsd: list[float] = []
        for row in prediction_rows:
            for interface in load_interfaces(row.get("ost_per_interface_json", "")):
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
                    "binder_self_aligned_backbone_rmsd_angstrom_median": aggregate(
                        binder_fold, statistics.median
                    ),
                    "lddt_mean": aggregate(lddt, statistics.mean),
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
                    "target_aligned_binder_backbone_rmsd_angstrom_median": aggregate(
                        target_aligned_binder, statistics.median
                    ),
                    "irmsd_angstrom_median": aggregate(irmsd, statistics.median),
                    "ilddt_mean": aggregate(ilddt, statistics.mean),
                },
            }
        )

    return {
        "schema_version": 1,
        "source_batch_csv": str(source.resolve()),
        "aggregation": {
            "prediction_level_scores": "Deduplicated by prediction path before aggregation",
            "rmsd": "Median across available prediction, binder-group, or interface values",
            "bounded_scores": "Arithmetic mean across available prediction values",
            "chirality_violation": (
                "Residue-pooled (inverted + below-threshold) / assessable residues; "
                "unassessable residues reported separately"
            ),
        },
        "methods": methods,
    }


def format_number(value: Any, digits: int = 3) -> str:
    return "NA" if value is None else f"{value:.{digits}f}"


def markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Batch metric summary",
        "",
        "## Folding",
        "",
        "| Method | Predictions | Target self-aligned backbone RMSD, median (Å) | "
        "Binder self-aligned backbone RMSD, median (Å) | Mean lDDT | "
        "Chirality violations |",
        "|---|---:|---:|---:|---:|---:|",
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
            "| {method} | {n} | {target} | {binder} | {lddt} | {violation} |".format(
                method=item["method"],
                n=item["n_predictions"],
                target=format_number(
                    folding["target_self_aligned_backbone_rmsd_angstrom_median"]["value"]
                ),
                binder=format_number(
                    folding["binder_self_aligned_backbone_rmsd_angstrom_median"]["value"]
                ),
                lddt=format_number(folding["lddt_mean"]["value"]),
                violation=violation_text,
            )
        )

    lines.extend(
        [
            "",
            "## Docking",
            "",
            "| Method | Predictions | Mean DockQ | Target-aligned binder backbone RMSD, "
            "median (Å) | iRMSD, median (Å) | Mean iLDDT |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for item in report["methods"]:
        docking = item["docking"]
        lines.append(
            "| {method} | {n} | {dockq} | {binder} | {irmsd} | {ilddt} |".format(
                method=item["method"],
                n=item["n_predictions"],
                dockq=format_number(docking["dockq_mean"]["value"]),
                binder=format_number(
                    docking["target_aligned_binder_backbone_rmsd_angstrom_median"]["value"]
                ),
                irmsd=format_number(docking["irmsd_angstrom_median"]["value"]),
                ilddt=format_number(docking["ilddt_mean"]["value"]),
            )
        )

    lines.extend(
        [
            "",
            "RMSDs are medians; lDDT, DockQ, and iLDDT are arithmetic means. "
            "Chirality violation is residue-pooled over assessable residues, with "
            "unassessable residues shown separately. Skipped predictions are excluded.",
            "",
        ]
    )
    return "\n".join(lines)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def main() -> int:
    args = parse_args()
    with args.batch_csv.open(encoding="utf-8", newline="") as handle:
        report = summarize(list(csv.DictReader(handle)), args.batch_csv)
    write_text(args.json_output, json.dumps(report, indent=2, allow_nan=False) + "\n")
    write_text(args.markdown_output, markdown(report))
    print(args.markdown_output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
