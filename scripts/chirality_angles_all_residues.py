#!/usr/bin/env python3
"""Score the C-alpha chirality of all standard chiral amino acids in one chain.

This companion to chirality_angles.py scores both the 19 standard L amino
acids and their 19 D-amino-acid CCD counterparts.  Glycine is omitted because
its C-alpha has two hydrogen atoms and is therefore achiral.

Requires numpy and gemmi.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import chirality_angles as chirality


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structure", type=Path,
                        help="PDB/mmCIF structure; all models are evaluated.")
    parser.add_argument("--chain", required=True, action="append",
                        help="Author chain ID to score; repeat for additional chains.")
    parser.add_argument("--out-dir", type=Path,
                        help="Output directory; all-residue files are added without replacing other reports.")
    parser.add_argument("--threshold-deg", type=float, default=15.0,
                        help="Minimum target-oriented improper angle for a pass (default: 15).")
    parser.add_argument("--detailed", action="store_true",
                        help="Also write all_residue_chirality_detailed.csv and summary.json.")
    args = parser.parse_args(argv)

    try:
        if not math.isfinite(args.threshold_deg) or not 0 < args.threshold_deg < 90:
            raise ValueError("--threshold-deg must be finite and between 0 and 90.")

        out = args.out_dir or Path.cwd() / (args.structure.stem + "_all_residue_chirality")
        outputs = [out / "all_residue_chirality.csv", out / "all_residue_chirality_report.txt"]
        if args.detailed:
            outputs += [out / "all_residue_chirality_detailed.csv",
                        out / "all_residue_chirality_summary.json"]
        existing = [path.name for path in outputs if path.exists()]
        if existing:
            raise ValueError("Refusing to overwrite existing all-residue output(s): " + ", ".join(existing))

        # include_laa=True changes the original script's DAA-only default to
        # score all supported chiral standard amino-acid residue codes.
        rows, report = chirality.evaluate(
            args.structure,
            threshold_deg=args.threshold_deg,
            include_laa=True,
            chains=args.chain,
        )
        report["target_source"] = "residue CCD names in requested chain(s); DAA and LAA enabled"

        out.mkdir(parents=True, exist_ok=True)
        headers, table, readable = chirality.simple_output(rows, report)
        with (out / "all_residue_chirality.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(headers)
            writer.writerows(table)
        (out / "all_residue_chirality_report.txt").write_text(readable, encoding="utf-8")
        if args.detailed:
            with (out / "all_residue_chirality_detailed.csv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=chirality.FIELDS)
                writer.writeheader()
                writer.writerows(rows)
            (out / "all_residue_chirality_summary.json").write_text(
                json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                encoding="utf-8",
            )
        print(readable, end="")
        print(f"Results: {out.resolve()}")
        if not rows:
            return 2
        return 2 if any(row["passed"] is None for row in rows) else 0
    except (ValueError, OSError, RuntimeError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
