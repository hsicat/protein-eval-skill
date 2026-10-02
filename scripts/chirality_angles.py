#!/usr/bin/env python3
"""Score CA chirality using target-oriented improper angles, never volume signs.

Requires numpy and gemmi; no model weights, GPU, JAX, or RDKit at runtime.
See README.md for target assumptions, atom order, and metric definitions.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

DAA_TO_LAA = dict(zip(
    "DAL DAR DAS DCY DGL DGN DHI DIL DLE DLY DPN DPR DSG DSN DTH DTR DTY DVA MED".split(),
    "ALA ARG ASP CYS GLU GLN HIS ILE LEU LYS PHE PRO ASN SER THR TRP TYR VAL MET".split(),
))
LAA = frozenset(DAA_TO_LAA.values())
FIELDS = [
    "model", "chain", "residue_id", "label_seq_id", "residue_name", "target_ccd",
    "target_dl", "center", "target_rs", "atom_order", "altloc", "status",
    "raw_angle_rad", "raw_angle_deg", "oriented_angle_rad", "oriented_angle_deg",
    "threshold_deg", "margin_deg", "deficit_deg", "loss_rad", "passed", "error",
]


def target_spec(code):
    """CIP descending heavy-neighbour order followed by CA, for 38 CCDs."""
    is_d = code in DAA_TO_LAA
    parent = DAA_TO_LAA.get(code, code)
    if parent not in LAA:
        raise ValueError(f"Unsupported target CCD {code!r}; supported: 19 DAA + 19 LAA.")
    # Cysteine's sulfur reverses C versus CB CIP priority and the CA R/S label.
    cysteine = parent == "CYS"
    order = ("N", "CB", "C", "CA") if cysteine else ("N", "C", "CB", "CA")
    target_r = is_d != cysteine
    return order, target_r, "D" if is_d else "L"


def improper_angle(positions):
    """Same atan2 convention as AF3 daa_bp.signed_improper_torsion (float64)."""
    p = np.asarray(positions, dtype=np.float64)
    if p.shape != (4, 3) or not np.isfinite(p).all():
        raise ValueError("Four finite XYZ coordinates are required.")
    a, b, c = p[0] - p[1], p[2] - p[1], p[2] - p[3]
    n1, n2 = np.cross(a, b), np.cross(b, c)
    if min(np.linalg.norm(b), np.linalg.norm(n1), np.linalg.norm(n2)) <= 1e-12:
        raise ValueError("Degenerate improper geometry (zero axis or collinear plane atoms).")
    return float(np.arctan2(np.dot(b / np.linalg.norm(b), np.cross(n1, n2)),
                           np.dot(n1, n2)))


def angle_metrics(phi, target_r, threshold_deg):
    oriented = phi if target_r else -phi
    threshold = math.radians(threshold_deg)
    loss = max(0.0, threshold - oriented)
    # An exactly coplanar configuration can occur at either end of atan2.
    # Flag +/-pi as well as 0 instead of treating +180 degrees as excellent.
    planar = abs(math.sin(phi)) < 1e-10
    if planar:
        raise ValueError("Coplanar improper geometry; handedness is unassessable.")
    return dict(raw_angle_rad=phi, raw_angle_deg=math.degrees(phi),
                oriented_angle_rad=oriented, oriented_angle_deg=math.degrees(oriented),
                margin_deg=math.degrees(oriented) - threshold_deg,
                deficit_deg=math.degrees(loss), loss_rad=loss, passed=loss == 0.0,
                status="pass" if loss == 0.0 else ("inverted" if oriented < 0 else "below_threshold"))


def residue_positions(residue, order):
    """Choose one coherent alternate conformer, shared blank atoms allowed."""
    atoms = [a for a in residue if a.name.strip() in order and a.occ > 0]
    def label(atom):
        return atom.altloc.strip("\x00 ")
    labels = sorted({label(a) for a in atoms if label(a)}) or [""]
    candidates = []
    for alt in labels:
        selected = []
        for name in order:
            matches = [a for a in atoms if a.name.strip() == name and label(a) in ("", alt)]
            if not matches:
                break
            selected.append(max(matches, key=lambda a: (a.occ, label(a) == alt)))
        if len(selected) == 4:
            candidates.append((sum(a.occ for a in selected), alt, selected))
    if not candidates:
        raise ValueError("Missing/zero-occupancy atoms or no complete coherent altloc: " + ",".join(order))
    _, alt, selected = sorted(candidates, key=lambda x: (-x[0], x[1]))[0]
    return [[a.pos.x, a.pos.y, a.pos.z] for a in selected], alt


def summarize(rows):
    valid = [r for r in rows if r.get("passed") is not None]
    passed = sum(r["passed"] for r in valid)
    n, nv = len(rows), len(valid)
    loss_sum = sum(r["loss_rad"] for r in valid)
    complete = n > 0 and n == nv
    return dict(
        n_sites=n, n_assessable=nv, n_unassessable=n - nv, n_pass=passed,
        n_inverted=sum(r["status"] == "inverted" for r in rows),
        n_below_threshold=sum(r["status"] == "below_threshold" for r in rows),
        # This conservative rate retains invalid sites in the denominator.
        chirality_correct_rate_pct=100 * passed / n if n else None,
        assessable_correct_rate_pct=100 * passed / nv if nv else None,
        coverage_pct=100 * nv / n if n else None,
        total_chirality_loss_rad=loss_sum if complete else None,
        mean_chirality_loss_rad=loss_sum / n if complete else None,
        assessed_loss_sum_rad=loss_sum if nv else None,
        assessed_loss_mean_rad=loss_sum / nv if nv else None,
        complete=complete,
    )


def evaluate(path, threshold_deg=15.0, include_laa=False, chains=None, sites=None):
    import gemmi
    structure = gemmi.read_structure(str(path))
    if not len(structure):
        raise ValueError("Structure contains no models.")
    rows, excluded = [], []
    for model_index, model in enumerate(structure, start=1):
        model_rows, found = [], set()
        present_chains = {chain.name for chain in model}
        if chains and not set(chains).issubset(present_chains):
            raise ValueError(f"Model {model_index}: unknown chains {set(chains) - present_chains}")
        for chain in model:
            if chains and chain.name not in chains:
                continue
            for residue in chain:
                code = residue.name.strip().upper()
                key = (chain.name, str(residue.seqid))
                if sites is not None:
                    if key not in sites:
                        continue
                    if key in found:
                        raise ValueError(f"Ambiguous repeated residue identity {key} in model {model_index}.")
                    target = sites[key]
                else:
                    target = code
                    if code not in DAA_TO_LAA and not (include_laa and code in LAA):
                        if code not in LAA and code not in {"GLY", "HOH", "WAT", "DOD"}:
                            if any(a.name.strip() == "CA" for a in residue):
                                excluded.append(dict(model=model_index, chain=chain.name,
                                                     residue_id=str(residue.seqid), residue_name=code,
                                                     reason="unsupported residue with CA; not scored"))
                        continue
                    if key in found:
                        raise ValueError(f"Ambiguous repeated residue identity {key} in model {model_index}.")
                found.add(key)
                order, target_r, dl = target_spec(target)
                row = dict(model=model_index, chain=chain.name, residue_id=str(residue.seqid),
                           label_seq_id=residue.label_seq, residue_name=code, target_ccd=target,
                           target_dl=dl, center="CA", target_rs="R" if target_r else "S",
                           atom_order=",".join(order), threshold_deg=threshold_deg,
                           passed=None, status="unassessable", error="")
                try:
                    positions, alt = residue_positions(residue, order)
                    row["altloc"] = alt
                    phi = improper_angle(positions)
                    row.update(raw_angle_rad=phi, raw_angle_deg=math.degrees(phi))
                    row.update(angle_metrics(phi, target_r, threshold_deg))
                except ValueError as exc:
                    row["error"] = str(exc)
                model_rows.append(row)
        if sites is not None:
            for key in sorted(set(sites) - found):
                order, target_r, dl = target_spec(sites[key])
                model_rows.append(dict(model=model_index, chain=key[0], residue_id=key[1],
                                       target_ccd=sites[key], target_dl=dl, center="CA",
                                       target_rs="R" if target_r else "S", atom_order=",".join(order),
                                       threshold_deg=threshold_deg, passed=None, status="unassessable",
                                       error="Requested residue is missing from structure."))
        rows.extend(model_rows)
    def grouped(selected):
        return {"DAA": summarize([r for r in selected if r["target_dl"] == "D"]),
                "LAA": summarize([r for r in selected if r["target_dl"] == "L"]),
                "all_selected": summarize(selected)}
    report = dict(
        schema_version=1, input_file=str(Path(path).resolve()),
        input_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        algorithm="CIP-ordered signed improper atan2; target R:+phi, S:-phi",
        center="CA only", threshold_deg=threshold_deg,
        target_source="explicit --site" if sites is not None else "structure CCD residue names",
        selection=dict(include_laa=include_laa, chains=chains),
        accuracy_definition="pass / all selected sites; pass requires oriented angle >= threshold",
        loss_definition="sum or mean of max(0, radians(threshold) - oriented_angle_rad)",
        invalid_policy="invalid sites retained in accuracy denominator; full loss null if any invalid",
        pooled_note="Each model-residue pair counts once; models also summarized separately.",
        pooled=grouped(rows),
        models=[dict(model=i, **grouped([r for r in rows if r["model"] == i]))
                for i in range(1, len(structure) + 1)],
        excluded_unsupported_residues=excluded,
    )
    return rows, report


def simple_output(rows, report):
    """English presentation only; all decisions use the original full precision."""
    multiple = len(report["models"]) > 1
    headers = (["Model"] if multiple else []) + [
        "Site", "Residue", "Target", "Angle (deg)", "Gap (deg)", "Result"]
    labels = {"pass": "Pass", "inverted": "Fail: wrong handedness",
              "below_threshold": "Fail: below threshold", "unassessable": "Cannot assess"}
    table = []
    for row in rows:
        number = lambda key: "N/A" if row.get(key) is None else f"{row[key]:.2f}"
        values = ([str(row["model"])] if multiple else []) + [
            f"{row['chain']}:{row['residue_id']}", row.get("residue_name") or "Missing",
            f"{row['target_dl']} ({row['target_ccd']})", number("oriented_angle_deg"),
            number("deficit_deg"), labels[row["status"]]]
        table.append(values)
    widths = [max(len(h), *(len(r[i]) for r in table)) if table else len(h)
              for i, h in enumerate(headers)]
    line = lambda values: "  ".join(v.ljust(w) for v, w in zip(values, widths)).rstrip()
    lines = ["Chirality report", f"File: {Path(report['input_file']).name}",
             f"Pass threshold: {report['threshold_deg']:g} deg (CA centers)", "",
             line(headers), line(["-" * w for w in widths])]
    lines.extend(line(r) for r in table)
    lines += ["", "Angle is oriented toward the target D/L chirality.",
              "Gap = max(0, threshold - angle). A zero gap passes.",
              "Values are rounded; pass/fail uses full precision.", ""]
    def summary_lines(groups):
        result = []
        for group in ("DAA", "LAA"):
            s = groups[group]
            if s["n_sites"]:
                result.append(f"{group} chirality accuracy: {s['chirality_correct_rate_pct']:.1f}% "
                              f"({s['n_pass']}/{s['n_sites']} pass)")
        s = groups["all_selected"]
        loss = s["total_chirality_loss_rad"]
        shown = "N/A" if loss is None else f"{loss:.5f} rad"
        result.append(f"Total chirality loss (sum): {shown}")
        if s["n_unassessable"]:
            result.append(f"Cannot assess: {s['n_unassessable']} site(s), counted as not passing; "
                          "total loss is unavailable.")
        return result
    if multiple:
        lines.append("Combined results (each model-site pair counts once):")
    lines.extend(summary_lines(report["pooled"]))
    if multiple:
        for model in report["models"]:
            lines += ["", f"Model {model['model']}:"] + summary_lines(model)
    if not rows:
        lines += ["No target sites found. If DAA names are absent, specify targets with --site."]
    for row in rows:
        if row.get("error"):
            lines.append(f"Note: model {row['model']}, {row['chain']}:{row['residue_id']}: {row['error']}")
    for row in report["excluded_unsupported_residues"]:
        lines.append(f"Excluded: model {row['model']}, {row['chain']}:{row['residue_id']} "
                     f"({row['residue_name']}) — unsupported residue.")
    return headers, table, "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structure", type=Path, help="PDB/mmCIF structure (all models are evaluated).")
    parser.add_argument("--out-dir", type=Path, help="New output directory; existing directory is refused.")
    parser.add_argument("--threshold-deg", type=float, default=15.0)
    parser.add_argument("--include-laa", action="store_true", help="Also score the 19 standard L amino acids.")
    parser.add_argument("--detailed", action="store_true",
                        help="Also write full diagnostics to detailed_per_residue.csv and summary.json.")
    parser.add_argument("--chain", action="append", help="Restrict to author chain ID; repeatable.")
    parser.add_argument("--site", action="append", metavar="CHAIN:RESID:CCD",
                        help="Explicit target site (author IDs, e.g. B:20:DAL or B:20:ALA); repeatable.")
    args = parser.parse_args(argv)
    try:
        if not math.isfinite(args.threshold_deg) or not 0 < args.threshold_deg < 90:
            raise ValueError("--threshold-deg must be finite and between 0 and 90.")
        if args.site and (args.chain or args.include_laa):
            raise ValueError("--site cannot be combined with --chain or --include-laa.")
        sites = None
        if args.site:
            sites = {}
            for value in args.site:
                parts = value.split(":")
                if len(parts) != 3 or not parts[1]:
                    raise ValueError(f"Invalid --site {value!r}; expected CHAIN:RESID:CCD.")
                chain, resid, code = parts
                code = code.upper()
                target_spec(code)
                if (chain, resid) in sites:
                    raise ValueError(f"Duplicate --site {chain}:{resid}.")
                sites[chain, resid] = code
        out = args.out_dir or Path.cwd() / (args.structure.stem + "_chirality_angles")
        if out.exists():
            raise ValueError(f"Output directory already exists: {out}; choose a new --out-dir.")
        rows, report = evaluate(args.structure, args.threshold_deg, args.include_laa, args.chain, sites)
        out.mkdir(parents=True, exist_ok=False)
        headers, table, readable = simple_output(rows, report)
        with (out / "per_residue.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(headers)
            writer.writerows(table)
        (out / "report.txt").write_text(readable, encoding="utf-8")
        if args.detailed:
            with (out / "detailed_per_residue.csv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS)
                writer.writeheader()
                writer.writerows(rows)
            (out / "summary.json").write_text(json.dumps(report, indent=2, ensure_ascii=False,
                                                          allow_nan=False) + "\n", encoding="utf-8")
        print(readable, end="")
        print(f"Results: {out.resolve()}")
        if not rows:
            return 2
        return 2 if any(r["passed"] is None for r in rows) else 0
    except (ValueError, OSError, RuntimeError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
