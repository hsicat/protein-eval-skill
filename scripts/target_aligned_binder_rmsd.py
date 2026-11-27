#!/usr/bin/env python3
"""Compute target-aligned and self-aligned binder RMSDs from an OST alignment.

The model target is fitted once to the reference target with matched N, CA, C,
and O atoms.  That transform is then applied without refitting to every binder.
Each binder is reported using both CA-only and N/CA/C/O atom selections.

Each binder group is also fitted to itself, independently for CA and backbone
atoms, to measure binder folding/internal-assembly accuracy without target pose.

Requires numpy and gemmi (the repository's ``urop`` conda environment).
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import gemmi
import numpy as np


BACKBONE = ("N", "CA", "C", "O")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--ost-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def load_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def chain_residues(path: Path, chain_id_type: str) -> dict[str, list[Any]]:
    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        raise ValueError(f"Expected exactly one model in {path}; found {len(structure)}")
    grouped: dict[str, list[Any]] = {}
    owners: dict[str, tuple[int, str]] = {}
    for index, chain in enumerate(structure[0]):
        for residue in chain:
            info = gemmi.find_tabulated_residue(residue.name)
            if not info.is_amino_acid():
                continue
            if chain_id_type == "auth":
                chain_id = chain.name
            elif chain_id_type == "label":
                chain_id = residue.subchain
            else:
                raise ValueError(f"Unsupported chain ID type: {chain_id_type!r}")
            if not chain_id:
                raise ValueError(
                    f"Residue {chain.name}:{residue.seqid} in {path} has no {chain_id_type} chain ID"
                )
            owner = (index, residue.subchain)
            if owners.setdefault(chain_id, owner) != owner:
                raise ValueError(f"Ambiguous {chain_id_type} chain {chain_id!r} in {path}")
            grouped.setdefault(chain_id, []).append(residue)
    return grouped


def parse_fasta_pair(text: str) -> tuple[str, str, str, str]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if len(lines) != 4 or not lines[0].startswith(">reference:") or not lines[2].startswith(">model:"):
        raise ValueError(f"Unexpected OpenStructure alignment format: {text!r}")
    reference_chain = lines[0].split(":", 1)[1]
    model_chain = lines[2].split(":", 1)[1]
    if len(lines[1]) != len(lines[3]):
        raise ValueError(f"Alignment lengths differ for {reference_chain}:{model_chain}")
    return reference_chain, model_chain, lines[1], lines[3]


def residue_pairs(
    ref_chain: str,
    mdl_chain: str,
    refs: list[Any],
    mdls: list[Any],
    alignment: tuple[str, str] | None,
) -> tuple[list[tuple[Any, Any]], dict[str, Any]]:
    if alignment is None:
        raise ValueError(f"OpenStructure JSON has no alignment for mapping {ref_chain}:{mdl_chain}")
    if not refs or not mdls:
        raise ValueError(f"No observed residues for mapping {ref_chain}:{mdl_chain}")
    ref_aln, mdl_aln = alignment
    for chain, aligned, residues in ((ref_chain, ref_aln, refs), (mdl_chain, mdl_aln, mdls)):
        observed = sum(char != "-" for char in aligned)
        if observed != len(residues):
            raise ValueError(
                f"OpenStructure alignment residue count mismatch for {chain}: "
                f"alignment={observed}, structure={len(residues)}"
            )
    refs_iter, mdls_iter = iter(refs), iter(mdls)
    pairs, mismatches, unmatched = [], [], []
    for ref_char, mdl_char in zip(ref_aln, mdl_aln):
        ref_residue = next(refs_iter) if ref_char != "-" else None
        mdl_residue = next(mdls_iter) if mdl_char != "-" else None
        if ref_residue is None:
            continue
        if mdl_residue is None:
            unmatched.append(str(ref_residue.seqid))
            continue
        pairs.append((ref_residue, mdl_residue))
        if ref_char.upper() != mdl_char.upper():
            mismatches.append({
                "reference_residue": str(ref_residue.seqid),
                "model_residue": str(mdl_residue.seqid),
                "reference_code": ref_char,
                "model_code": mdl_char,
            })
    return pairs, {
        "model_chain": mdl_chain,
        "reference_residue_count": len(refs),
        "model_residue_count": len(mdls),
        "aligned_residue_pairs": len(pairs),
        "unmatched_reference_residues": unmatched,
        "sequence_mismatches": mismatches,
    }


def atom_for(residue: Any, name: str) -> Any | None:
    candidates = [
        atom for atom in residue
        if atom.name.strip() == name and atom.occ > 0
        and all(math.isfinite(coord) for coord in (atom.pos.x, atom.pos.y, atom.pos.z))
    ]
    if not candidates:
        return None
    blank = [atom for atom in candidates if not atom.altloc.strip("\x00 ")]
    pool = blank or candidates
    return max(pool, key=lambda atom: atom.occ)


def coordinates(
    pairs_by_ref_chain: dict[str, list[tuple[Any, Any]]],
    chains: list[str],
    atom_names: tuple[str, ...],
    alignment_details: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    reference_xyz: list[list[float]] = []
    model_xyz: list[list[float]] = []
    omitted: list[dict[str, str]] = []
    aligned_residue_pairs = 0
    reference_residue_count = sum(alignment_details[c]["reference_residue_count"] for c in chains)
    for chain in chains:
        for ref_residue, mdl_residue in pairs_by_ref_chain.get(chain, []):
            aligned_residue_pairs += 1
            for atom_name in atom_names:
                ref_atom = atom_for(ref_residue, atom_name)
                mdl_atom = atom_for(mdl_residue, atom_name)
                if ref_atom is None or mdl_atom is None:
                    omitted.append(
                        {
                            "reference_chain": chain,
                            "reference_residue": str(ref_residue.seqid),
                            "model_residue": str(mdl_residue.seqid),
                            "atom": atom_name,
                            "reason": "missing_reference" if ref_atom is None else "missing_model",
                        }
                    )
                    continue
                reference_xyz.append([ref_atom.pos.x, ref_atom.pos.y, ref_atom.pos.z])
                model_xyz.append([mdl_atom.pos.x, mdl_atom.pos.y, mdl_atom.pos.z])
    possible = aligned_residue_pairs * len(atom_names)
    used = len(reference_xyz)
    details = {
        "chains": chains,
        "atom_names": list(atom_names),
        "aligned_residue_pairs": aligned_residue_pairs,
        "reference_residue_count": reference_residue_count,
        "reference_residue_coverage": (
            aligned_residue_pairs / reference_residue_count if reference_residue_count else None
        ),
        "possible_atom_pairs": possible,
        "used_atom_pairs": used,
        "atom_pair_coverage": (used / possible) if possible else None,
        "omitted_atoms": omitted,
    }
    return np.asarray(reference_xyz, dtype=float), np.asarray(model_xyz, dtype=float), details


def kabsch(model_xyz: np.ndarray, reference_xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if len(model_xyz) < 3:
        raise ValueError("Fewer than three atom pairs remain")
    model_centered = model_xyz - model_xyz.mean(axis=0)
    reference_centered = reference_xyz - reference_xyz.mean(axis=0)
    if np.linalg.matrix_rank(model_centered) < 2 or np.linalg.matrix_rank(reference_centered) < 2:
        raise ValueError("Atom pairs are collinear")
    u, _, vt = np.linalg.svd(model_centered.T @ reference_centered)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    translation = reference_xyz.mean(axis=0) - model_xyz.mean(axis=0) @ rotation
    return rotation, translation


def rmsd(reference_xyz: np.ndarray, transformed_model_xyz: np.ndarray) -> float:
    return math.sqrt(float(np.mean(np.sum((transformed_model_xyz - reference_xyz) ** 2, axis=1))))


def measure(
    pairs_by_ref_chain: dict[str, list[tuple[Any, Any]]],
    chains: list[str],
    atom_names: tuple[str, ...],
    alignment_details: dict[str, Any],
    fit: tuple[np.ndarray, np.ndarray] | None = None,
    fit_failure: str | None = None,
) -> dict[str, Any]:
    reference_xyz, model_xyz, details = coordinates(
        pairs_by_ref_chain, chains, atom_names, alignment_details
    )
    try:
        errors = [alignment_details[chain]["reason"] for chain in chains
                  if alignment_details[chain].get("reason")]
        if errors:
            raise ValueError("; ".join(errors))
        if fit_failure:
            raise ValueError(fit_failure)
        if not len(reference_xyz):
            raise ValueError("No matched atom pairs remain")
        rotation, translation = fit if fit is not None else kabsch(model_xyz, reference_xyz)
    except ValueError as exc:
        return {"status": "UNAVAILABLE", "reason": str(exc), **details}
    value = rmsd(reference_xyz, model_xyz @ rotation + translation)
    return {
        "status": "SUCCESS",
        "rmsd_angstrom": value,
        "fit_method": ("target fit" if fit is not None
                       else "independent unweighted Kabsch fit; no outlier rejection"),
        **details,
    }


def compute(config: dict[str, Any], ost: dict[str, Any]) -> dict[str, Any]:
    reference = Path(config["reference"])
    model = Path(config["model"])
    mapping = config["chain_mapping"]
    target_chains = config["target_chains"]
    binders = config["binders"]
    reference_id_type = config["reference_chain_id_type"]
    model_id_type = config["model_chain_id_type"]

    reference_residues = chain_residues(reference, reference_id_type)
    model_residues = chain_residues(model, model_id_type)
    needed = set(target_chains)
    for binder in binders:
        needed.update(binder["chains"])
    alignments = {}
    malformed_alignments = []
    for text in ost.get("aln", []):
        try:
            ref_chain, mdl_chain, ref_aln, mdl_aln = parse_fasta_pair(text)
        except ValueError as exc:
            malformed_alignments.append(str(exc))
        else:
            alignments[(ref_chain, mdl_chain)] = (ref_aln, mdl_aln)
    pairs, alignment_details = {}, {}
    for chain in sorted(needed):
        model_chain = mapping[chain]
        residues = reference_residues.get(chain, [])
        try:
            pairs[chain], alignment_details[chain] = residue_pairs(
                chain, model_chain, residues, model_residues.get(model_chain, []),
                alignments.get((chain, model_chain)),
            )
        except ValueError as exc:
            alignment_details[chain] = {
                "model_chain": model_chain,
                "reference_residue_count": len(residues),
                "aligned_residue_pairs": 0,
                "unmatched_reference_residues": [str(residue.seqid) for residue in residues],
                "reason": str(exc),
            }

    ref_target, mdl_target, target_details = coordinates(
        pairs, target_chains, BACKBONE, alignment_details
    )
    target_fit, target_failure, target_rmsd = None, None, None
    try:
        errors = [alignment_details[chain]["reason"] for chain in target_chains
                  if alignment_details[chain].get("reason")]
        if errors:
            raise ValueError("; ".join(errors))
        target_fit = kabsch(mdl_target, ref_target)
        rotation, translation = target_fit
        target_rmsd = rmsd(ref_target, mdl_target @ rotation + translation)
    except ValueError as exc:
        target_failure = str(exc)

    target_aligned_binders = []
    self_aligned_binders = []
    for binder in binders:
        chains = binder["chains"]
        identity = {
            "id": binder["id"],
            "reference_chains": chains,
            "model_chains": [mapping[chain] for chain in chains],
        }
        target_aligned_binders.append(
            {
                **identity,
                "target_aligned_binder_ca_rmsd_angstrom": measure(
                    pairs, chains, ("CA",), alignment_details, target_fit, target_failure
                ),
                "target_aligned_binder_backbone_rmsd_angstrom": measure(
                    pairs, chains, BACKBONE, alignment_details, target_fit, target_failure
                ),
            }
        )
        self_aligned_binders.append(
            {
                **identity,
                "binder_self_aligned_ca_rmsd_angstrom": measure(
                    pairs, chains, ("CA",), alignment_details
                ),
                "binder_self_aligned_backbone_rmsd_angstrom": measure(
                    pairs, chains, BACKBONE, alignment_details
                ),
            }
        )

    target_status = (
        "SUCCESS"
        if target_failure is None and all(
            result.get("status") == "SUCCESS"
            for binder in target_aligned_binders
            for result in (
                binder["target_aligned_binder_ca_rmsd_angstrom"],
                binder["target_aligned_binder_backbone_rmsd_angstrom"],
            )
        )
        else "PARTIAL"
    )
    self_status = (
        "SUCCESS"
        if all(
            result.get("status") == "SUCCESS"
            for binder in self_aligned_binders
            for result in (
                binder["binder_self_aligned_ca_rmsd_angstrom"],
                binder["binder_self_aligned_backbone_rmsd_angstrom"],
            )
        )
        else "PARTIAL"
    )

    return {
        "schema_version": 2,
        "status": "SUCCESS" if target_status == self_status == "SUCCESS" else "PARTIAL",
        "metric_family": "binder_rmsd",
        "units": "angstrom",
        "reference": str(reference),
        "model": str(model),
        "chain_mapping": mapping,
        "target_aligned": {
            "status": target_status,
            "alignment": {
                "status": "UNAVAILABLE" if target_failure else "SUCCESS",
                "reason": target_failure,
                "reference_chains": target_chains,
                "model_chains": [mapping[chain] for chain in target_chains],
                "atoms": list(BACKBONE),
                "method": "single unweighted Kabsch fit; no outlier rejection",
                "target_backbone_rmsd_angstrom": target_rmsd,
                **target_details,
            },
            "binders": target_aligned_binders,
        },
        "self_aligned": {
            "status": self_status,
            "scope": "Each annotated binder group is fitted independently of the target",
            "multichain_policy": "All chains in one binder group are fitted jointly",
            "atom_set_policy": "CA and backbone metrics each use their own optimal fit",
            "binders": self_aligned_binders,
        },
        "residue_correspondence": {
            "source": "OpenStructure pairwise alignments",
            "sequence_direction_reversal_allowed": False,
            "malformed_alignments": malformed_alignments,
            "chains": alignment_details,
        },
    }


def main() -> int:
    args = parse_args()
    try:
        result = compute(load_json(args.config), load_json(args.ost_json))
    except Exception as exc:  # Preserve a machine-readable failure artifact.
        result = {
            "schema_version": 1,
            "status": "UNAVAILABLE",
            "reason": f"{type(exc).__name__}: {exc}",
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return 0 if result["status"] == "SUCCESS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
