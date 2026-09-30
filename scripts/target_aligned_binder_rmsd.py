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
    for chain in structure[0]:
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


def residue_pairs_from_ost(
    ost: dict[str, Any],
    reference_residues: dict[str, list[Any]],
    model_residues: dict[str, list[Any]],
    required_mapping: dict[str, str],
) -> tuple[dict[str, list[tuple[Any, Any]]], dict[str, Any]]:
    alignments: dict[tuple[str, str], tuple[str, str]] = {}
    for text in ost.get("aln", []):
        ref_chain, mdl_chain, ref_aln, mdl_aln = parse_fasta_pair(text)
        alignments[(ref_chain, mdl_chain)] = (ref_aln, mdl_aln)

    pairs_by_ref_chain: dict[str, list[tuple[Any, Any]]] = {}
    diagnostics: dict[str, Any] = {}
    for ref_chain, mdl_chain in required_mapping.items():
        key = (ref_chain, mdl_chain)
        if key not in alignments:
            raise ValueError(f"OpenStructure JSON has no alignment for mapping {ref_chain}:{mdl_chain}")
        if ref_chain not in reference_residues:
            raise ValueError(f"Reference chain {ref_chain!r} is absent")
        if mdl_chain not in model_residues:
            raise ValueError(f"Model chain {mdl_chain!r} is absent")
        ref_aln, mdl_aln = alignments[key]
        refs = reference_residues[ref_chain]
        mdls = model_residues[mdl_chain]
        if sum(char != "-" for char in ref_aln) != len(refs):
            raise ValueError(
                f"OpenStructure alignment/reference residue count mismatch for {ref_chain}: "
                f"alignment={sum(char != '-' for char in ref_aln)}, structure={len(refs)}"
            )
        if sum(char != "-" for char in mdl_aln) != len(mdls):
            raise ValueError(
                f"OpenStructure alignment/model residue count mismatch for {mdl_chain}: "
                f"alignment={sum(char != '-' for char in mdl_aln)}, structure={len(mdls)}"
            )

        ref_index = model_index = 0
        pairs: list[tuple[Any, Any]] = []
        mismatches: list[dict[str, str]] = []
        for ref_char, mdl_char in zip(ref_aln, mdl_aln):
            ref_residue = refs[ref_index] if ref_char != "-" else None
            mdl_residue = mdls[model_index] if mdl_char != "-" else None
            if ref_char != "-":
                ref_index += 1
            if mdl_char != "-":
                model_index += 1
            if ref_residue is not None and mdl_residue is not None:
                pairs.append((ref_residue, mdl_residue))
                if ref_char.upper() != mdl_char.upper():
                    mismatches.append(
                        {
                            "reference_residue": str(ref_residue.seqid),
                            "model_residue": str(mdl_residue.seqid),
                            "reference_code": ref_char,
                            "model_code": mdl_char,
                        }
                    )
        pairs_by_ref_chain[ref_chain] = pairs
        diagnostics[ref_chain] = {
            "model_chain": mdl_chain,
            "reference_residue_count": len(refs),
            "model_residue_count": len(mdls),
            "aligned_residue_pairs": len(pairs),
            "sequence_mismatches": mismatches,
        }
    return pairs_by_ref_chain, diagnostics


def atom_for(residue: Any, name: str) -> Any | None:
    candidates = [atom for atom in residue if atom.name.strip() == name and atom.occ > 0]
    if not candidates:
        return None
    blank = [atom for atom in candidates if not atom.altloc.strip("\x00 ")]
    pool = blank or candidates
    return max(pool, key=lambda atom: atom.occ)


def coordinates(
    pairs_by_ref_chain: dict[str, list[tuple[Any, Any]]],
    chains: list[str],
    atom_names: tuple[str, ...],
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    reference_xyz: list[list[float]] = []
    model_xyz: list[list[float]] = []
    omitted: list[dict[str, str]] = []
    aligned_residue_pairs = 0
    for chain in chains:
        if chain not in pairs_by_ref_chain:
            raise ValueError(f"No residue pairs for reference chain {chain!r}")
        for ref_residue, mdl_residue in pairs_by_ref_chain[chain]:
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


def measured_result(
    pairs_by_ref_chain: dict[str, list[tuple[Any, Any]]],
    chains: list[str],
    atom_names: tuple[str, ...],
    rotation: np.ndarray,
    translation: np.ndarray,
) -> dict[str, Any]:
    reference_xyz, model_xyz, details = coordinates(pairs_by_ref_chain, chains, atom_names)
    if not len(reference_xyz):
        return {"status": "UNAVAILABLE", "reason": "No matched atom pairs remain", **details}
    value = rmsd(reference_xyz, model_xyz @ rotation + translation)
    return {"status": "SUCCESS", "rmsd_angstrom": value, **details}


def self_aligned_result(
    pairs_by_ref_chain: dict[str, list[tuple[Any, Any]]],
    chains: list[str],
    atom_names: tuple[str, ...],
) -> dict[str, Any]:
    reference_xyz, model_xyz, details = coordinates(pairs_by_ref_chain, chains, atom_names)
    try:
        rotation, translation = kabsch(model_xyz, reference_xyz)
    except ValueError as exc:
        return {"status": "UNAVAILABLE", "reason": str(exc), **details}
    value = rmsd(reference_xyz, model_xyz @ rotation + translation)
    return {
        "status": "SUCCESS",
        "rmsd_angstrom": value,
        "fit_method": "independent unweighted Kabsch fit; no outlier rejection",
        **details,
    }


def compute(config: dict[str, Any], ost: dict[str, Any]) -> dict[str, Any]:
    if ost.get("status") != "SUCCESS":
        raise ValueError("OpenStructure result is not successful")
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
    needed_mapping = {chain: mapping[chain] for chain in needed}
    pairs, alignment_details = residue_pairs_from_ost(
        ost, reference_residues, model_residues, needed_mapping
    )

    ref_target, mdl_target, target_details = coordinates(pairs, target_chains, BACKBONE)
    rotation, translation = kabsch(mdl_target, ref_target)
    transformed_target = mdl_target @ rotation + translation
    target_rmsd = rmsd(ref_target, transformed_target)

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
                "target_aligned_binder_ca_rmsd_angstrom": measured_result(
                    pairs, chains, ("CA",), rotation, translation
                ),
                "target_aligned_binder_backbone_rmsd_angstrom": measured_result(
                    pairs, chains, BACKBONE, rotation, translation
                ),
            }
        )
        self_aligned_binders.append(
            {
                **identity,
                "binder_self_aligned_ca_rmsd_angstrom": self_aligned_result(
                    pairs, chains, ("CA",)
                ),
                "binder_self_aligned_backbone_rmsd_angstrom": self_aligned_result(
                    pairs, chains, BACKBONE
                ),
            }
        )

    target_status = (
        "SUCCESS"
        if all(
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
