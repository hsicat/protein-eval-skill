#!/usr/bin/env python3
"""Score annotated polymer roles with OpenStructure 2.12.

The full-complex compare-structures output remains the source of chain mapping
and DockQ components. This script scores role-specific lDDT and pairwise
target-binder iLDDT using the same OpenStructure Scorer defaults.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from ost import io, mol
from ost.mol.alg.scoring import Scorer


def selected(entity: Any, chains: list[str]) -> Any:
    query = " or ".join(f"cname={mol.QueryQuoteName(chain)}" for chain in chains)
    result = entity.Select(query)
    present = {chain.GetName() for chain in result.chains}
    if present != set(chains):
        raise ValueError(f"Requested chains {chains}, found {sorted(present)}")
    return result


def scorer_for(
    reference: Any, model: Any, mapping: dict[str, str], min_pep_length: int
) -> Scorer:
    return Scorer(
        selected(model, list(mapping.values())),
        selected(reference, list(mapping)),
        custom_mapping=mapping,
        min_pep_length=min_pep_length,
    )


def pair_distance_counts(scorer: Scorer) -> tuple[int, int]:
    """Use the exact inputs of Scorer._compute_ilddt, retaining its counts."""
    model_chains = {chain.GetName() for chain in scorer.stereochecked_model.chains}
    alignments = {
        alignment.GetSequence(1).GetName(): alignment
        for alignment in scorer.stereochecked_aln
        if alignment.GetSequence(1).GetName() in model_chains
    }
    mapping = {
        model_chain: reference_chain
        for model_chain, reference_chain in scorer.mapping.GetFlatMapping(
            mdl_as_key=True
        ).items()
        if model_chain in alignments
    }
    score, _, expected, conserved, *_ = scorer.lddt_scorer.lDDT(
        scorer.stereochecked_model,
        chain_mapping=mapping,
        residue_mapping=alignments,
        check_resnames=False,
        no_intrachain=True,
        return_dist_test=True,
    )
    public_score = scorer.ilddt
    if score is None and public_score is None:
        return expected, conserved
    if score is None or public_score is None or not math.isclose(
        score, public_score, rel_tol=1e-10, abs_tol=1e-12
    ):
        raise RuntimeError("OpenStructure iLDDT score/count mismatch")
    return expected, conserved


def compute(config: dict[str, Any]) -> dict[str, Any]:
    reference = io.LoadMMCIF(config["reference"])
    model = io.LoadMMCIF(config["model"])
    mapping = config["chain_mapping"]
    targets = config["target_chains"]
    binders = config["binders"]
    min_pep_length = config["min_pep_length"]

    target_scorer = scorer_for(
        reference, model, {chain: mapping[chain] for chain in targets}, min_pep_length
    )
    target_lddt = target_scorer.lddt
    binder_lddt: list[dict[str, Any]] = []
    for binder in binders:
        group_mapping = {chain: mapping[chain] for chain in binder["chains"]}
        group_scorer = scorer_for(reference, model, group_mapping, min_pep_length)
        binder_lddt.append(
            {
                "binder_id": binder["id"],
                "reference_chains": binder["chains"],
                "lddt": group_scorer.lddt,
            }
        )

    pairs: list[dict[str, Any]] = []
    for target_chain in targets:
        for binder in binders:
            for binder_chain in binder["chains"]:
                pair_mapping = {
                    target_chain: mapping[target_chain],
                    binder_chain: mapping[binder_chain],
                }
                pair_scorer = scorer_for(
                    reference, model, pair_mapping, min_pep_length
                )
                expected, conserved = pair_distance_counts(pair_scorer)
                pairs.append(
                    {
                        "reference_chains": [target_chain, binder_chain],
                        "binder_id": binder["id"],
                        "ilddt": pair_scorer.ilddt,
                        "distance_tests_expected": expected,
                        "distance_tests_conserved": conserved,
                    }
                )
    eligible = [pair for pair in pairs if pair["distance_tests_expected"] > 0]
    total_expected = sum(pair["distance_tests_expected"] for pair in eligible)
    total_conserved = sum(pair["distance_tests_conserved"] for pair in eligible)
    return {
        "schema_version": 1,
        "status": "SUCCESS",
        "target_only_lddt": target_lddt,
        "binder_only_lddt": binder_lddt,
        "target_binder_ilddt": {
            "per_chain_pair": pairs,
            "macro": (
                sum(pair["ilddt"] for pair in eligible) / len(eligible)
                if eligible else None
            ),
            "micro": total_conserved / total_expected if total_expected else None,
            "n_eligible_pairs": len(eligible),
            "distance_tests_expected": total_expected,
            "distance_tests_conserved": total_conserved,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    result = compute(config)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
