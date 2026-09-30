#!/usr/bin/env python3
"""Print chain-ID aliases, sequences, and non-polymers for structure auditing."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import gemmi


def inspect(path: Path) -> dict:
    structure = gemmi.read_structure(str(path))
    models = []
    for model_index, model in enumerate(structure, start=1):
        chains = []
        nonpolymers = []
        for chain in model:
            by_label: dict[str, list] = {}
            for residue in chain:
                info = gemmi.find_tabulated_residue(residue.name)
                if info.is_amino_acid():
                    by_label.setdefault(residue.subchain or "", []).append(residue)
                elif residue.name not in {"HOH", "WAT", "DOD"}:
                    nonpolymers.append(
                        {
                            "auth_chain": chain.name,
                            "label_chain": residue.subchain,
                            "residue": str(residue.seqid),
                            "name": residue.name,
                        }
                    )
            for label_chain, residues in by_label.items():
                chains.append(
                    {
                        "auth_chain": chain.name,
                        "label_chain": label_chain,
                        "observed_residue_count": len(residues),
                        "sequence": "".join(
                            gemmi.find_tabulated_residue(residue.name).one_letter_code
                            for residue in residues
                        ),
                        "first_residue": str(residues[0].seqid),
                        "last_residue": str(residues[-1].seqid),
                    }
                )
        models.append({"model": model_index, "chains": chains, "nonpolymers": nonpolymers})
    return {
        "path": str(path.resolve()),
        "format": path.suffix.lower(),
        "model_count": len(structure),
        "embedded_assemblies": [assembly.name for assembly in structure.assemblies],
        "models": models,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structures", type=Path, nargs="+")
    args = parser.parse_args()
    result = [inspect(path) for path in args.structures]
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
