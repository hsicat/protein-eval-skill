# Annotation schema

Store annotations at `evaluation/annotations/samples.json` unless the user chooses another persistent path. All stored paths are relative to the `evaluation/` directory; absolute paths are allowed but reduce portability. Prediction keys must be exact paths discovered beneath `predicted/`.

```json
{
  "schema_version": 1,
  "samples": {
    "SAMPLE_ID": {
      "reference": "../GT/SAMPLE_ID.cif",
      "reference_chain_id_type": "label",
      "reference_assembly": {
        "selection": "biological_assembly",
        "id": "1",
        "notes": "Materialized assembly file from the project GT manifest"
      },
      "target_chains": ["T1", "T2"],
      "binders": [
        {"id": "binder_1", "chains": ["B1", "B2"]}
      ],
      "context_chains": ["X"],
      "chirality_chains": ["B1", "B2"],
      "nonpolymer_ligands": [],
      "skip": false,
      "predictions": {
        "predicted/method/SAMPLE_ID/model_0.cif": {
          "model_chain_id_type": "label",
          "chain_mapping": {
            "T1": "A",
            "T2": "B",
            "B1": "C",
            "B2": "D"
          },
          "mapping_notes": "Validated from sequences and structure; equivalent copies resolved as documented"
        }
      }
    }
  }
}
```

## Required decisions

- `reference` must identify the exact file to score. Do not replace it with an asymmetric unit or another biological assembly silently.
- `reference_assembly.selection` is one of `biological_assembly`, `asymmetric_unit`, or `prepared_subset`. Record an ID and provenance/notes when available. Prefer a materialized exact assembly or prepared subset so every tool evaluates identical coordinates.
- `reference_chain_id_type` and `model_chain_id_type` are `label` or `auth`. Inspect both namespaces. The mapping keys and values use these declared namespaces and must also match the IDs accepted by OpenStructure.
- `target_chains` form one alignment group.
- Each `binders` entry is evaluated as one binder group and may contain multiple chains. Multiple binders require multiple entries.
- `context_chains` are documented but not target or binder chains.
- `chirality_chains` explicitly defines the biological scope of chirality scoring; do not silently default it to all chains.
- `nonpolymer_ligands` remains empty until the user explicitly identifies ligand instances for separate ligand evaluation.
- Set `skip` to `true` only for an intentional, user-confirmed exclusion and add a
  nonempty `skip_reason`. Skipped predictions remain visible in the audit and batch
  summary but no metric commands run for them. A prediction-level `skip` and
  `skip_reason` may be used when only one model file should be excluded.

The same reference chain cannot occur in multiple biological roles. Chain mappings must be one-to-one for the scored chains. If symmetry-equivalent chains permit multiple score-changing mappings, document the selected mapping or ask the user.
