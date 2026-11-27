# RCG eval_v01 reference-to-prediction chain mappings

These files are optional examples of the mappings used for this benchmark.
For regular evaluations, use your evaluation directory's own
`annotations/samples.json` or explicitly choose a file with `--annotations`.
Audit the annotations against your own inputs.

Each method directory contains the byte-for-byte `annotations/samples.json`
used in its recorded `eval_v01` batch and an SHA-256 checklist for the predicted
structures and ground-truth CIFs it references. The JSON preserves the chain ID
namespaces, exact assembly choices, target and binder groups, per-prediction
chain mappings, and skips. Each `chain_mapping` maps reference chain IDs to
prediction chain IDs.

| Method directory | Mapped predictions | Intentional skips | Source batch ID |
|---|---:|---:|---|
| `af3` | 195 | 0 | `84d8b0948925c4d9` |
| `af3_rcg` | 195 | 0 | `1a832230377560ac` |
| `af3_rcg_steering` | 195 | 0 | `5c3cc8c15db4aa00` |
| `af3_steering` | 195 | 5 | `463964f611eb4188` |
| `boltz1` | 195 | 0 | `8f800986102ace22` |
| `boltz1x` | 195 | 0 | `e3552ca31acc6918` |
| `boltz2` | 195 | 0 | `5c03eb9c36c45e8d` |
| `boltz2x` | 195 | 0 | `28db41a60b36e370` |

The source batch IDs identify `<method>/results/batches/<id>/run_manifest.json`
in `eval_v01`. The five skipped AF3 steering predictions are for `3V0X` and
retain their original skip reasons.

## Optional reuse for exact benchmark inputs

Use the exact `predicted/` layout and structure files from the matching method
in `eval_v01`. In that method directory, `references/GT/` must resolve to the
frozen GT CIFs. The checksum paths and the paths inside `samples.json` are both
relative to the method directory.

For example, from the `af3_rcg` evaluation directory, with `skill_dir` set to the
absolute path of your skill installation or checkout:

```bash
method=af3_rcg
mapping_dir="$skill_dir/ref_mapping_1127/$method"
shasum -a 256 -c "$mapping_dir/inputs.sha256" && \
  conda run -n urop python "$skill_dir/scripts/batch_evaluate.py" "$PWD" \
    --annotations "$mapping_dir/samples.json" --method-name "$method"
```

On Linux, `sha256sum -c` is equivalent to `shasum -a 256 -c`. Proceed with
`--execute` after checksum verification succeeds, the audit reports `READY`,
and any warnings have been reviewed. If a checksum differs, locate the exact
input used for the source run before applying the frozen mapping.

For historical score reproduction, follow `eval_v01/REPRODUCE.md` for the six
original methods and `eval_v01/af3_rcg_evaluation/REPRODUCE.md` for the RCG
extension. Each uses its own archived scoring protocol.
