# v0.6.0 protocol update — 2026-11-11

## Scope

This update changes how selected metrics are **defined inside each predicted
complex**. It does not change the existing dataset-wise aggregation operators
or units: prediction-level bounded scores are arithmetic means across
predictions, RMSDs retain their existing medians across prediction, binder
group, or interface values, and chirality remains residue-pooled. A metric's
reported dataset value can still change because its per-complex input changed.

## Four metric changes

1. Docking metrics use annotated target–binder chain interfaces. Target–target
   and binder–binder interfaces remain in the untouched OpenStructure JSON but
   no longer contribute to headline DockQ or iRMSD. The original full-subset
   OpenStructure scores remain in `summary.json` for diagnostics.
2. DockQ success applies strict `DockQ > 0.23` separately to each eligible
   target–binder interface. Each prediction contributes its fraction of
   successful interfaces to the method mean. This replaces testing whether a
   prediction's mean DockQ exceeds the cutoff. Continuous per-interface DockQ,
   equal interface mean DockQ, and native-contact-weighted DockQ-wave are
   retained.
3. Folding lDDT is computed separately for all annotated target chains and
   for each annotated binder group. Eligible binder-group lDDTs are averaged
   once per prediction before the existing method-level arithmetic mean. The
   original target-plus-binder global lDDT is supplemental.
4. Target–binder iLDDT is computed with OpenStructure for each annotated
   target-chain/binder-chain pair that has native inter-chain distance tests.
   `macro` is the equal mean of those pair scores inside the prediction;
   `micro` is conserved distance tests divided by expected distance tests
   pooled over those pairs inside the prediction. Macro is the headline; micro
   and the original full-subset iLDDT are supplemental. A pair with no native
   distance tests is unavailable, not zero.

## Dataset aggregation deliberately unchanged

| Metric | Existing method-level aggregation retained |
|---|---|
| Target self-aligned backbone RMSD | Median across available predictions |
| Binder self-aligned and target-aligned RMSD | Median across available binder groups |
| Target-only, binder-only, and global lDDT | Arithmetic mean of available prediction values |
| Target–binder DockQ mean, DockQ-wave, iLDDT macro/micro | Arithmetic mean of available prediction values |
| Target–binder iRMSD | Median across available interface values |
| DockQ success | Arithmetic mean of per-prediction success fractions |
| Chirality violation | Pooled assessable residues across the dataset |

The new DockQ success fraction is one value per prediction; complex A with two
interfaces and complex B with one interface each contributes one value to the
method mean. This does **not** implement a general change to equal-complex
weighting for binder RMSD or iRMSD. That possible later revision remains open.

## Output and reproducibility

- New per-run `role_metrics_config.json` and `role_metrics.json` contain the
  role-specific OpenStructure calculations. `summary.json` contains those
  results plus the target–binder DockQ subset and retains raw OST scores.
- `batch_summary.csv` adds target/binder lDDT, target–binder DockQ and iLDDT,
  per-interface JSON, distance-aware results, and the per-prediction DockQ
  success fraction. The old raw `ost_*` fields remain diagnostic.
- `method_summary.md` now has Folding, Docking, and Supplemental tables.
  `method_summary.json` uses schema version 4.
- The role-scoring script hash and changed settings enter the prediction-run
  fingerprint. Prior successful runs remain intact; a new execution creates a
  distinct run rather than silently rewriting old scores.
- The locally installed version had `--method-name` and independent-copy
  binder guidance that the GitHub checkout lacked. Those existing changes are
  synchronized into this update without altering their behavior.

OpenStructure remains version 2.12.0. RMSD fitting, chain mapping, chirality,
and the four numeric threshold defaults are unchanged.
