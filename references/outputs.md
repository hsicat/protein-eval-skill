# Outputs and interpretation

Each prediction path is preserved below `results/`, including its filename, followed by `runs/<fingerprint-prefix>/`. Thus predictions with identical sample IDs or basenames cannot collide.

Each batch has one content-addressed manifest at
`results/batches/<batch-fingerprint-prefix>/run_manifest.json`. It records the
annotation file and hash, complete annotations, shared tool versions, and each
planned or skipped prediction with its input hashes, settings, run fingerprint,
and result directory. An identical batch reuses this manifest instead of
rewriting it.

Each run contains:

- `reference_scored_chains.cif`, `model_scored_chains.cif`, and `subset_manifest.json`: exact polymer-only tool inputs containing target and binder chains but excluding context and non-polymer residues.
- `commands.txt`: always retained with exact commands, paths, environments/versions, and exit codes.
- Diagnostic stream files are conditional: failed commands retain both
  `<tool>.stdout.txt` and `<tool>.stderr.txt`, including empty files; successful
  commands retain only a nonempty stderr file and discard stdout.
- `ost_scores.json`: untouched OpenStructure JSON.
- `role_metrics_config.json` and `role_metrics.json`: annotated-role
  OpenStructure lDDT and target-binder iLDDT results, including per-pair
  distance-test counts. The role-specific target-binder DockQ subset is stored
  in `summary.json` under `role_metrics.target_binder_dockq`.
- `custom_rmsd_config.json` and `binder_rmsd.json`.
- `chirality/`: the chirality tool's readable report, simple CSV, detailed CSV, and raw summary JSON.
- `summary.json`: parsed per-prediction summary retaining local lDDT, mappings,
  alignments, residue inconsistencies, cleanup/stereochemical diagnostics,
  per-interface identities and scores, raw and role-scoped DockQ, custom
  RMSDs, chirality summaries, statuses, warnings, and unavailable reasons.

`results/batch_summary.csv` has one row per prediction/binder for the fingerprints
selected by the latest batch execution, plus explicit skipped rows. Historical
content-addressed runs remain on disk. Per-interface OpenStructure results remain
present as JSON in the CSV and as structured data in each summary; never collapse
them to only one DockQ value.
The CSV exposes pooled all-residue chirality as the primary result and also retains
separate DAA and LAA counts and accuracies. It also includes the applicable
target self-aligned, binder self-aligned, and target-aligned binder-backbone
RMSD cutoffs plus the per-interface DockQ cutoff. Its DockQ success field is a
nullable fraction per prediction; the older raw `ost_dockq_ave_above_threshold`
Boolean is retained only as a diagnostic.

After the CSV is written, `scripts/summarize_batch_metrics.py` produces:

- `results/method_summary.json`: machine-readable method-level aggregations.
- `results/method_summary.md`: Folding and Docking Markdown tables that must be
  included in the final user reply.

The Folding table reports median target self-aligned backbone RMSD, median
binder self-aligned backbone RMSD, mean target-only and binder-only lDDT, and residue-pooled chirality
violations. Target self-aligned backbone RMSD is sourced from
`target_backbone_fit_rmsd_angstrom`: the target N/CA/C/O residual after the
joint all-target fit. It also reports independent success proportions for the
two folding RMSDs. Binder-only lDDT is averaged across eligible binder groups
once within a prediction. The Docking table reports mean target-binder DockQ,
median target-aligned binder backbone RMSD, median target-binder per-interface
iRMSD, mean macro target-binder iLDDT, and independent DockQ and binder RMSD
success rates. A supplemental table retains mean full-subset lDDT,
target-binder DockQ-wave, and micro target-binder iLDDT.
Prediction-level scores are deduplicated by prediction path before aggregation;
binder RMSDs use one value per annotated binder group, and iRMSD uses one value
per eligible target-binder interface. DockQ success is first classified per
target-binder interface, then averaged within each prediction and across
predictions. Each prediction has equal weight in that dataset mean. Target
self-aligned RMSD success is counted once per prediction. Each binder-RMSD success
is counted once per annotated binder group. Unavailable scores are excluded
from the corresponding denominator. Explicitly skipped predictions are excluded.

## Metric names

- `openstructure.lddt`: global all-atom lDDT.
- `openstructure.ilddt`: global all-atom lDDT restricted to inter-chain contacts.
- `role_metrics.target_only_lddt`: all-atom lDDT on all annotated target chains.
- `role_metrics.binder_only_lddt[*].lddt`: all-atom lDDT on one annotated
  binder group; eligible groups are averaged once per prediction for the method
  summary.
- `role_metrics.target_binder_dockq.per_interface`: only native interfaces
  crossing annotated target and binder chains, with original DockQ components.
  `mean` is the equal per-interface mean; `wave` uses native-contact counts.
- `role_metrics.target_binder_ilddt.per_chain_pair`: pairwise OpenStructure
  inter-chain lDDT and exact expected/conserved distance-test counts. `macro`
  is the equal mean of eligible pair scores, and `micro` is the ratio of pooled
  conserved to expected tests within that prediction. Pairs with zero expected
  tests are excluded from both means and reported with `ilddt: null`.
- `openstructure.local_lddt`: reference-residue keyed all-atom local lDDT.
- `openstructure.per_interface[*].lrmsd`: OpenStructure's DockQ ligand RMSD for
  one chain-pair interface. The larger chain is treated as receptor and fitted;
  RMSD is measured on the smaller chain using mapped backbone atoms. This is not
  the custom grouped-binder metric.
- `target_backbone_fit_rmsd_angstrom`: target RMSD after the joint all-target N/CA/C/O fit.
- `target_self_aligned_backbone_rmsd_below_threshold`: nullable Boolean from
  the strict target-folding RMSD cutoff, which defaults to `< 2 Å`.
- `target_aligned_binder_ca_rmsd_angstrom`: mapped binder CA displacement after the target fit, without binder refitting.
- `target_aligned_binder_backbone_rmsd_angstrom`: mapped binder N/CA/C/O displacement after the target fit, without binder refitting.
- `target_aligned_binder_backbone_rmsd_below_threshold`: nullable Boolean from
  the strict RMSD cutoff, which defaults to `< 2 Å`.
- `target_binder_dockq_success_fraction`: fraction of scored target-binder
  interfaces passing strict `DockQ > 0.23`; null when none can be scored.
- `ost_dockq_ave_above_threshold`: legacy diagnostic Boolean using the raw
  all-scored-chain average; it is not the headline DockQ success field.
- `binder_self_aligned_ca_rmsd_angstrom`: CA RMSD after independently fitting
  that binder group on its CA atoms.
- `binder_self_aligned_backbone_rmsd_angstrom`: N/CA/C/O RMSD after independently
  fitting that binder group on its N/CA/C/O atoms.
- `binder_self_aligned_backbone_rmsd_below_threshold`: nullable Boolean from
  the strict binder-folding RMSD cutoff, which defaults to `< 2 Å`.

Missing residues or atoms are excluded only through the recorded aligned intersection. Always report used/possible atom pairs, coverage, and the omission list alongside a custom RMSD. If the target fit has fewer than three non-collinear atom pairs or a binder measurement has no atom pairs, that metric is unavailable with a reason.

A score exactly equal to its cutoff is a failure because both comparisons are
strict. A missing or non-finite score is unavailable, not a failure. Custom
cutoffs supplied to the batch runner are stored in the batch fingerprint and
manifest, CSV, and method-level summary. They are intentionally excluded from
prediction-run fingerprints so changing only a reporting cutoff can reuse the
existing structure scores.

## Reruns

The fingerprint includes reference and prediction file hashes, complete sample/prediction annotations, settings, inspected script hashes, and tool versions. An existing successful fingerprint is reused. A changed fingerprint creates another run and preserves the earlier result. An incomplete fingerprint is not resumed unless the operator reviews it and passes `--resume-incomplete`.

The batch fingerprint covers the complete annotation document, folding and
docking classification cutoffs, shared tool versions, and the ordered
planned/skipped prediction records. A changed batch creates a new directory
below `results/batches/`; historical batch manifests remain available.
