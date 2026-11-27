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
- `ost_scores_no_stereochecks.json`: untouched OpenStructure JSON for global
  all-atom lDDT, local lDDT, and global iLDDT with the checks disabled.
- `role_metrics_config.json` and `role_metrics.json`: annotated-role
  OpenStructure lDDT and target-binder iLDDT results, including per-pair
  distance-test counts. The role-specific target-binder DockQ subset is stored
  in `summary.json` under `role_metrics.target_binder_dockq`.
- `role_metrics_no_stereochecks.json`: target-only and binder-group all-atom
  lDDT, target–binder chain-pair iLDDT with distance-test counts, and the
  macro and micro iLDDT results with checks disabled.
- `custom_rmsd_config.json` and `binder_rmsd.json`.
- `chirality/`: the chirality tool's readable report, simple CSV, detailed CSV, and raw summary JSON.
- `summary.json`: parsed per-prediction summary retaining local lDDT, mappings,
  alignments, residue inconsistencies, cleanup/stereochemical diagnostics,
  per-interface identities and scores, raw and role-scoped DockQ, custom
  RMSDs, chirality summaries, statuses, warnings, and unavailable reasons.

`results/batch_summary.csv` has one row per prediction/binder from the active
execution's explicit result list, plus rows for predictions skipped in that
execution. Historical
content-addressed runs remain on disk. Per-interface OpenStructure results remain
present as JSON in the CSV and as structured data in each summary; never collapse
them to only one DockQ value.
The CSV exposes pooled All (D+L) chirality as the primary result and retains
separate DAA and LAA counts and accuracies. For each group, `*_eligible_count`
includes unassessable sites, `*_skipped_count` records those sites, and
`*_accuracy_pct` divides `*_correct_count` by `*_assessable_count`. Zero
assessable sites yields NA. The chirality summary JSON includes chain-wise
counts. It also includes the applicable
target self-aligned, binder self-aligned, and target-aligned binder-backbone
RMSD cutoffs plus the per-interface DockQ cutoff. Its DockQ success field is a
nullable fraction per prediction; the older raw `ost_dockq_ave_above_threshold`
Boolean is retained only as a diagnostic.

After the CSV is written, `scripts/summarize_batch_metrics.py` produces:

- `results/method_summary.json`: machine-readable method-level aggregations.
- `results/method_summary.md`: Folding and Docking Markdown tables that must be
  included in the final user reply, along with the Supplemental table.

The Folding table reports median target self-aligned backbone RMSD, median
per-prediction mean binder self-aligned backbone RMSD, mean target-only and
binder-only lDDT, and pooled All (D+L), D-only and L-only chirality accuracy.
Accuracy is passed / assessable, pooling site counts across included
complex–seed predictions. The JSON also retains pooled violations and
unassessable counts. Target self-aligned
backbone RMSD is sourced from
`target_backbone_fit_rmsd_angstrom`: the target N/CA/C/O residual after the
joint all-target fit. It also reports independent success proportions for the
two folding RMSDs. Binder-only lDDT is averaged across eligible binder groups
once within a prediction. The Docking table reports mean target-binder DockQ,
median per-prediction mean L-RMSD, median
per-prediction mean target-binder iRMSD, mean macro target-binder iLDDT, and
independent DockQ and L-RMSD success rates. A supplemental table retains
mean full-subset lDDT,
target-binder DockQ-wave, and micro target-binder iLDDT.
Prediction-level scores are deduplicated by prediction path before aggregation.
Each binder-group RMSD and target-binder interface iRMSD remains available in
the per-prediction outputs. Within each prediction, eligible binder-group or
interface RMSDs are averaged once; the method headline is the median of those
per-prediction means. The earlier group- and interface-pooled medians are
retained under `method_summary.json.methods[*].supplemental_pooled`.
DockQ success is first classified per target-binder interface, then averaged
within each prediction and across
predictions. Each prediction has equal weight in that dataset mean. Target
self-aligned RMSD success is counted once per prediction. Each binder-group RMSD
success is classified separately, averaged to a pass fraction within its
prediction, then averaged across predictions; the original group-pooled rate
is supplemental. Unavailable scores are excluded from the corresponding
denominator; a prediction with no available value for one metric is omitted
only from that metric. Explicitly skipped predictions are excluded. Chirality
violation continues to pool assessable residues across the dataset, with no
complex-equal averaging.

The existing lDDT and iLDDT fields use OpenStructure's default stereochemistry
checks. The disabled-mode scores are stored under
`summary.json.openstructure_no_stereochecks` and
`summary.json.role_metrics_no_stereochecks`, with the same inner score names.
In the CSV and method JSON, append `_no_stereochecks` to the original score
field name (before `_json` for the chain-pair list). Available data and
pair eligibility are evaluated independently in each mode. The tables label
both modes and show NA for missing disabled-mode scores in older results.
Disabling stereochemistry checks preserves OpenStructure's standard input cleanup.

## Metric names

Use **L-RMSD** as the figure/table display name for
`target_aligned_binder_backbone_rmsd_angstrom`, formerly “target-aligned binder
backbone RMSD” or “Binder pose RMSD”. Measure binder N/CA/C/O displacement after
a joint N/CA/C/O fit of all annotated target chains, without refitting the binder.
Use **L-RMSD success** for its threshold rate and **L-RMSD (Å)** on axes.
Keep stored field names, calculations, aggregation, and cutoffs unchanged.
Keep CA-only, self-aligned, and OpenStructure per-interface `lrmsd` results distinct.

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
- `target_aligned_binder_backbone_rmsd_angstrom`: **L-RMSD**; mapped binder N/CA/C/O displacement after the target fit, without binder refitting.
- `target_aligned_binder_backbone_rmsd_below_threshold`: nullable Boolean from
  the strict RMSD cutoff, which defaults to `< 2 Å`.
- `target_binder_dockq_success_fraction`: fraction of scored target-binder
  interfaces passing the batch cutoff, which defaults to strict `DockQ > 0.23`.
  It is null when no interfaces can be scored.
- `ost_dockq_ave_above_threshold`: legacy diagnostic Boolean using the raw
  all-scored-chain average; it is not the headline DockQ success field.
- `binder_self_aligned_ca_rmsd_angstrom`: CA RMSD after independently fitting
  that binder group on its CA atoms.
- `binder_self_aligned_backbone_rmsd_angstrom`: N/CA/C/O RMSD after independently
  fitting that binder group on its N/CA/C/O atoms.
- `binder_self_aligned_backbone_rmsd_below_threshold`: nullable Boolean from
  the strict binder-folding RMSD cutoff, which defaults to `< 2 Å`.

The reference-residue coverage denominator is the number of observed
reference residues in the selected chains. Atom-pair coverage uses the
possible atom pairs within the matched residues as its denominator. The
correspondence diagnostics list unmatched reference residues, and each
measurement lists omitted atoms. A fit needs at least three non-collinear
atom pairs. A measurement without matched atom pairs is unavailable. Binder
self-aligned results can remain available when the target fit fails.

A score exactly equal to its cutoff is a failure because both comparisons are
strict. A missing or non-finite score is unavailable, not a failure. Custom
cutoffs supplied to the batch runner are stored in the batch fingerprint and
manifest, CSV, and method-level summary. They are intentionally excluded from
prediction-run fingerprints so changing only a reporting cutoff can reuse the
existing structure scores.

To change only the DockQ cutoff when summarizing an existing CSV, pass
`--dockq-threshold` to `scripts/summarize_batch_metrics.py`. It recalculates
each prediction's pass fraction from its saved target–binder interface scores.
For older CSVs with counts alone, it uses those counts at the saved cutoff
(or the historical default of 0.23 if the CSV omits the cutoff). Success at a
different cutoff is unavailable for predictions without interface scores.

## Reruns

The fingerprint includes input hashes, complete annotations, scoring settings,
tool versions, and script hashes for both the chirality wrapper and its
`chirality_angles.py` companion. Matching successful results can be reused.
Changed fingerprints create separate runs. `--resume-incomplete` creates a
fresh directory under the original run's `attempts/` for each retry.
The CSV's `summary_path` identifies the actual attempt or reused result.
An empty discovery blocks the batch before any summaries are generated.

The batch fingerprint covers the complete annotation document, folding and
docking classification cutoffs, shared tool versions, and the ordered
planned/skipped prediction records. A changed batch creates a new directory
below `results/batches/`; historical batch manifests remain available.
