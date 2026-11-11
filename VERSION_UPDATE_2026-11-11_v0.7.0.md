# v0.7.0 dataset aggregation update — 2026-11-11

## Scope

This release changes **method/dataset aggregation** for two binder backbone
RMSDs, their strict threshold success rates, and target–binder iRMSD. It does
not change the v0.6.0 per-complex structure scoring, atom alignment, chain
roles, thresholds, or OpenStructure calculations. In particular, chirality
remains residue-pooled across all assessable D- and L-amino-acid residues.

## Changed method-level headline values

| Metric | v0.6.0 | v0.7.0 |
|---|---|---|
| Binder self-aligned backbone RMSD | Median across available binder-group RMSDs | Mean available binder-group RMSDs once per prediction; median across available prediction means |
| Target-aligned binder backbone RMSD | Median across available binder-group RMSDs | Mean available binder-group RMSDs once per prediction; median across available prediction means |
| Both binder RMSD success rates | Fraction of available binder groups below their respective strict `< 2 Å` cutoffs | Fraction passing within each prediction, then arithmetic mean of available prediction fractions |
| Target–binder iRMSD | Median across available target–binder interface values | Mean available interface iRMSDs once per prediction; median across available prediction means |

The old group-pooled and interface-pooled results remain in
`method_summary.json.methods[*].supplemental_pooled`, with their original units
and denominators. The primary method summary uses one available value per
prediction for each changed metric. A prediction with no available binder
group or interface value is excluded from that metric only; missing values
are never substituted with zero. A prediction with multiple eligible groups
or interfaces is no longer overrepresented in those method headlines.

The RMSD cutoff is applied to each available binder-group RMSD before
averaging pass fractions. The complex-level success fraction is **not** the
Boolean result of comparing its mean RMSD with the threshold, and it does not
mean every binder group passed. The method-level success percentage is the
mean of these fractions. Counts of available complexes and successful/available
groups accompany the percentage; the pooled group ratio is diagnostic and
may differ from that percentage.

## Unchanged method-level rules

- Target self-aligned backbone RMSD: one value per prediction, method median;
  strict threshold success is per prediction.
- Target-only, binder-only, and original global lDDT: arithmetic means of
  available prediction values. Binder-only lDDT still averages eligible binder
  groups inside each prediction.
- Target–binder DockQ and DockQ-wave: arithmetic means of available prediction
  values. DockQ success remains the mean of per-prediction target–binder
  interface pass fractions.
- Target–binder iLDDT macro and micro: arithmetic means of available prediction
  values. Micro pools distance tests only within each prediction.
- **Chirality:** correct, assessable, and unassessable residue counts retain
  the existing residue-pooled rule. DAA and LAA outputs retain their existing
  per-residue definitions. No complex-equal chirality headline is added.

## Figures and compatibility

The `protein-eval-graph` skill and its selected seven primary manuscript
metric categories are unchanged: target and binder self-aligned backbone RMSD,
global lDDT, D-Cα chirality correctness, target-aligned binder backbone RMSD,
DockQ, and interface lDDT. iRMSD remains an optional eighth panel. Paper
figure fields must still be checked against the requested metric definitions;
the graph skill currently has its own score-column and chirality conventions.

`method_summary.json` moves from schema version 4 to 5. Its existing primary
RMSD and success field names now contain the complex-balanced values;
`supplemental_pooled` supplies the prior comparison values. The Markdown
headings explicitly say “median of complex means” and identify binder success
as a mean of per-complex group pass fractions. Existing per-run `summary.json`
and `batch_summary.csv` data already contain the necessary binder-group and
per-interface values; a v0.6.0 batch can be summarized again without
re-running structure scoring. Existing scored runs are not rewritten.

This release still treats each **prediction** as one complex observation. If
one biological sample has more seeds than another, it still contributes more
prediction observations. Equal weighting of unique biological samples would
be a separate future policy change.
