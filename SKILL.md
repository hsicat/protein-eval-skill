---
name: protein-eval-skill
description: Audit and batch-evaluate Boltz protein/peptide complex predictions against annotated ground truth using OpenStructure, target-aligned and self-aligned binder RMSDs, and residue chirality. Use for evaluation directories containing recursive predicted structures and results; do not use for de novo structure prediction or unannotated automatic chain-role guessing.
---

# Protein Evaluation Skill

Evaluate only after biological roles, the exact reference assembly/prepared structure, sample matching, and every reference-to-prediction chain mapping are explicit. Never infer these from matching chain letters alone.

## Before evaluation

1. Locate the evaluation directory, its references, and any existing annotations. Check the bundled scoring scripts if reviewing an unfamiliar installation or a proposed override.
2. Verify `conda run -n ost212 ost --version` reports OpenStructure 2.12.0 and inspect `ost compare-structures --help`. Verify the bundled chirality CLI in `urop`.
   If either environment is absent, use `environment-urop.yml` and
   `environment-ost212.yml` from the skill root to prepare it after obtaining
   the user's permission to install packages.
3. Recursively enumerate only `.cif`, `.mmcif`, and `.pdb` prediction files beneath `predicted/`.
4. Inspect candidate reference and prediction structures with `scripts/inspect_structure.py`. Compare label and author chain IDs, ordered sequences, residue modifications, structure metadata, and available manifests. OpenStructure commonly exposes mmCIF label asym IDs; record the ID namespace explicitly.
5. Create or update the persistent annotations described in [references/annotations.md](references/annotations.md). Ask the user about every unresolved sample match, chain role, assembly, equivalent-chain choice, or chirality-chain scope before proceeding.
6. Run `scripts/batch_evaluate.py` without `--execute`. Do not execute until its audit is `READY` and the user has resolved any biological ambiguity.

Saved mappings are not authority by themselves. The audit must confirm that paths and chain IDs still exist; the single content-addressed batch manifest then binds all mappings to input hashes, settings, and tool versions.

The batch scores the one-to-one mapping in each prediction's annotation. If the
benchmark calls for the highest mean DockQ among sequence-equivalent mappings,
compare the role-preserving candidates in a separate pre-batch step. Record
the selected mapping in the annotation and use it for all metrics.

Use the evaluation directory's own `annotations/samples.json` by default,
or a user-selected file passed with `--annotations`. The bundled
[RCG eval_v01 mappings](ref_mapping_1127/README.md) are optional reference
material. Use them for scoring only when the user opts in and the input
checksums match the corresponding archived benchmark files.

## Execute

Run the orchestrator in `urop` with the bundled chirality scorer:

```bash
conda run -n urop python <skill-dir>/scripts/batch_evaluate.py \
  <evaluation-dir> \
  --execute
```

The RMSD classifications default to strict `< 2 Å` for target self-aligned
backbone RMSD, binder self-aligned backbone RMSD, and L-RMSD. DockQ defaults
to strict `> 0.23`. Override them for a batch with
`--target-self-aligned-rmsd-threshold-angstrom <value>`,
`--binder-self-aligned-rmsd-threshold-angstrom <value>`,
`--binder-rmsd-threshold-angstrom <value>`, and
`--dockq-threshold <value>`.
The thresholds are fingerprinted batch/report settings and therefore remain
traceable without forcing the underlying structure scores to be recomputed.

For a single-method evaluation directory whose `predicted/` subdirectories
are seed names (such as AF3 baseline or AF3 steering), pass
`--method-name af3_baseline` or `--method-name af3_steering` so the
method summary combines all seeds under the actual method name.

The runner uses the chirality scripts bundled in `<skill-dir>/scripts/` by
default. Pass `--chirality-script <path>` only to test an explicitly inspected
replacement. It must support `--chain-id-type auth|label` and have an
importable `chirality_angles.py` alongside it. For a whole-example
chirality result, pass `--chirality-scope whole-example` to pool all scored
target and binder chains. The default uses the annotated `chirality_chains`.

The orchestrator invokes OpenStructure in `ost212` with only the annotated
target and binder polymer chains. It supplies their `-c` pairs in
reference-structure chain order and adds `--min-pep-length 1` when an annotated
binder chain has fewer than six observed residues. Review very-short-chain
alignments and flag mismatches.

## Fixed metric definitions

Use **L-RMSD** in figures, tables, axes, panel titles, captions, and threshold-success
labels for `target_aligned_binder_backbone_rmsd_angstrom`. Treat
“target-aligned binder backbone RMSD” and “Binder pose RMSD” as legacy aliases.
Define L-RMSD as binder backbone (`N`, `CA`, `C`, `O`) RMSD after jointly fitting
all annotated target chains on their backbone atoms, without refitting the
binder. Preserve the calculation, mappings, grouping, aggregation, thresholds,
CLI options, and stored JSON/CSV field names. Keep target/binder self-aligned
RMSDs and the CA-only diagnostic under their own names.

The customized RMSDs are distinct from OpenStructure per-interface `lrmsd`:

- Establish residue correspondence from OpenStructure's validated, direction-preserving pairwise alignments.
- Use the aligned atom intersection. Fit predicted target to reference target once, across all annotated target chains together, with backbone `N`, `CA`, `C`, and `O` atoms and an unweighted Kabsch fit. Do not reject outliers.
- Require at least three non-collinear matched target atom pairs. Report the target backbone fit RMSD on those same atoms.
- Apply that transform to each annotated binder group without fitting the binder again.
- When sequence-identical binder chains are independent copies occupying distinct
  target pockets, annotate one single-chain binder group per copy. Measure each
  copy's target-aligned RMSD using the same target fit, and fit each copy
  independently for binder self-aligned RMSD. Preserve the validated one-to-one
  chain mapping. A binder whose chains
  jointly form one functional unit remains one multichain binder group.
- Report the binder group count and aggregate scores consistently across methods.
  For method-level binder RMSD headlines, average available copy scores once per
  prediction, then take the median across predictions. Retain the binder-group
  pooled median only as a labeled diagnostic.
- Report both `target_aligned_binder_ca_rmsd_angstrom` and `target_aligned_binder_backbone_rmsd_angstrom` (`N`, `CA`, `C`, `O`). Report matched reference-residue coverage separately from atom-pair coverage within those matched residues.
- Never reverse residue correspondence to improve a score; a directionally reversed pose must remain penalized.
- Separately fit each annotated binder group to itself, independently for CA and
  backbone atoms, and report `binder_self_aligned_ca_rmsd_angstrom` and
  `binder_self_aligned_backbone_rmsd_angstrom`. For a multichain binder group,
  fit all member chains jointly. These self-aligned metrics measure folding or
  internal-assembly accuracy while deliberately discarding target-relative pose.
- Every fit requires at least three non-collinear matched atom pairs. Keep an
  unavailable self-aligned metric unavailable rather than substituting zero.
  If the target fit fails, report target RMSD and L-RMSD as unavailable and
  retain binder self-aligned RMSDs wherever the binder correspondence supports a fit.

Run OpenStructure with global all-atom lDDT, local lDDT, interface lDDT
(`--ilddt`), and DockQ. OpenStructure `lrmsd` is a per-interface DockQ
component: it treats the larger chain as receptor, fits that receptor using
mapped backbone atoms, and reports RMSD on the smaller chain. It may numerically
match a target-aligned binder-backbone RMSD for a simple two-chain complex, but
it is not a substitute for the explicitly grouped custom metric.

For folding, additionally compute all-atom lDDT on all annotated target
chains alone and on each annotated binder group alone. The target-only score
includes target-target contacts; a multichain binder group's score includes
its internal contacts. Average eligible binder-group lDDTs once within each
prediction for the method-level folding headline. Keep the original full
target-plus-binder lDDT as a supplemental diagnostic.

Compute all-atom lDDT and iLDDT with OpenStructure stereochemistry checks both
enabled and disabled. Use the same scoring subsets, chain mappings, and
aggregation rules for each mode. Store both versions, with explicit
`_no_stereochecks` fields or files for the disabled mode. OpenStructure's
standard input cleanup applies in both modes.

For docking, use only native interfaces between an annotated target chain and
an annotated binder chain. Retain the unfiltered OpenStructure output as a
diagnostic. Report per-interface DockQ, the equal mean of eligible target-binder
DockQs within a prediction, and their native-contact-weighted DockQ-wave.
Apply the strict DockQ cutoff to each target-binder interface first; the
prediction's success value is the fraction of eligible interfaces that pass.
Average those fractions across available predictions. For iRMSD, average
available target-binder interface values once per prediction, then take the
median across predictions. Retain the interface-pooled median as a diagnostic.

Compute target-binder iLDDT separately for each target-chain/binder-chain
pair with native inter-chain distance tests. The macro value is the equal mean
of eligible pair scores within a prediction; the micro value pools conserved
and expected OpenStructure distance tests across those pairs within that
prediction. Use macro for the docking headline and retain micro, pair scores,
and the original global inter-chain iLDDT as diagnostics. Pairs with no native
distance tests are unavailable, not zero. Read [references/outputs.md](references/outputs.md)
for exact fields and aggregation units.

Chirality passes when the oriented Cα improper angle is at least 15°.
Compute accuracy as passed / assessable residues for All (D+L), D-only, and
L-only. Pool counts across the selected chains for an example and across the
included complex–seed predictions for a method. Record chain-wise results,
per-residue classifications, and unassessable counts. Use NA when there are
zero assessable residues. Interpret chain IDs in the model's declared namespace.

For every available target-binder interface DockQ, independently classify it
with a strict `>` comparison. Use strict `<` comparisons for the target
self-aligned backbone RMSD, each annotated binder group's self-aligned backbone
RMSD, and each binder group's target-aligned backbone RMSD. A value exactly at
the cutoff does not pass. Keep unavailable scores unavailable rather than
counting them as failures.

Read [references/outputs.md](references/outputs.md) when interpreting, validating, or reporting results.

## Report the completed batch

After every executed batch, read `results/method_summary.md` and reproduce all
three Markdown tables in the final reply. Do not report only a prose conclusion.
The required charts are:

- **Folding:** target self-aligned backbone RMSD and threshold success, binder
  self-aligned backbone RMSD and threshold success, target-only and binder-only
  mean lDDT, and pooled All (D+L), D-only, and L-only chirality accuracy.
- **Docking:** mean target-binder DockQ, prediction-balanced interface DockQ
  success, L-RMSD and L-RMSD success,
  target-binder iRMSD, and mean target-binder macro iLDDT.
- **Supplemental:** original global lDDT, target-binder DockQ-wave, and
  target-binder micro iLDDT.

Show adjacent, explicitly labeled checks-on and checks-off columns for
target-only, binder-only and global lDDT, and macro and micro iLDDT.
Use NA for the checks-off columns when reading older results without those scores.

Use the deterministic aggregation generated by `scripts/summarize_batch_metrics.py`:
target RMSD is a median across available predictions. Binder RMSDs and
target-binder iRMSD first average eligible binder groups or interfaces within
each prediction, then take a median across predictions. Classify each available
binder group with the strict RMSD threshold, average the pass fraction within
each prediction, then average those fractions across predictions. Bounded
scores are arithmetic means across predictions. Each included complex–seed
prediction contributes one structural value per available metric in the
all-seed summary, so complexes with more seeds contribute more observations.
The batch produces the all-seed tables. For a separate metric-specific oracle
analysis, use the per-prediction results to select one prediction per complex
for each metric, then summarize the selected values. Chirality accuracy pools
residue counts across the included predictions. Report unassessable chirality
residues separately. Exclude explicitly skipped predictions and show the number
of predictions included for each method. Link `results/method_summary.md`,
`results/method_summary.json`, and `results/batch_summary.csv` in the reply.
Threshold success rates use only available values in their denominators. Target
self-aligned RMSD is counted once per prediction. Both binder RMSD success rates
are means of per-prediction binder-group pass fractions. DockQ success is an
equal mean of the per-prediction target-binder interface success fractions.
State all four cutoffs beside the reported rates. The older group-pooled RMSD
and success results, and interface-pooled iRMSD, remain in
`method_summary.json` under `supplemental_pooled` for comparison.

For manuscript figures or tables, use `protein-eval-graph` and preserve its
seven selected primary metric categories: target and binder self-aligned
backbone RMSD, global lDDT, D-Cα chirality correctness, L-RMSD,
DockQ, and interface lDDT. Interface RMSD is optional. Verify
that plotted lDDT, DockQ, and chirality fields match the figure's stated
definitions; this skill's additional diagnostic metrics do not add primary
figure panels automatically.
When invoking `protein-eval-graph` for the standard seven-panel comparison,
explicitly select `--metrics chirality,target,binder,lddt,placed,ilddt,dockq`;
its general default also includes optional iRMSD.
The graph skill uses the batch summary's within-prediction reductions and
target–binder interface fields by default. Global lDDT comes from the
Supplemental table, and the D-Cα panel uses D-only pooled accuracy. Pass
`--stereochecks off` to plot the saved checks-off lDDT and iLDDT scores.

## Boundaries

- Polymer peptide/protein binders always use `ost compare-structures`, never `compare-ligand-structures`.
- If a non-polymer ligand is present, ask which ligand instances to evaluate before adding a separate `compare-ligand-structures` branch. Keep lDDT-PLI and pocket-aligned ligand RMSD separate from polymer metrics. The bundled runner intentionally does not execute this branch.
- Keep all threshold results independent; do not silently combine them into a
  joint folding- or docking-success label.
- Do not call OpenStructure per-interface `lrmsd` an RMSD for an entire multichain binder against an entire multichain target.
- Do not score interfaces involving context chains. Treat any such interface reported by OpenStructure as a warning requiring review.
- Reuse a successful result only when its full fingerprint matches. A changed
  fingerprint gets a new run directory. Pass `--resume-incomplete` to retry an
  incomplete fingerprint in a fresh attempt directory. The batch reports only
  the results produced or explicitly reused during that execution.
