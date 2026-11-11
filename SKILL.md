---
name: protein-eval-skill
description: Audit and batch-evaluate Boltz protein/peptide complex predictions against annotated ground truth using OpenStructure, target-aligned and self-aligned binder RMSDs, and residue chirality. Use for evaluation directories containing recursive predicted structures and results; do not use for de novo structure prediction or unannotated automatic chain-role guessing.
---

# Protein Evaluation Skill

Evaluate only after biological roles, the exact reference assembly/prepared structure, sample matching, and every reference-to-prediction chain mapping are explicit. Never infer these from matching chain letters alone.

## Before evaluation

1. Inspect the repository's current evaluation scripts, especially `tools/evaluation.py` when present, and the bundled `scripts/chirality_angles_all_residues.py` and `scripts/chirality_angles.py`; do not assume evaluation-repository scripts are unchanged.
2. Verify `conda run -n ost212 ost --version` reports OpenStructure 2.12.0 and inspect `ost compare-structures --help`. Verify the bundled chirality CLI in `urop`.
   If either environment is absent, use `environment-urop.yml` and
   `environment-ost212.yml` from the skill root to prepare it after obtaining
   the user's permission to install packages.
3. Recursively enumerate only `.cif`, `.mmcif`, and `.pdb` prediction files beneath `predicted/`.
4. Inspect candidate reference and prediction structures with `scripts/inspect_structure.py`. Compare label and author chain IDs, ordered sequences, residue modifications, structure metadata, and available manifests. OpenStructure commonly exposes mmCIF label asym IDs; record the ID namespace explicitly.
5. Create or update the persistent annotations described in [references/annotations.md](references/annotations.md). Ask the user about every unresolved sample match, chain role, assembly, equivalent-chain choice, or chirality-chain scope before proceeding.
6. Run `scripts/batch_evaluate.py` without `--execute`. Do not execute until its audit is `READY` and the user has resolved any biological ambiguity.

Saved mappings are not authority by themselves. The audit must confirm that paths and chain IDs still exist; the single content-addressed batch manifest then binds all mappings to input hashes, settings, and tool versions.

## Execute

Run the orchestrator in `urop`, passing the repository's inspected chirality script:

```bash
conda run -n urop python <skill-dir>/scripts/batch_evaluate.py \
  <evaluation-dir> \
  --execute
```

The RMSD classifications default to strict `< 2 Å` for target self-aligned,
binder self-aligned, and target-aligned binder backbone RMSDs. DockQ defaults
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
replacement; its adjacent `chirality_angles.py` must be importable.

The orchestrator invokes OpenStructure in `ost212`. It supplies `-c` pairs in reference-structure chain order and adds `--min-pep-length 1` when an annotated binder chain has fewer than six observed residues. Review very-short-chain alignments and flag mismatches rather than trusting them automatically.

## Fixed metric definitions

The customized RMSDs are distinct from OpenStructure per-interface `lrmsd`:

- Establish residue correspondence from OpenStructure's validated, direction-preserving pairwise alignments.
- Use the aligned atom intersection. Fit predicted target to reference target once, across all annotated target chains together, with backbone `N`, `CA`, `C`, and `O` atoms and an unweighted Kabsch fit. Do not reject outliers.
- Require at least three non-collinear matched target atom pairs. Report the target backbone fit RMSD on those same atoms.
- Apply that transform to each annotated binder group without fitting the binder again.
- When sequence-identical binder chains are independent copies occupying distinct
  target pockets, annotate one single-chain binder group per copy. Measure each
  copy's target-aligned RMSD using the same target fit, and fit each copy
  independently for binder self-aligned RMSD. Preserve the validated one-to-one
  chain mapping. For example, in binder_target_5's 1PCG and 7LL8 use separate
  C and D groups; in 4GLN use separate A and D groups. A binder whose chains
  jointly form one functional unit remains one multichain binder group.
- Report the binder group count and aggregate scores consistently across methods;
  if the headline statistic is per structure, explicitly combine copy scores
  once per structure rather than counting a two-copy structure twice.
- Report both `target_aligned_binder_ca_rmsd_angstrom` and `target_aligned_binder_backbone_rmsd_angstrom` (`N`, `CA`, `C`, `O`), including atom/residue coverage and omissions.
- Never reverse residue correspondence to improve a score; a directionally reversed pose must remain penalized.
- Separately fit each annotated binder group to itself, independently for CA and
  backbone atoms, and report `binder_self_aligned_ca_rmsd_angstrom` and
  `binder_self_aligned_backbone_rmsd_angstrom`. For a multichain binder group,
  fit all member chains jointly. These self-aligned metrics measure folding or
  internal-assembly accuracy while deliberately discarding target-relative pose.
- Every fit requires at least three non-collinear matched atom pairs. Keep an
  unavailable self-aligned metric unavailable rather than substituting zero.

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

For docking, use only native interfaces between an annotated target chain and
an annotated binder chain. Retain the unfiltered OpenStructure output as a
diagnostic. Report per-interface DockQ, the equal mean of eligible target-binder
DockQs within a prediction, and their native-contact-weighted DockQ-wave.
Apply the strict DockQ cutoff to each target-binder interface first; the
prediction's success value is the fraction of eligible interfaces that pass.
Average those fractions across available predictions. Filter iRMSD to those
same target-binder interfaces without changing its existing dataset-level
median across available interface values.

Compute target-binder iLDDT separately for each target-chain/binder-chain
pair with native inter-chain distance tests. The macro value is the equal mean
of eligible pair scores within a prediction; the micro value pools conserved
and expected OpenStructure distance tests across those pairs within that
prediction. Use macro for the docking headline and retain micro, pair scores,
and the original global inter-chain iLDDT as diagnostics. Pairs with no native
distance tests are unavailable, not zero. Read [references/outputs.md](references/outputs.md)
for exact fields and aggregation units.

Chirality uses the bundled `scripts/chirality_angles_all_residues.py` with detailed output. The primary accuracy is pooled over every eligible D- and L-amino-acid residue in the explicitly annotated chirality chains. Preserve separate DAA and LAA summaries, per-residue classifications, correct count, eligible count, and skipped/unassessable count. Unavailable values stay unavailable, never zero-filled.

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
  mean lDDT, and chirality violations.
- **Docking:** mean target-binder DockQ, complex-balanced interface DockQ
  success, target-aligned binder backbone RMSD and threshold success,
  target-binder iRMSD, and mean target-binder macro iLDDT.
- **Supplemental:** original global lDDT, target-binder DockQ-wave, and
  target-binder micro iLDDT.

Use the deterministic aggregation generated by `scripts/summarize_batch_metrics.py`:
RMSDs are medians, bounded scores are arithmetic means, and chirality violation
is the residue-pooled count and percentage of inverted plus below-threshold
residues among assessable residues. Report unassessable chirality residues
separately. Exclude explicitly skipped predictions and show the number of
predictions included for each method. Link `results/method_summary.md`,
`results/method_summary.json`, and `results/batch_summary.csv` in the reply.
Threshold success rates use only available values in their denominators. Target
self-aligned RMSD is counted once per prediction; both binder RMSDs are counted
once per annotated binder group. DockQ success is an equal mean of the
per-prediction target-binder interface success fractions. State all four cutoffs
beside the reported rates. The dataset-level RMSD medians, bounded-score means,
and residue-pooled chirality aggregation remain unchanged.

## Boundaries

- Polymer peptide/protein binders always use `ost compare-structures`, never `compare-ligand-structures`.
- If a non-polymer ligand is present, ask which ligand instances to evaluate before adding a separate `compare-ligand-structures` branch. Keep lDDT-PLI and pocket-aligned ligand RMSD separate from polymer metrics. The bundled runner intentionally does not execute this branch.
- Keep all threshold results independent; do not silently combine them into a
  joint folding- or docking-success label.
- Do not call OpenStructure per-interface `lrmsd` an RMSD for an entire multichain binder against an entire multichain target.
- Do not score interfaces involving context chains. Treat any such interface reported by OpenStructure as a warning requiring review.
- Do not overwrite successful results. Identical fingerprints reuse their run; changed inputs, annotations, settings, scripts, or versions receive a new run directory.
