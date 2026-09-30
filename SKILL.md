---
name: protein-eval-skill
description: Audit and batch-evaluate Boltz protein/peptide complex predictions against annotated ground truth using OpenStructure, target-aligned and self-aligned binder RMSDs, and residue chirality. Use for evaluation directories containing recursive predicted structures and results; do not use for de novo structure prediction or unannotated automatic chain-role guessing.
---

# Protein Evaluation Skill

Evaluate only after biological roles, the exact reference assembly/prepared structure, sample matching, and every reference-to-prediction chain mapping are explicit. Never infer these from matching chain letters alone.

## Before evaluation

1. Inspect the repository's current evaluation scripts, especially `tools/evaluation.py` and `tools/chirality_angles_all_residues.py`; do not assume they are unchanged.
2. Verify `conda run -n ost212 ost --version` reports OpenStructure 2.12.0 and inspect `ost compare-structures --help`. Verify the chirality CLI in `urop`.
3. Recursively enumerate only `.cif`, `.mmcif`, and `.pdb` prediction files beneath `predicted/`.
4. Inspect candidate reference and prediction structures with `scripts/inspect_structure.py`. Compare label and author chain IDs, ordered sequences, residue modifications, structure metadata, and available manifests. OpenStructure commonly exposes mmCIF label asym IDs; record the ID namespace explicitly.
5. Create or update the persistent annotations described in [references/annotations.md](references/annotations.md). Ask the user about every unresolved sample match, chain role, assembly, equivalent-chain choice, or chirality-chain scope before proceeding.
6. Run `scripts/batch_evaluate.py` without `--execute`. Do not execute until its audit is `READY` and the user has resolved any biological ambiguity.

Saved mappings are not authority by themselves. The audit must confirm that paths and chain IDs still exist; the content-addressed run manifest then binds the mapping to input hashes, settings, and tool versions.

## Execute

Run the orchestrator in `urop`, passing the repository's inspected chirality script:

```bash
conda run -n urop python <skill-dir>/scripts/batch_evaluate.py \
  <evaluation-dir> \
  --chirality-script <repo-root>/tools/chirality_angles_all_residues.py \
  --execute
```

The orchestrator invokes OpenStructure in `ost212`. It supplies `-c` pairs in reference-structure chain order and adds `--min-pep-length 1` when an annotated binder chain has fewer than six observed residues. Review very-short-chain alignments and flag mismatches rather than trusting them automatically.

## Fixed metric definitions

The customized RMSDs are distinct from OpenStructure per-interface `lrmsd`:

- Establish residue correspondence from OpenStructure's validated, direction-preserving pairwise alignments.
- Use the aligned atom intersection. Fit predicted target to reference target once, across all annotated target chains together, with backbone `N`, `CA`, `C`, and `O` atoms and an unweighted Kabsch fit. Do not reject outliers.
- Require at least three non-collinear matched target atom pairs. Report the target backbone fit RMSD on those same atoms.
- Apply that transform to each annotated multichain binder group without fitting the binder again.
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

Chirality uses the inspected `chirality_angles_all_residues.py` with detailed output. The primary accuracy is pooled over every eligible D- and L-amino-acid residue in the explicitly annotated chirality chains. Preserve separate DAA and LAA summaries, per-residue classifications, correct count, eligible count, and skipped/unassessable count. Unavailable values stay unavailable, never zero-filled.

Read [references/outputs.md](references/outputs.md) when interpreting, validating, or reporting results.

## Boundaries

- Polymer peptide/protein binders always use `ost compare-structures`, never `compare-ligand-structures`.
- If a non-polymer ligand is present, ask which ligand instances to evaluate before adding a separate `compare-ligand-structures` branch. Keep lDDT-PLI and pocket-aligned ligand RMSD separate from polymer metrics. The bundled runner intentionally does not execute this branch.
- Do not calculate dataset-level docking success, thresholded accuracy, or success proportions.
- Do not call OpenStructure per-interface `lrmsd` an RMSD for an entire multichain binder against an entire multichain target.
- Do not score interfaces involving context chains. Treat any such interface reported by OpenStructure as a warning requiring review.
- Do not overwrite successful results. Identical fingerprints reuse their run; changed inputs, annotations, settings, scripts, or versions receive a new run directory.
