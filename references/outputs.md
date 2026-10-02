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
- `custom_rmsd_config.json` and `binder_rmsd.json`.
- `chirality/`: the chirality tool's readable report, simple CSV, detailed CSV, and raw summary JSON.
- `summary.json`: parsed per-prediction summary retaining local lDDT, mappings, alignments, residue inconsistencies, cleanup/stereochemical diagnostics, per-interface identities and scores, DockQ averages/full averages, custom RMSDs, chirality summaries, statuses, warnings, and unavailable reasons.

`results/batch_summary.csv` has one row per prediction/binder for the fingerprints
selected by the latest batch execution, plus explicit skipped rows. Historical
content-addressed runs remain on disk. Per-interface OpenStructure results remain
present as JSON in the CSV and as structured data in each summary; never collapse
them to only one DockQ value.
The CSV exposes pooled all-residue chirality as the primary result and also retains
separate DAA and LAA counts and accuracies.

## Metric names

- `openstructure.lddt`: global all-atom lDDT.
- `openstructure.ilddt`: global all-atom lDDT restricted to inter-chain contacts.
- `openstructure.local_lddt`: reference-residue keyed all-atom local lDDT.
- `openstructure.per_interface[*].lrmsd`: OpenStructure's DockQ ligand RMSD for
  one chain-pair interface. The larger chain is treated as receptor and fitted;
  RMSD is measured on the smaller chain using mapped backbone atoms. This is not
  the custom grouped-binder metric.
- `target_backbone_fit_rmsd_angstrom`: target RMSD after the joint all-target N/CA/C/O fit.
- `target_aligned_binder_ca_rmsd_angstrom`: mapped binder CA displacement after the target fit, without binder refitting.
- `target_aligned_binder_backbone_rmsd_angstrom`: mapped binder N/CA/C/O displacement after the target fit, without binder refitting.
- `binder_self_aligned_ca_rmsd_angstrom`: CA RMSD after independently fitting
  that binder group on its CA atoms.
- `binder_self_aligned_backbone_rmsd_angstrom`: N/CA/C/O RMSD after independently
  fitting that binder group on its N/CA/C/O atoms.

Missing residues or atoms are excluded only through the recorded aligned intersection. Always report used/possible atom pairs, coverage, and the omission list alongside a custom RMSD. If the target fit has fewer than three non-collinear atom pairs or a binder measurement has no atom pairs, that metric is unavailable with a reason.

## Reruns

The fingerprint includes reference and prediction file hashes, complete sample/prediction annotations, settings, inspected script hashes, and tool versions. An existing successful fingerprint is reused. A changed fingerprint creates another run and preserves the earlier result. An incomplete fingerprint is not resumed unless the operator reviews it and passes `--resume-incomplete`.

The batch fingerprint covers the complete annotation document, shared tool
versions, and the ordered planned/skipped prediction records. A changed batch
creates a new directory below `results/batches/`; historical batch manifests
remain available.
