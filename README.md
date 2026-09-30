# protein-eval-skill

A Codex skill for audited batch evaluation of Boltz protein/peptide complex
predictions against annotated ground truth. It combines OpenStructure scores,
target-aligned and self-aligned binder RMSDs, and residue chirality metrics.

## Install from GitHub

After replacing `OWNER` with the GitHub owner, install the standalone repository
with Codex's skill installer:

```bash
python ~/.codex/skills/.system/skill-installer/scripts/install-skill-from-github.py \
  --repo OWNER/protein-eval-skill \
  --path . \
  --name protein-eval-skill
```

The skill becomes available as `$protein-eval-skill` on the next Codex turn.

## Runtime requirements

- Conda environment `urop` with Python, NumPy, and Gemmi.
- Conda environment `ost212` with OpenStructure 2.12.0.
- The evaluation repository's `tools/chirality_angles_all_residues.py` and
  adjacent `tools/chirality_angles.py`.
- An evaluation directory with persistent biological-role and chain-mapping
  annotations as described in `references/annotations.md`.

The environment names are currently fixed in the runner. The chirality script
is supplied at execution time with `--chirality-script`.

## Validate

```bash
conda run -n urop python \
  ~/.codex/skills/.system/skill-creator/scripts/quick_validate.py .
```

## Use

```text
Use $protein-eval-skill to audit and evaluate <evaluation-directory>.
```

The skill audits before execution and asks for unresolved biological roles,
assemblies, or chain mappings rather than guessing them.
