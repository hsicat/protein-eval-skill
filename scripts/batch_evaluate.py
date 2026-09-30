#!/usr/bin/env python3
"""Audit or execute an annotated batch of Boltz structure evaluations.

Run this script in the ``urop`` conda environment. OpenStructure is invoked in
the separate ``ost212`` environment. The script never infers biological roles,
sample matches, assemblies, or chain mappings.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
from typing import Any

import gemmi

from target_aligned_binder_rmsd import chain_residues


STRUCTURE_SUFFIXES = {".cif", ".mmcif", ".pdb"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evaluation_dir", type=Path)
    parser.add_argument(
        "--annotations",
        type=Path,
        help="Default: <evaluation_dir>/annotations/samples.json",
    )
    parser.add_argument(
        "--chirality-script",
        type=Path,
        help="Path to tools/chirality_angles_all_residues.py (required with --execute)",
    )
    parser.add_argument("--execute", action="store_true", help="Run tools after a clean audit")
    parser.add_argument(
        "--audit-output", type=Path, help="Optionally save the audit JSON; it is always printed"
    )
    parser.add_argument(
        "--resume-incomplete",
        action="store_true",
        help="Allow tools to be rerun in an existing non-successful content-addressed run",
    )
    return parser.parse_args()


def load_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def discover_predictions(evaluation_dir: Path) -> list[str]:
    predicted = evaluation_dir / "predicted"
    if not predicted.is_dir():
        return []
    return sorted(
        path.relative_to(evaluation_dir).as_posix()
        for path in predicted.rglob("*")
        if path.is_file() and path.suffix.lower() in STRUCTURE_SUFFIXES
    )


def resolved(evaluation_dir: Path, relative: str) -> Path:
    path = Path(relative)
    return path if path.is_absolute() else (evaluation_dir / path).resolve()


def chain_inventory(path: Path, chain_id_type: str) -> dict[str, int]:
    return {chain: len(residues) for chain, residues in chain_residues(path, chain_id_type).items()}


def write_polymer_subset(
    input_path: Path, output_path: Path, chain_id_type: str, included_chains: set[str]
) -> dict[str, Any]:
    """Materialize exactly the selected polymer chains while preserving coordinates/IDs."""
    structure = gemmi.read_structure(str(input_path)).clone()
    if len(structure) != 1:
        raise ValueError(f"Expected exactly one model in {input_path}; found {len(structure)}")
    model = structure[0]
    removed_polymer_residues: list[dict[str, str]] = []
    removed_nonpolymer_count = 0
    for chain_index in range(len(model) - 1, -1, -1):
        chain = model[chain_index]
        for residue_index in range(len(chain) - 1, -1, -1):
            residue = chain[residue_index]
            info = gemmi.find_tabulated_residue(residue.name)
            chain_id = chain.name if chain_id_type == "auth" else residue.subchain
            if not info.is_amino_acid() or chain_id not in included_chains:
                if info.is_amino_acid():
                    removed_polymer_residues.append(
                        {
                            "auth_chain": chain.name,
                            "label_chain": residue.subchain,
                            "residue": str(residue.seqid),
                            "name": residue.name,
                        }
                    )
                else:
                    removed_nonpolymer_count += 1
                del chain[residue_index]
        if not len(chain):
            model.remove_chain(chain.name)
    present = set(chain_inventory_from_structure(structure, chain_id_type))
    if present != included_chains:
        raise ValueError(
            f"Subset chain mismatch for {input_path}: requested={sorted(included_chains)}, "
            f"written={sorted(present)}"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    structure.make_mmcif_document().write_file(str(output_path))
    return {
        "input": str(input_path),
        "output": str(output_path),
        "chain_id_type": chain_id_type,
        "included_chains": sorted(included_chains),
        "removed_polymer_residues": removed_polymer_residues,
        "removed_nonpolymer_residue_count": removed_nonpolymer_count,
    }


def chain_inventory_from_structure(structure: Any, chain_id_type: str) -> dict[str, int]:
    grouped: dict[str, int] = {}
    for chain in structure[0]:
        for residue in chain:
            if not gemmi.find_tabulated_residue(residue.name).is_amino_acid():
                continue
            chain_id = chain.name if chain_id_type == "auth" else residue.subchain
            grouped[chain_id] = grouped.get(chain_id, 0) + 1
    return grouped


def prediction_records(annotations: dict[str, Any]) -> dict[str, tuple[str, dict[str, Any], dict[str, Any]]]:
    records: dict[str, tuple[str, dict[str, Any], dict[str, Any]]] = {}
    for sample_id, sample in annotations.get("samples", {}).items():
        for prediction_path, prediction in sample.get("predictions", {}).items():
            if prediction_path in records:
                raise ValueError(f"Prediction is annotated more than once: {prediction_path}")
            records[prediction_path] = (sample_id, sample, prediction)
    return records


def audit(evaluation_dir: Path, annotation_path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    report: dict[str, Any] = {
        "schema_version": 1,
        "evaluation_dir": str(evaluation_dir),
        "annotations": str(annotation_path),
        "discovered_predictions": discover_predictions(evaluation_dir),
        "errors": [],
        "warnings": [],
        "predictions": [],
    }
    if not annotation_path.is_file():
        report["errors"].append(f"Annotation file is missing: {annotation_path}")
        report["status"] = "BLOCKED"
        return report, {}
    try:
        annotations = load_json(annotation_path)
        if annotations.get("schema_version") != 1:
            raise ValueError("annotations schema_version must be 1")
        records = prediction_records(annotations)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        report["errors"].append(str(exc))
        report["status"] = "BLOCKED"
        return report, {}

    discovered = set(report["discovered_predictions"])
    annotated = set(records)
    for path in sorted(discovered - annotated):
        report["errors"].append(
            f"Unannotated prediction; sample matching and chain roles must be resolved: {path}"
        )
    for path in sorted(annotated - discovered):
        report["errors"].append(f"Annotated prediction path is missing: {path}")

    for prediction_path in sorted(discovered & annotated):
        sample_id, sample, prediction = records[prediction_path]
        item: dict[str, Any] = {
            "sample_id": sample_id,
            "prediction": prediction_path,
            "errors": [],
            "warnings": [],
        }
        if sample.get("skip") or prediction.get("skip"):
            skip_reason = prediction.get("skip_reason") or sample.get("skip_reason")
            if not isinstance(skip_reason, str) or not skip_reason.strip():
                item["errors"].append("Skipped predictions require a nonempty skip_reason")
            else:
                item.update({"status": "SKIPPED", "skip_reason": skip_reason.strip()})
            report["predictions"].append(item)
            report["errors"].extend(
                f"{prediction_path}: {message}" for message in item["errors"]
            )
            continue
        try:
            reference_text = sample["reference"]
            reference = resolved(evaluation_dir, reference_text)
            model = resolved(evaluation_dir, prediction_path)
            if not reference.is_file():
                raise ValueError(f"Reference file is missing: {reference}")
            assembly = sample.get("reference_assembly")
            if not isinstance(assembly, dict) or assembly.get("selection") not in {
                "biological_assembly",
                "asymmetric_unit",
                "prepared_subset",
            }:
                raise ValueError(
                    "reference_assembly.selection must explicitly be biological_assembly, "
                    "asymmetric_unit, or prepared_subset"
                )
            reference_id_type = sample["reference_chain_id_type"]
            model_id_type = prediction["model_chain_id_type"]
            if reference_id_type not in {"auth", "label"} or model_id_type not in {"auth", "label"}:
                raise ValueError("chain ID types must be auth or label")
            ref_inventory = chain_inventory(reference, reference_id_type)
            mdl_inventory = chain_inventory(model, model_id_type)
            target_chains = sample["target_chains"]
            binders = sample["binders"]
            if not target_chains:
                raise ValueError("target_chains must not be empty")
            if not binders or any(not binder.get("id") or not binder.get("chains") for binder in binders):
                raise ValueError("binders must contain at least one named, nonempty chain group")
            binder_ids = [binder["id"] for binder in binders]
            if len(binder_ids) != len(set(binder_ids)):
                raise ValueError("binder IDs must be unique within a sample")
            context = sample.get("context_chains", [])
            roles = target_chains + [c for binder in binders for c in binder["chains"]] + context
            if len(roles) != len(set(roles)):
                raise ValueError("A reference chain occurs in more than one biological role")
            unknown_roles = sorted(set(roles) - set(ref_inventory))
            if unknown_roles:
                raise ValueError(f"Annotated reference chains are absent: {unknown_roles}")
            if "chirality_chains" not in sample:
                raise ValueError(
                    "chirality_chains must explicitly identify which reference chains are scored"
                )
            unknown_chirality = sorted(set(sample["chirality_chains"]) - set(ref_inventory))
            if unknown_chirality:
                raise ValueError(f"Chirality reference chains are absent: {unknown_chirality}")
            mapping = prediction["chain_mapping"]
            scored_chains = target_chains + [c for binder in binders for c in binder["chains"]]
            missing_mapping = sorted(set(scored_chains + sample["chirality_chains"]) - set(mapping))
            if missing_mapping:
                raise ValueError(f"Required reference chains lack mappings: {missing_mapping}")
            unknown_ref_mapping = sorted(set(mapping) - set(ref_inventory))
            unknown_model_mapping = sorted(set(mapping.values()) - set(mdl_inventory))
            if unknown_ref_mapping or unknown_model_mapping:
                raise ValueError(
                    f"Mapping uses absent chains: reference={unknown_ref_mapping}, "
                    f"model={unknown_model_mapping}"
                )
            if len(mapping.values()) != len(set(mapping.values())):
                raise ValueError("chain_mapping must be one-to-one")
            short_binders = [
                chain
                for binder in binders
                for chain in binder["chains"]
                if ref_inventory[chain] < 6
            ]
            if short_binders:
                item["warnings"].append(
                    "Short intended binder chain(s) require --min-pep-length 1 and mapping review: "
                    + ", ".join(short_binders)
                )
            if sample.get("nonpolymer_ligands"):
                item["errors"].append(
                    "Non-polymer ligands are annotated, but their inclusion requires explicit user "
                    "confirmation and a separate compare-ligand-structures branch"
                )
            item.update(
                {
                    "reference": reference_text,
                    "reference_chain_inventory": ref_inventory,
                    "model_chain_inventory": mdl_inventory,
                    "chain_mapping": mapping,
                    "short_binder_chains": short_binders,
                }
            )
        except (KeyError, OSError, RuntimeError, ValueError) as exc:
            item["errors"].append(str(exc))
        report["predictions"].append(item)
        report["errors"].extend(
            f"{prediction_path}: {message}" for message in item["errors"]
        )
        report["warnings"].extend(
            f"{prediction_path}: {message}" for message in item["warnings"]
        )
    report["status"] = "READY" if not report["errors"] else "BLOCKED"
    return report, annotations


def write_skipped_summary(
    evaluation_dir: Path,
    sample_id: str,
    prediction_path: str,
    sample: dict[str, Any],
    prediction: dict[str, Any],
) -> Path:
    relative_prediction = Path(prediction_path).relative_to("predicted")
    output = evaluation_dir / "results" / relative_prediction / "skipped.json"
    reason = prediction.get("skip_reason") or sample.get("skip_reason")
    write_json(
        output,
        {
            "schema_version": 1,
            "sample_id": sample_id,
            "prediction_path": prediction_path,
            "status": "SKIPPED",
            "skip_reason": reason,
            "warnings": [f"Intentionally skipped: {reason}"],
            "unavailable_metrics": [
                {"metric": "all", "reason": f"Intentionally skipped: {reason}"}
            ],
        },
    )
    return output


def command_version(command: list[str]) -> str:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    text = (result.stdout + result.stderr).strip()
    if result.returncode:
        raise RuntimeError(f"Version command failed ({shlex.join(command)}): {text}")
    return text


def tool_versions(chirality_script: Path, custom_script: Path) -> dict[str, Any]:
    batch_script = Path(__file__).resolve()
    return {
        "openstructure": command_version(["conda", "run", "-n", "ost212", "ost", "--version"]),
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "gemmi": gemmi.__version__,
        "chirality_script": str(chirality_script.resolve()),
        "chirality_script_sha256": sha256(chirality_script),
        "custom_rmsd_script": str(custom_script.resolve()),
        "custom_rmsd_script_sha256": sha256(custom_script),
        "batch_script": str(batch_script),
        "batch_script_sha256": sha256(batch_script),
    }


def run_command(
    command: list[str], run_dir: Path, label: str, command_log: Path
) -> subprocess.CompletedProcess[str]:
    with command_log.open("a", encoding="utf-8") as handle:
        handle.write(f"\n[{label}]\n{shlex.join(command)}\n")
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    (run_dir / f"{label}.stdout.txt").write_text(result.stdout, encoding="utf-8")
    (run_dir / f"{label}.stderr.txt").write_text(result.stderr, encoding="utf-8")
    with command_log.open("a", encoding="utf-8") as handle:
        handle.write(f"exit_code={result.returncode}\n")
    return result


def ordered_mapping(reference: Path, id_type: str, mapping: dict[str, str]) -> list[tuple[str, str]]:
    inventory = chain_inventory(reference, id_type)
    return [(chain, mapping[chain]) for chain in inventory if chain in mapping]


def ost_summary(raw: dict[str, Any], scored_reference_chains: set[str]) -> dict[str, Any]:
    interfaces = []
    mapped_interfaces = raw.get("dockq_interfaces", [])
    keys = ("dockq", "fnat", "irmsd", "lrmsd", "fnonnat", "nnat", "nmdl")
    for index, identity in enumerate(mapped_interfaces):
        row: dict[str, Any] = {"identity": identity}
        for key in keys:
            values = raw.get(key)
            row[key] = values[index] if isinstance(values, list) and index < len(values) else None
        interfaces.append(row)
    reference_interfaces = raw.get("dockq_reference_interfaces", [])
    extra = [
        pair for pair in reference_interfaces if not set(pair).issubset(scored_reference_chains)
    ]
    return {
        "status": raw.get("status"),
        "ost_version": raw.get("ost_version"),
        "lddt": raw.get("lddt"),
        "ilddt": raw.get("ilddt"),
        "local_lddt": raw.get("local_lddt"),
        "chain_mapping": raw.get("chain_mapping"),
        "alignments": raw.get("aln"),
        "inconsistent_residues": raw.get("inconsistent_residues"),
        "cleanup": {
            key: raw.get(key)
            for key in (
                "model_clashes",
                "model_bad_bonds",
                "model_bad_angles",
                "reference_clashes",
                "reference_bad_bonds",
                "reference_bad_angles",
            )
        },
        "dockq_reference_interfaces": reference_interfaces,
        "per_interface": interfaces,
        "dockq_ave": raw.get("dockq_ave"),
        "dockq_wave": raw.get("dockq_wave"),
        "dockq_ave_full": raw.get("dockq_ave_full"),
        "dockq_wave_full": raw.get("dockq_wave_full"),
        "extra_reference_interfaces": extra,
        "arguments": raw.get("arguments"),
        "log": raw.get("log"),
        "traceback": raw.get("traceback"),
    }


def chirality_summary(path: Path, exit_code: int) -> dict[str, Any]:
    if not path.is_file():
        return {"status": "UNAVAILABLE", "reason": f"Tool exited {exit_code} without summary JSON"}
    raw = load_json(path)
    pooled = raw.get("pooled", {})
    all_selected = pooled.get("all_selected", {})
    status = "SUCCESS"
    if not all_selected.get("n_sites"):
        status = "UNAVAILABLE"
    elif all_selected.get("n_unassessable"):
        status = "PARTIAL"
    return {
        "status": status,
        "primary_scope": "all eligible D- and L-amino-acid residues in annotated chirality chains",
        "primary": all_selected,
        "daa": pooled.get("DAA"),
        "laa": pooled.get("LAA"),
        "excluded_unsupported_residues": raw.get("excluded_unsupported_residues"),
        "tool_exit_code": exit_code,
    }


def run_one(
    evaluation_dir: Path,
    annotation_path: Path,
    sample_id: str,
    sample: dict[str, Any],
    prediction_path: str,
    prediction: dict[str, Any],
    audit_item: dict[str, Any],
    chirality_script: Path,
    custom_script: Path,
    versions: dict[str, Any],
    resume_incomplete: bool,
) -> dict[str, Any]:
    reference = resolved(evaluation_dir, sample["reference"])
    model = resolved(evaluation_dir, prediction_path)
    mapping = prediction["chain_mapping"]
    mapping_pairs = ordered_mapping(reference, sample["reference_chain_id_type"], mapping)
    min_pep_length = 1 if audit_item["short_binder_chains"] else 6
    settings = {
        "ost": {
            "lddt": True,
            "local_lddt": True,
            "ilddt": True,
            "dockq": True,
            "min_pep_length": min_pep_length,
        },
        "custom_rmsd": {
            "alignment_atoms": ["N", "CA", "C", "O"],
            "target_aligned_binder_measurements": [["CA"], ["N", "CA", "C", "O"]],
            "alignment_scope": "all annotated target chains",
            "binder_self_aligned_measurements": [["CA"], ["N", "CA", "C", "O"]],
            "binder_self_alignment_scope": "each annotated binder group independently",
            "missing_policy": "aligned atom intersection; each fit requires 3 non-collinear atom pairs",
        },
        "chirality": {"threshold_deg": 15.0, "primary": "pooled all_selected"},
    }
    fingerprint_payload = {
        "sample_id": sample_id,
        "prediction_path": prediction_path,
        "reference_sha256": sha256(reference),
        "prediction_sha256": sha256(model),
        "sample_annotation": sample,
        "prediction_annotation": prediction,
        "settings": settings,
        "tool_versions": versions,
    }
    fingerprint = hashlib.sha256(
        json.dumps(fingerprint_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    run_dir = evaluation_dir / "results" / Path(prediction_path).relative_to("predicted") / "runs" / fingerprint[:16]
    existing_summary = run_dir / "summary.json"
    if existing_summary.is_file() and load_json(existing_summary).get("status") == "SUCCESS":
        return load_json(existing_summary)
    if run_dir.exists() and not resume_incomplete:
        raise RuntimeError(
            f"Incomplete run already exists: {run_dir}; inspect it and use --resume-incomplete if safe"
        )
    run_dir.mkdir(parents=True, exist_ok=True)
    command_log = run_dir / "commands.txt"
    command_log.write_text(
        "# Exact commands for this prediction\n"
        f"# sample_id={sample_id}\n# prediction={model}\n# reference={reference}\n"
        f"# annotations={annotation_path.resolve()}\n"
        f"# tool_versions={json.dumps(versions, sort_keys=True)}\n",
        encoding="utf-8",
    )
    manifest = {
        "schema_version": 1,
        "fingerprint": fingerprint,
        "inputs": fingerprint_payload,
        "run_dir": str(run_dir.resolve()),
    }
    write_json(run_dir / "run_manifest.json", manifest)

    scored_reference_chains = set(sample["target_chains"])
    scored_reference_chains.update(
        chain for binder in sample["binders"] for chain in binder["chains"]
    )
    scored_model_chains = {mapping[chain] for chain in scored_reference_chains}
    reference_subset = run_dir / "reference_scored_chains.cif"
    model_subset = run_dir / "model_scored_chains.cif"
    subset_manifest = {
        "reference": write_polymer_subset(
            reference,
            reference_subset,
            sample["reference_chain_id_type"],
            scored_reference_chains,
        ),
        "model": write_polymer_subset(
            model,
            model_subset,
            prediction["model_chain_id_type"],
            scored_model_chains,
        ),
    }
    write_json(run_dir / "subset_manifest.json", subset_manifest)
    with command_log.open("a", encoding="utf-8") as handle:
        handle.write(
            "\n[custom_polymer_subset]\n"
            f"reference_input={reference}\nreference_output={reference_subset}\n"
            f"reference_chains={','.join(sorted(scored_reference_chains))}\n"
            f"model_input={model}\nmodel_output={model_subset}\n"
            f"model_chains={','.join(sorted(scored_model_chains))}\n"
        )

    ost_path = run_dir / "ost_scores.json"
    ost_command = [
        "conda", "run", "-n", "ost212", "ost", "compare-structures",
        "-m", str(model_subset), "-r", str(reference_subset), "-c",
        *[f"{ref}:{mdl}" for ref, mdl in mapping_pairs],
        "--lddt", "--local-lddt", "--ilddt", "--dockq", "-o", str(ost_path),
    ]
    if min_pep_length == 1:
        ost_command.extend(["--min-pep-length", "1"])
    ost_run = run_command(ost_command, run_dir, "openstructure", command_log)
    raw_ost = load_json(ost_path) if ost_path.is_file() else {
        "status": "FAILURE", "traceback": f"No output JSON; exit code {ost_run.returncode}"
    }
    parsed_ost = ost_summary(raw_ost, scored_reference_chains)

    custom_config = {
        "reference": str(reference_subset),
        "model": str(model_subset),
        "reference_chain_id_type": sample["reference_chain_id_type"],
        "model_chain_id_type": prediction["model_chain_id_type"],
        "chain_mapping": mapping,
        "target_chains": sample["target_chains"],
        "binders": sample["binders"],
    }
    custom_config_path = run_dir / "custom_rmsd_config.json"
    write_json(custom_config_path, custom_config)
    custom_path = run_dir / "binder_rmsd.json"
    custom_command = [
        sys.executable, str(custom_script), "--config", str(custom_config_path),
        "--ost-json", str(ost_path), "--output", str(custom_path),
    ]
    custom_run = run_command(custom_command, run_dir, "custom_rmsd", command_log)
    custom = load_json(custom_path) if custom_path.is_file() else {
        "status": "UNAVAILABLE", "reason": f"No output JSON; exit code {custom_run.returncode}"
    }
    target_aligned = custom.get("target_aligned") or {
        "status": "UNAVAILABLE", "reason": custom.get("reason", "No target-aligned result")
    }
    self_aligned = custom.get("self_aligned") or {
        "status": "UNAVAILABLE", "reason": custom.get("reason", "No self-aligned result")
    }

    chirality_dir = run_dir / "chirality"
    predicted_chirality_chains = [mapping[chain] for chain in sample["chirality_chains"]]
    chirality_command = [
        sys.executable, str(chirality_script), str(model_subset),
        *[value for chain in predicted_chirality_chains for value in ("--chain", chain)],
        "--out-dir", str(chirality_dir), "--threshold-deg", "15", "--detailed",
    ]
    chirality_run = run_command(chirality_command, run_dir, "chirality", command_log)
    chirality = chirality_summary(
        chirality_dir / "all_residue_chirality_summary.json", chirality_run.returncode
    )

    warnings = list(audit_item["warnings"])
    if parsed_ost.get("chain_mapping") != mapping:
        warnings.append(
            f"OpenStructure mapping differs from annotation: {parsed_ost.get('chain_mapping')}"
        )
    if parsed_ost["extra_reference_interfaces"]:
        warnings.append(
            "OpenStructure reported interfaces involving unscored reference chains: "
            + repr(parsed_ost["extra_reference_interfaces"])
        )
    if audit_item["short_binder_chains"] and custom.get("status") == "SUCCESS":
        questionable = []
        for chain in audit_item["short_binder_chains"]:
            details = custom.get("residue_correspondence", {}).get("chains", {}).get(chain, {})
            if details.get("sequence_mismatches") or not details.get("aligned_residue_pairs"):
                questionable.append(chain)
        if questionable:
            warnings.append("Questionable very-short-chain match: " + ", ".join(questionable))

    statuses = [parsed_ost.get("status"), custom.get("status"), chirality.get("status")]
    status = "SUCCESS" if statuses == ["SUCCESS", "SUCCESS", "SUCCESS"] else "PARTIAL"
    if parsed_ost.get("status") != "SUCCESS":
        status = "FAILURE"
    unavailable_metrics = [
        {"metric": "nonpolymer_ligand", "reason": "Not requested/confirmed for this batch"}
    ]
    for section_name, section in (
        ("target_aligned_binder_rmsd", target_aligned),
        ("binder_self_aligned_rmsd", self_aligned),
    ):
        for binder in section.get("binders", []):
            for metric_name, metric in binder.items():
                if isinstance(metric, dict) and metric.get("status") == "UNAVAILABLE":
                    unavailable_metrics.append(
                        {
                            "metric": f"{section_name}.{binder.get('id')}.{metric_name}",
                            "reason": metric.get("reason", "Unavailable"),
                        }
                    )
    summary = {
        "schema_version": 2,
        "status": status,
        "sample_id": sample_id,
        "prediction_path": prediction_path,
        "reference_path": sample["reference"],
        "reference_assembly": sample["reference_assembly"],
        "reference_chain_id_type": sample["reference_chain_id_type"],
        "model_chain_id_type": prediction["model_chain_id_type"],
        "chain_mapping_used": mapping,
        "target_chains": sample["target_chains"],
        "binders": sample["binders"],
        "context_chains": sample.get("context_chains", []),
        "tool_versions": versions,
        "run_fingerprint": fingerprint,
        "run_dir": str(run_dir.resolve()),
        "warnings": warnings,
        "openstructure": parsed_ost,
        "target_aligned_binder_rmsd": target_aligned,
        "binder_self_aligned_rmsd": self_aligned,
        "rmsd_residue_correspondence": custom.get("residue_correspondence"),
        "chirality": chirality,
        "unavailable_metrics": unavailable_metrics,
    }
    write_json(existing_summary, summary)
    return summary


def regenerate_batch_csv(
    evaluation_dir: Path, active_run_fingerprints: set[str] | None = None
) -> Path:
    rows: list[dict[str, Any]] = []
    results_dir = evaluation_dir / "results"
    summary_paths = list(results_dir.glob("**/runs/*/summary.json"))
    summary_paths.extend(results_dir.glob("**/skipped.json"))
    for path in sorted(summary_paths):
        summary = load_json(path)
        fingerprint = summary.get("run_fingerprint")
        if fingerprint and active_run_fingerprints is not None and fingerprint not in active_run_fingerprints:
            continue
        target_aligned = summary.get("target_aligned_binder_rmsd", {})
        target_binders = {item["id"]: item for item in target_aligned.get("binders", [])}
        self_aligned = summary.get("binder_self_aligned_rmsd", {})
        self_binders = {item["id"]: item for item in self_aligned.get("binders", [])}
        binder_defs = summary.get("binders") or [{"id": "", "chains": []}]
        for binder_def in binder_defs:
            target_binder = target_binders.get(binder_def["id"], {})
            target_ca = target_binder.get("target_aligned_binder_ca_rmsd_angstrom", {})
            target_backbone = target_binder.get(
                "target_aligned_binder_backbone_rmsd_angstrom", {}
            )
            self_binder = self_binders.get(binder_def["id"], {})
            self_ca = self_binder.get("binder_self_aligned_ca_rmsd_angstrom", {})
            self_backbone = self_binder.get("binder_self_aligned_backbone_rmsd_angstrom", {})
            ost = summary.get("openstructure", {})
            chirality = summary.get("chirality", {})
            primary = chirality.get("primary") or {}
            daa = chirality.get("daa") or {}
            laa = chirality.get("laa") or {}
            rows.append(
                {
                    "sample_id": summary.get("sample_id"),
                    "prediction_path": summary.get("prediction_path"),
                    "run_fingerprint": summary.get("run_fingerprint"),
                    "status": summary.get("status"),
                    "ost_status": ost.get("status"),
                    "target_aligned_binder_rmsd_status": target_aligned.get("status"),
                    "target_aligned_binder_ca_rmsd_status": target_ca.get("status"),
                    "target_aligned_binder_backbone_rmsd_status": target_backbone.get("status"),
                    "binder_self_aligned_rmsd_status": self_aligned.get("status"),
                    "binder_self_aligned_ca_rmsd_status": self_ca.get("status"),
                    "binder_self_aligned_backbone_rmsd_status": self_backbone.get("status"),
                    "chirality_status": chirality.get("status"),
                    "binder_id": binder_def["id"],
                    "binder_reference_chains": "+".join(binder_def["chains"]),
                    "ost_lddt": ost.get("lddt"),
                    "ost_ilddt": ost.get("ilddt"),
                    "ost_dockq_ave": ost.get("dockq_ave"),
                    "ost_dockq_wave": ost.get("dockq_wave"),
                    "ost_dockq_ave_full": ost.get("dockq_ave_full"),
                    "ost_dockq_wave_full": ost.get("dockq_wave_full"),
                    "ost_per_interface_json": json.dumps(ost.get("per_interface"), separators=(",", ":")),
                    "target_backbone_fit_rmsd_angstrom": target_aligned.get("alignment", {}).get(
                        "target_backbone_rmsd_angstrom"
                    ),
                    "target_aligned_binder_ca_rmsd_angstrom": target_ca.get("rmsd_angstrom"),
                    "target_aligned_binder_backbone_rmsd_angstrom": target_backbone.get(
                        "rmsd_angstrom"
                    ),
                    "binder_self_aligned_ca_rmsd_angstrom": self_ca.get("rmsd_angstrom"),
                    "binder_self_aligned_backbone_rmsd_angstrom": self_backbone.get(
                        "rmsd_angstrom"
                    ),
                    "chirality_all_correct_count": primary.get("n_pass"),
                    "chirality_all_eligible_count": primary.get("n_sites"),
                    "chirality_all_skipped_count": primary.get("n_unassessable"),
                    "chirality_all_accuracy_pct": primary.get("chirality_correct_rate_pct"),
                    "chirality_daa_correct_count": daa.get("n_pass"),
                    "chirality_daa_eligible_count": daa.get("n_sites"),
                    "chirality_daa_skipped_count": daa.get("n_unassessable"),
                    "chirality_daa_accuracy_pct": daa.get("chirality_correct_rate_pct"),
                    "chirality_laa_correct_count": laa.get("n_pass"),
                    "chirality_laa_eligible_count": laa.get("n_sites"),
                    "chirality_laa_skipped_count": laa.get("n_unassessable"),
                    "chirality_laa_accuracy_pct": laa.get("chirality_correct_rate_pct"),
                    "warnings": " | ".join(summary.get("warnings", [])),
                    "unavailable_metrics_json": json.dumps(
                        summary.get("unavailable_metrics", []), separators=(",", ":")
                    ),
                    "summary_path": str(path.resolve()),
                }
            )
    output = results_dir / "batch_summary.csv"
    if not rows:
        return output
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, output)
    return output


def main() -> int:
    args = parse_args()
    evaluation_dir = args.evaluation_dir.resolve()
    annotation_path = (args.annotations or evaluation_dir / "annotations" / "samples.json").resolve()
    report, annotations = audit(evaluation_dir, annotation_path)
    print(json.dumps(report, indent=2, allow_nan=False))
    if args.audit_output:
        write_json(args.audit_output.resolve(), report)
    if not args.execute:
        return 0 if report["status"] == "READY" else 2
    if report["status"] != "READY":
        print("Audit is blocked; no evaluation commands were run.", file=sys.stderr)
        return 2
    if not args.chirality_script or not args.chirality_script.is_file():
        print("--chirality-script must name the inspected repository script", file=sys.stderr)
        return 2
    custom_script = Path(__file__).with_name("target_aligned_binder_rmsd.py")
    try:
        versions = tool_versions(args.chirality_script.resolve(), custom_script)
        records = prediction_records(annotations)
        audit_by_prediction = {item["prediction"]: item for item in report["predictions"]}
        active_run_fingerprints: set[str] = set()
        for prediction_path in report["discovered_predictions"]:
            sample_id, sample, prediction = records[prediction_path]
            if sample.get("skip") or prediction.get("skip"):
                write_skipped_summary(
                    evaluation_dir, sample_id, prediction_path, sample, prediction
                )
                continue
            result = run_one(
                evaluation_dir,
                annotation_path,
                sample_id,
                sample,
                prediction_path,
                prediction,
                audit_by_prediction[prediction_path],
                args.chirality_script.resolve(),
                custom_script,
                versions,
                args.resume_incomplete,
            )
            if result.get("run_fingerprint"):
                active_run_fingerprints.add(result["run_fingerprint"])
        output = regenerate_batch_csv(evaluation_dir, active_run_fingerprints)
        print(f"Batch summary: {output}")
    except Exception as exc:
        print(f"Error: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
