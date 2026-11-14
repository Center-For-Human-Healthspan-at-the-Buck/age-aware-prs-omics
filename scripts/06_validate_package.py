#!/usr/bin/env python3
"""Validate figure paths, Arivale cohort size and package hashes."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARIVALE_N = 1491
REQUIRED = (
    "assets/fig06a_disease_association_workflow.pdf",
    "figures/Extended_Data_Fig4_AoU_directional_concordance.pdf",
    "figures/Extended_Data_Fig4_AoU_directional_concordance_coverage.tsv",
    "config/vae_age__run__id_assoc_bayes.yaml",
    "config/vae_resid_run__id_assoc_bayes.yaml",
    "scripts/aou_workbench/src/aou_replication_utils.py",
)
OBSOLETE = (
    "assets/fig07a_validation_workflow.pdf",
    "figures/Extended_Data_Fig3_AoU_replication_comparison.pdf",
)
EXPECTED = {
    ("OLS", "proteomics"): (16, 12),
    ("OLS", "chemistries"): (26, 16),
    ("VAE", "proteomics"): (81, 21),
    ("VAE", "chemistries"): (43, 21),
}


def check_inputs() -> None:
    for relative in REQUIRED:
        if not (ROOT / relative).is_file():
            raise ValueError(f"Missing required file: {relative}")
    for relative in OBSOLETE:
        if (ROOT / relative).exists():
            raise ValueError(f"Obsolete figure path: {relative}")
    for name in ("model_comparison_hits.tsv", "ols_positive_controls.tsv"):
        path = ROOT / "data" / "manuscript" / name
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if "N" not in rows[0]:
            raise ValueError(f"{name} has no N column")
        incorrect = [row["N"] for row in rows if row["N"] and row["N"] != str(ARIVALE_N)]
        if incorrect:
            raise ValueError(f"{name}: discovery cohort N differs from {ARIVALE_N}")

    path = ROOT / REQUIRED[2]
    with path.open(newline="", encoding="utf-8") as handle:
        coverage = list(csv.DictReader(handle, delimiter="\t"))
    observed = {
        (row["method"], row["omic_type"]): (
            int(row["n_modelled"]), int(row["n_concordant"])
        )
        for row in coverage
    }
    if observed != EXPECTED:
        raise ValueError(f"Extended Data Fig. 4 counts changed: {observed}")
    print("PASS: figure paths, Arivale cohort N, and AoU concordance counts")


def package_files() -> list[Path]:
    excluded = {"outputs", ".git", "__pycache__", ".ipynb_checkpoints"}
    files = []
    for path in ROOT.rglob("*"):
        if any(part in excluded or part.startswith("__MACOSX") for part in path.relative_to(ROOT).parts):
            continue
        if path.is_symlink():
            raise ValueError(f"Unexpected symlink: {path}")
        if path.is_file() and path.name not in {"MANIFEST.tsv", ".DS_Store"}:
            files.append(path)
    return sorted(files, key=lambda path: path.relative_to(ROOT).as_posix())


def generate_manifest() -> list[tuple[str, str, str]]:
    return [
        (
            path.relative_to(ROOT).as_posix(),
            str(path.stat().st_size),
            hashlib.sha256(path.read_bytes()).hexdigest(),
        )
        for path in package_files()
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-manifest", action="store_true")
    options = parser.parse_args()
    check_inputs()
    path = ROOT / "MANIFEST.tsv"
    rows = generate_manifest()
    if options.write_manifest:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(("path", "bytes", "sha256"))
            writer.writerows(rows)
        print(f"Wrote {len(rows)} file hashes to {path.name}")
    else:
        with path.open(newline="", encoding="utf-8") as handle:
            saved = list(csv.reader(handle, delimiter="\t"))
        if saved != [["path", "bytes", "sha256"], *[list(row) for row in rows]]:
            raise ValueError("MANIFEST.tsv does not match the package; regenerate it with --write-manifest")
        print(f"PASS: {len(rows)} package file hashes match MANIFEST.tsv")


if __name__ == "__main__":
    main()
