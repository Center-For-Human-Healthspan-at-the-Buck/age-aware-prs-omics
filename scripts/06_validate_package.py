#!/usr/bin/env python3
"""Validate figure paths, Arivale cohort size and package hashes."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

import fitz


ROOT = Path(__file__).resolve().parents[1]
ARIVALE_N = 1491
COVERAGE_PATH = "figures/Extended_Data_Fig4_AoU_directional_concordance_coverage.tsv"
FINAL_PDFS = {
    "figures/Fig_1.pdf": (183.0, 170.0),
    "figures/Fig_2.pdf": (183.0, 170.0),
    "figures/Fig_3.pdf": (183.0, 170.0),
    "figures/Fig_4.pdf": (183.0, 170.0),
    "figures/Fig_5.pdf": (183.0, 170.0),
    "figures/Fig_6.pdf": (183.0, 170.0),
    "figures/Extended_Data_Fig1_Domain_level_model_comparison.pdf": (180.0, 170.0),
    "figures/Extended_Data_Fig2_Vae_Positive_Control_PRS_Associations.pdf": (180.0, 170.0),
    "figures/Extended_Data_Fig3_VAE_longevity_associations.pdf": (180.0, 170.0),
    "figures/Extended_Data_Fig4_AoU_directional_concordance.pdf": (180.0, 170.0),
}
REQUIRED = (
    "assets/fig01_workflow_schematic.pdf",
    "assets/fig06a_disease_association_workflow.pdf",
    COVERAGE_PATH,
    "config/vae_age__run__id_assoc_bayes.yaml",
    "config/vae_resid_run__id_assoc_bayes.yaml",
    "scripts/aou_workbench/src/aou_replication_utils.py",
    *FINAL_PDFS,
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

    path = ROOT / COVERAGE_PATH
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

    points_per_mm = 72.0 / 25.4
    for relative, (max_width_mm, max_height_mm) in FINAL_PDFS.items():
        document = fitz.open(ROOT / relative)
        if document.page_count != 1:
            document.close()
            raise ValueError(f"{relative}: expected one page")
        page = document[0]
        width_mm = page.rect.width / points_per_mm
        height_mm = page.rect.height / points_per_mm
        has_editable_content = bool(page.get_text().strip()) or bool(page.get_drawings())
        document.close()
        if width_mm > max_width_mm + 0.1 or height_mm > max_height_mm + 0.1:
            raise ValueError(
                f"{relative}: {width_mm:.2f} x {height_mm:.2f} mm exceeds "
                f"{max_width_mm:g} x {max_height_mm:g} mm"
            )
        if not has_editable_content:
            raise ValueError(f"{relative}: page appears to contain only raster artwork")
    print("PASS: figure paths, sizes, Arivale cohort N, and AoU concordance counts")


def package_files() -> list[Path]:
    excluded = {"outputs", ".git", "__pycache__", ".ipynb_checkpoints"}
    files = []
    for path in ROOT.rglob("*"):
        if any(part in excluded or part.startswith("__MACOSX") for part in path.relative_to(ROOT).parts):
            continue
        if path.is_symlink():
            raise ValueError(f"Unexpected symlink: {path}")
        if path.is_file() and path.name not in {
            "MANIFEST.tsv",
            ".DS_Store",
            "Rplots.pdf",
        }:
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
