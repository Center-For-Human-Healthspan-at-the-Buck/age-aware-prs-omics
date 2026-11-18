#!/usr/bin/env python3
"""Assemble vector panel PDFs into final Main and Extended Data figures."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import fitz
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PT_PER_MM = 72.0 / 25.4
MAIN_MAX_MM = (183.0, 170.0)
ED_MAX_MM = (180.0, 170.0)


@dataclass(frozen=True)
class Panel:
    source: str
    slot_mm: tuple[float, float, float, float]
    label: str | None = None
    label_mm: tuple[float, float] | None = None
    trim: bool = True


@dataclass(frozen=True)
class Layout:
    output: str
    page_mm: tuple[float, float]
    source_dir: str
    panels: tuple[Panel, ...]
    labels: tuple[tuple[str, float, float], ...] = ()
    kind: str = "main"


LAYOUTS = (
    Layout(
        output="Fig_1.pdf",
        page_mm=(183.0, 92.5),
        source_dir="outputs/main",
        panels=(Panel("fig01_workflow_schematic.pdf", (0, 0, 183, 92.5), trim=False),),
    ),
    Layout(
        output="Fig_2.pdf",
        page_mm=(183.0, 170.0),
        source_dir="outputs/main",
        panels=(
            Panel("fig02a_age_residualized_chord.pdf", (3, 2, 87, 55), "a", (2, 5)),
            Panel("fig02b_age_included_chord.pdf", (93, 2, 87, 55), "b", (92, 5)),
            Panel("fig02c_age_residualized_top_features.pdf", (3, 58, 177, 54), "c", (2, 61)),
            Panel("fig02d_age_included_top_features.pdf", (3, 113, 177, 54), "d", (2, 116)),
        ),
    ),
    Layout(
        output="Fig_3.pdf",
        page_mm=(183.0, 170.0),
        source_dir="outputs/main",
        panels=(
            Panel("fig03a_neurodegeneration_fingerprint.pdf", (3, 2, 177, 75), "a", (2, 5)),
            Panel("fig03b_pd_feature_by_model_matrix.pdf", (3, 79, 87, 76), "b", (2, 82)),
            Panel("fig03c_als_feature_by_model_matrix.pdf", (93, 79, 87, 76), "c", (92, 82)),
            Panel("fig03bc_shared_legend.pdf", (20, 157, 143, 11)),
        ),
    ),
    Layout(
        output="Fig_4.pdf",
        page_mm=(183.0, 145.0),
        source_dir="outputs/main",
        panels=(
            Panel("fig04a_ad_lollipop.pdf", (3, 7, 104, 119), "a", (2, 7)),
            Panel("fig04b_ad_age_context.pdf", (110, 4, 70, 124), "b", (109, 7)),
            Panel("fig04b_legend.pdf", (24, 132, 135, 11)),
        ),
    ),
    Layout(
        output="Fig_5.pdf",
        page_mm=(183.0, 160.0),
        source_dir="outputs/main",
        panels=(
            Panel("fig05a_t2d_lollipop.pdf", (3, 2, 87, 84), "a", (2, 5)),
            Panel("fig05b_cvd_lollipop.pdf", (93, 2, 87, 84), "b", (92, 5)),
            Panel("fig05c_cvd_recurring_associations.pdf", (42, 88, 99, 69), "c", (41, 91)),
        ),
    ),
    Layout(
        output="Fig_6.pdf",
        page_mm=(183.0, 170.0),
        source_dir="outputs/main",
        panels=(
            Panel("fig06a1_disease_association_workflow.pdf", (2, 2, 101, 38)),
            Panel("fig06a2_chemistry_ancestry.pdf", (105, 2, 36, 29)),
            Panel("fig06a3_proteomics_ancestry.pdf", (144, 2, 36, 29)),
            Panel("fig06a_ancestry_legend.pdf", (105, 31, 75, 9)),
            Panel("fig06b1_clinical_directional.pdf", (3, 43, 87, 60)),
            Panel("fig06b2_protein_directional.pdf", (93, 43, 87, 60)),
            Panel("fig06c1_clinical_discordant.pdf", (3, 107, 87, 60)),
            Panel("fig06c2_protein_discordant.pdf", (93, 107, 87, 60)),
        ),
        labels=(("a", 2, 5),("b", 2, 46), ("c", 2, 110)),
    ),
    Layout(
        output="Extended_Data_Fig1_Domain_level_model_comparison.pdf",
        page_mm=(180.0, 157.7),
        source_dir="outputs/extended_data",
        panels=(
            Panel("extended_data_fig01a_domain_model_counts.pdf", (3, 2, 102, 75), "a", (2, 5)),
            Panel("extended_data_fig01b_pairwise_overlap.pdf", (108, 2, 69, 75), "b", (107, 5)),
            Panel("extended_data_fig01c_age_context.pdf", (3, 80, 174, 74), "c", (2, 83)),
        ),
        kind="extended",
    ),
    Layout(
        output="Extended_Data_Fig2_Vae_Positive_Control_PRS_Associations.pdf",
        page_mm=(180.0, 82.0),
        source_dir="outputs/extended_data",
        panels=(
            Panel("extended_data_fig02a_positive_controls_age_residualized.pdf", (3, 2, 86, 77), "a", (2, 5)),
            Panel("extended_data_fig02b_positive_controls_age_included.pdf", (92, 2, 85, 77), "b", (91, 5)),
        ),
        kind="extended",
    ),
    Layout(
        output="Extended_Data_Fig3_VAE_longevity_associations.pdf",
        page_mm=(180.0, 170.0),
        source_dir="outputs/extended_data",
        panels=(Panel("extended_data_fig03_longevity_lollipop.pdf", (3, 2, 174, 166)),),
        kind="extended",
    ),
)


def mm_rect(values: tuple[float, float, float, float]) -> fitz.Rect:
    x, y, width, height = values
    return fitz.Rect(
        x * PT_PER_MM,
        y * PT_PER_MM,
        (x + width) * PT_PER_MM,
        (y + height) * PT_PER_MM,
    )


def content_clip(page: fitz.Page, padding_pt: float = 2.0) -> fitz.Rect:
    scale = 2.0
    pixmap = page.get_pixmap(matrix=fitz.Matrix(scale, scale), alpha=False)
    pixels = np.frombuffer(pixmap.samples, dtype=np.uint8).reshape(
        pixmap.height, pixmap.width, pixmap.n
    )[..., :3]
    mask = np.any(pixels < 248, axis=2)
    if not mask.any():
        return page.rect
    rows, columns = np.where(mask)
    clip = fitz.Rect(
        columns.min() / scale - padding_pt,
        rows.min() / scale - padding_pt,
        (columns.max() + 1) / scale + padding_pt,
        (rows.max() + 1) / scale + padding_pt,
    )
    return clip & page.rect


def add_label(page: fitz.Page, label: str, x_mm: float, y_mm: float) -> None:
    page.insert_text(
        fitz.Point(x_mm * PT_PER_MM, y_mm * PT_PER_MM),
        label,
        fontname="hebo",
        fontsize=8,
        color=(0, 0, 0),
        overlay=True,
    )


def assemble(layout: Layout, output_dir: Path) -> Path:
    width_mm, height_mm = layout.page_mm
    max_width, max_height = MAIN_MAX_MM if layout.kind == "main" else ED_MAX_MM
    if width_mm > max_width + 0.01 or height_mm > max_height + 0.01:
        raise ValueError(
            f"{layout.output}: {width_mm:g} x {height_mm:g} mm exceeds "
            f"{max_width:g} x {max_height:g} mm"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    result = fitz.open()
    target = result.new_page(width=width_mm * PT_PER_MM, height=height_mm * PT_PER_MM)
    source_root = ROOT / layout.source_dir

    for panel in layout.panels:
        source_path = source_root / panel.source
        if not source_path.is_file() or source_path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing panel PDF: {source_path.relative_to(ROOT)}")
        source = fitz.open(source_path)
        if source.page_count != 1:
            source.close()
            raise ValueError(f"Panel PDF must contain one page: {source_path.relative_to(ROOT)}")
        source_page = source[0]
        clip = content_clip(source_page) if panel.trim else source_page.rect
        target.show_pdf_page(
            mm_rect(panel.slot_mm),
            source,
            0,
            clip=clip,
            keep_proportion=True,
            overlay=True,
        )
        source.close()
        if panel.label and panel.label_mm:
            add_label(target, panel.label, *panel.label_mm)

    for label, x_mm, y_mm in layout.labels:
        add_label(target, label, x_mm, y_mm)

    result.set_metadata({
        "title": layout.output.removesuffix(".pdf").replace("_", " "),
        "author": "Nogal et al.",
    })
    destination = output_dir / layout.output
    result.save(destination, garbage=4, deflate=True, clean=True)
    result.close()
    return destination


def validate_pdf(path: Path, layout: Layout) -> None:
    document = fitz.open(path)
    if document.page_count != 1:
        document.close()
        raise ValueError(f"{path.name}: expected one page")
    page = document[0]
    width_mm = page.rect.width / PT_PER_MM
    height_mm = page.rect.height / PT_PER_MM
    expected_width, expected_height = layout.page_mm
    document.close()
    if abs(width_mm - expected_width) > 0.05 or abs(height_mm - expected_height) > 0.05:
        raise ValueError(
            f"{path.name}: expected {expected_width:g} x {expected_height:g} mm, "
            f"found {width_mm:.2f} x {height_mm:.2f} mm"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--main-only", action="store_true")
    group.add_argument("--extended-only", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "figures")
    options = parser.parse_args()

    selected = LAYOUTS
    if options.main_only:
        selected = tuple(layout for layout in LAYOUTS if layout.kind == "main")
    elif options.extended_only:
        selected = tuple(layout for layout in LAYOUTS if layout.kind == "extended")

    for layout in selected:
        destination = assemble(layout, options.output_dir)
        validate_pdf(destination, layout)
        print(f"Saved {destination.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
