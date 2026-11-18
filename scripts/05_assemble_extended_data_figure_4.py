#!/usr/bin/env python3
"""Assemble Extended Data Fig. 4 from the locked, aggregate AoU results.

This is an independent rendering of the same feature-level summary used by
R/aou_plots.R. It does not require or export participant-level All of Us data.
The upstream file retains its original "replication_*" column names for schema
compatibility; the figure describes directional concordance across different
discovery and external disease-control estimands.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DOMAINS = (
    "T2D",
    "Heart_failure",
    "Parkinsons",
    "CAD_coronary",
    "CVD_composite",
    "AD_dementia",
    "ALS_MND",
)
DOMAIN_LABELS = {
    "T2D": "Type 2 diabetes",
    "Heart_failure": "Heart failure",
    "Parkinsons": "Parkinson's disease",
    "CAD_coronary": "Coronary artery disease",
    "CVD_composite": "Cardiovascular composite",
    "AD_dementia": "Alzheimer's disease / dementia",
    "ALS_MND": "ALS / motor-neuron disease",
}
DOMAIN_COLORS = {
    "T2D": "#F0772D",
    "Heart_failure": "#EB6575",
    "Parkinsons": "#A53379",
    "CAD_coronary": "#278537",
    "CVD_composite": "#C9BC45",
    "AD_dementia": "#4779A9",
    "ALS_MND": "#65C5E4",
}
EXPECTED_COVERAGE = {
    ("OLS", "proteomics"): (16, 12),
    ("OLS", "chemistries"): (26, 16),
    ("VAE", "proteomics"): (81, 21),
    ("VAE", "chemistries"): (43, 21),
}


def canonical_feature(label: str) -> str:
    label = re.sub(r"[\u2010-\u2015\u2212]", "-", str(label))
    return re.sub(r"\s+", " ", label.replace("_", " ")).strip().lower()


def summarize(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {
        "replication_test_family", "discovery_framework", "omic_type",
        "analysis_domain", "uniprot", "feature_display",
        "discovery_concordant_status", "discovery_concordant_q",
        "direction_relation",
    }
    missing = required.difference(raw)
    if missing:
        raise ValueError(f"Missing columns in aggregate AoU data: {sorted(missing)}")

    rows = raw.loc[
        raw["replication_test_family"].eq("disease main effect")
        & raw["discovery_framework"].isin(("VAE", "Univariate"))
        & raw["omic_type"].isin(("proteomics", "chemistries"))
        & raw["analysis_domain"].isin(DOMAINS)
    ].copy()
    rows["method"] = rows["discovery_framework"].replace({"Univariate": "OLS"})
    rows["feature_identity"] = rows.apply(
        lambda row: "uniprot:" + row["uniprot"].strip()
        if row["omic_type"] == "proteomics" and row["uniprot"].strip()
        else "label:" + canonical_feature(row["feature_display"]),
        axis=1,
    )
    rows["modelled"] = rows["discovery_concordant_status"].eq("ok")
    q_value = pd.to_numeric(rows["discovery_concordant_q"], errors="coerce")
    rows["concordant"] = (
        rows["modelled"] & q_value.lt(0.10)
        & rows["direction_relation"].eq("same_direction")
    )
    pairs = rows.groupby(
        ["method", "omic_type", "analysis_domain", "feature_identity"],
        as_index=False,
    )[["modelled", "concordant"]].max()
    coverage = pairs.groupby(["method", "omic_type"], as_index=False).agg(
        n_nominated=("modelled", "size"),
        n_modelled=("modelled", "sum"),
        n_concordant=("concordant", "sum"),
    )
    coverage["concordant_fraction"] = (
        coverage["n_concordant"] / coverage["n_modelled"]
    )
    counts = (
        pairs.loc[pairs["concordant"]]
        .groupby(["method", "omic_type", "analysis_domain"], as_index=False)
        .size()
        .rename(columns={"size": "n_concordant"})
    )
    found = {
        (row.method, row.omic_type): (int(row.n_modelled), int(row.n_concordant))
        for row in coverage.itertuples(index=False)
    }
    if found != EXPECTED_COVERAGE:
        raise ValueError(
            "Locked Fig. 4 totals have changed; investigate before plotting. "
            f"Expected {EXPECTED_COVERAGE}, found {found}."
        )
    return coverage, counts


def plot_domain_bars(
    ax: plt.Axes, counts: pd.DataFrame, layer: str, title: str
) -> None:
    methods = ("OLS", "VAE")
    bottoms = [0, 0]
    for domain in DOMAINS:
        heights = [
            int(counts.loc[
                counts["method"].eq(method)
                & counts["omic_type"].eq(layer)
                & counts["analysis_domain"].eq(domain),
                "n_concordant",
            ].sum())
            for method in methods
        ]
        ax.bar(
            methods, heights, bottom=bottoms, width=0.63,
            color=DOMAIN_COLORS[domain], linewidth=0.5, edgecolor="white",
        )
        bottoms = [bottom + height for bottom, height in zip(bottoms, heights)]
    for i, total in enumerate(bottoms):
        ax.text(i, total + 0.3, str(total), ha="center", va="bottom",
                fontsize=9, fontweight="bold")
    ax.set(title=title, ylim=(0, 25), yticks=(0, 5, 10, 15, 20, 25))
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=9)
    ax.tick_params(axis="y", labelsize=8)


def draw(coverage: pd.DataFrame, counts: pd.DataFrame, destination: Path) -> None:
    fig = plt.figure(figsize=(180 / 25.4, 156.25 / 25.4), facecolor="white")
    top = fig.add_axes((0.14, 0.645, 0.81, 0.255))
    left = fig.add_axes((0.14, 0.155, 0.35, 0.275))
    right = fig.add_axes((0.57, 0.155, 0.38, 0.275), sharey=left)
    fig.text(0.035, 0.945, "a", fontsize=12, fontweight="bold")
    fig.text(0.035, 0.456, "b", fontsize=12, fontweight="bold")
    fig.text(0.535, 0.456, "c", fontsize=12, fontweight="bold")
    fig.text(0.14, 0.945, "Candidate-set size and directional concordance",
             fontsize=9.2, fontweight="bold")

    for row in coverage.itertuples(index=False):
        color = "#E97727" if row.method == "OLS" else "#257EB5"
        marker = "^" if row.omic_type == "proteomics" else "o"
        top.scatter(row.n_modelled, row.concordant_fraction,
                    s=35 * row.n_concordant, marker=marker,
                    color=color, edgecolor="black", linewidth=0.55, zorder=3)
        offsets = {
            ("OLS", "proteomics"): (4.0, 0.015),
            ("OLS", "chemistries"): (4.0, 0.015),
            ("VAE", "chemistries"): (4.0, 0.015),
            ("VAE", "proteomics"): (-17.5, 0.025),
        }
        dx, dy = offsets[(row.method, row.omic_type)]
        top.text(row.n_modelled + dx, row.concordant_fraction + dy,
                 f"{row.method} – {row.omic_type.capitalize()}",
                 fontsize=8.1, va="center")
    top.set(
        xlim=(8, 94), ylim=(0, 0.86), xticks=(20, 40, 60, 80),
        yticks=(0, 0.2, 0.4, 0.6, 0.8),
        xlabel="Unique disease–feature associations modelled in All of Us (n)",
        ylabel="Directionally concordant fraction\n(q < 0.10)",
    )
    top.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    top.xaxis.label.set_size(7.2)
    top.yaxis.label.set_size(7.2)
    top.yaxis.labelpad = 3
    top.spines[["top", "right"]].set_visible(False)
    top.tick_params(labelsize=8)
    fig.text(0.54, 0.557, "Directionally concordant associations (n)",
             fontsize=7.0, ha="center")
    size_ax = fig.add_axes((0.37, 0.477, 0.37, 0.07))
    for index, number in enumerate((12, 16, 21)):
        size_ax.scatter(index * 1.3, 0, s=35 * number,
                        color="#c7c7c7", edgecolor="black", linewidth=0.45)
        size_ax.text(index * 1.3 + 0.42, 0, str(number), va="center", fontsize=7.5)
    size_ax.set(xlim=(-0.5, 3.6), ylim=(-1, 1))
    size_ax.axis("off")

    plot_domain_bars(
        left, counts, "proteomics",
        "Directionally concordant proteins"
    )
    plot_domain_bars(
        right, counts, "chemistries",
        "Directionally concordant chemistries"
    )
    left.set_ylabel("Directionally concordant associations (n)", fontsize=7.2)
    right.tick_params(axis="y", left=False, labelleft=False)
    for ax in (left, right):
        ax.title.set_fontsize(7.5)
        ax.title.set_fontweight("bold")
        ax.title.set_ha("left")
        ax.title.set_position((0, 1.02))
    legend = [
        Line2D([0], [0], marker="s", linestyle="", color="none",
               markerfacecolor=DOMAIN_COLORS[d], markeredgecolor="none",
               markersize=7, label=DOMAIN_LABELS[d])
        for d in DOMAINS
    ]
    fig.legend(handles=legend, title="Disease domain", loc="lower center",
               bbox_to_anchor=(0.54, 0.0), ncol=4, frameon=False,
               fontsize=6.4, title_fontsize=7.2, handlelength=1.2)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination.with_suffix(".pdf"), bbox_inches=None)
    fig.savefig(destination.with_suffix(".png"), dpi=300, bbox_inches=None)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data/aou/replication_summary.tsv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "figures")
    options = parser.parse_args()
    coverage, counts = summarize(options.data)
    prefix = options.output_dir / "Extended_Data_Fig4_AoU_directional_concordance"
    draw(coverage, counts, prefix)
    coverage.to_csv(prefix.with_name(prefix.name + "_coverage.tsv"), sep="\t", index=False)
    counts.to_csv(prefix.with_name(prefix.name + "_domains.tsv"), sep="\t", index=False)
    print(coverage.to_string(index=False))
    print(f"Saved {prefix.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
