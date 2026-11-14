"""Reusable utilities for the combined VAE + univariate AoU replication workflow.

The module contains only project-generic logic.  BigQuery execution remains in the
notebooks so that SQL is visible and easy to audit.  The functions here cover:

* discovery-target harmonization;
* feature and PRS/domain normalization;
* disease-registry SQL generation;
* ancestry-PC parsing;
* proteomics target mapping;
* chemistry feature assembly;
* a common HC3-robust OLS engine for both omic types;
* Benjamini-Hochberg correction on unique external-replication tests.
"""

from __future__ import annotations

import ast
import hashlib
import json
import math
import re
import unicodedata
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

# Optional analysis dependencies are imported lazily by the regression functions.


# -----------------------------------------------------------------------------
# Labels and shared constants
# -----------------------------------------------------------------------------

CASE_CONTROL_DOMAINS = [
    "AD_dementia",
    "Parkinsons",
    "ALS_MND",
    "CAD_coronary",
    "Heart_failure",
    "T2D",
    "CVD_composite",
]

DOMAIN_LABELS = {
    "AD_dementia": "Alzheimer's disease / dementia",
    "Parkinsons": "Parkinson's disease",
    "ALS_MND": "ALS / motor-neuron disease",
    "CAD_coronary": "Coronary artery disease",
    "Heart_failure": "Heart failure",
    "T2D": "Type 2 diabetes",
    "CVD_composite": "Cardiovascular disease composite",
    "Longevity_healthspan": "Longevity / healthspan",
    "unmapped": "Unmapped",
}

DISCOVERY_MODEL_LABELS = {
    "VAE_age": "VAE: age-included",
    "VAE_resid": "VAE: age-residualized",
    "OLS_age": "Univariate OLS: age-included",
    "OLS_resid": "Univariate OLS: age-residualized",
    "Spline": "Univariate spline age-interaction",
}

OMIC_LABELS = {
    "proteomics": "Proteomics",
    "chemistries": "Clinical chemistries",
}


# -----------------------------------------------------------------------------
# Basic normalization helpers
# -----------------------------------------------------------------------------


def _is_missing(value: Any) -> bool:
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def normalize_text(value: Any) -> str:
    """Return a conservative lowercase comparison key."""

    if _is_missing(value):
        return ""
    text = unicodedata.normalize("NFKC", str(value)).strip().lower()
    text = text.replace("≥", ">=").replace("≤", "<=")
    text = re.sub(r"[\s\-_]+", " ", text)
    text = re.sub(r"[^a-z0-9<>+=./ ]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def slugify(value: Any) -> str:
    """Return a stable ASCII identifier suitable for column names."""

    text = normalize_text(value)
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return re.sub(r"_+", "_", text).strip("_")


def normalize_gene(value: Any) -> str:
    if _is_missing(value):
        return ""
    text = str(value).strip().upper()
    text = re.sub(r"\s+", "", text)
    return text


def normalize_uniprot(value: Any) -> str:
    """Normalize a UniProt accession, removing an optional isoform suffix."""

    if _is_missing(value):
        return ""
    text = str(value).strip().upper()
    # The Olink file normally has one accession.  If several are present, retain
    # the first canonical-looking token for deterministic matching.
    tokens = re.split(r"[;,|\s]+", text)
    token = next((t for t in tokens if t), "")
    return re.sub(r"-\d+$", "", token)


def normalize_assay_symbol(value: Any) -> str:
    if _is_missing(value):
        return ""
    text = str(value).strip().upper()
    text = re.sub(r"[^A-Z0-9]+", "_", text)
    return re.sub(r"_+", "_", text).strip("_")


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def first_existing_column(
    df: pd.DataFrame,
    candidates: Sequence[str],
    required: bool = False,
) -> Optional[str]:
    for column in candidates:
        if column in df.columns:
            return column
    if required:
        raise ValueError(
            "None of the expected columns was found: " + ", ".join(candidates)
            + ". Available columns: " + ", ".join(map(str, df.columns))
        )
    return None


def stable_id(prefix: str, values: Iterable[Any], length: int = 16) -> str:
    payload = "||".join("" if _is_missing(v) else str(v) for v in values)
    digest = hashlib.sha1(payload.encode("utf-8")).hexdigest()[:length]
    return f"{prefix}_{digest}"


def join_unique(values: Iterable[Any], sep: str = "; ") -> str:
    cleaned = sorted(
        {
            str(v).strip()
            for v in values
            if not _is_missing(v) and str(v).strip() not in {"", "nan", "None"}
        }
    )
    return sep.join(cleaned)


def consensus_direction(values: Iterable[Any]) -> str:
    directions = {
        canonicalize_direction(v)
        for v in values
        if canonicalize_direction(v) in {"positive", "negative"}
    }
    if len(directions) == 1:
        return next(iter(directions))
    if len(directions) > 1:
        return "mixed"
    return "not_defined"


def canonicalize_direction(value: Any) -> str:
    if _is_missing(value):
        return "not_defined"
    text = normalize_text(value)
    if text in {"positive", "pos", "+", "higher", "increase", "increased"}:
        return "positive"
    if text in {"negative", "neg", "-", "lower", "decrease", "decreased"}:
        return "negative"
    if text in {"mixed", "conflicting", "mixed or unknown", "mixed_or_unknown"}:
        return "mixed"
    return "not_defined"


def direction_from_effect(value: Any) -> str:
    try:
        effect = float(value)
    except (TypeError, ValueError):
        return "not_defined"
    if not np.isfinite(effect):
        return "not_defined"
    if effect > 0:
        return "positive"
    if effect < 0:
        return "negative"
    return "zero"


# -----------------------------------------------------------------------------
# PRS-to-clinical-domain mapping
# -----------------------------------------------------------------------------


def map_prs_to_domain(prs_name: Any) -> str:
    """Map a discovery PRS label to the AoU case-control phenotype.

    The order matters: coronary and heart-failure labels are resolved before the
    broader cardiovascular pattern.
    """

    s = normalize_text(prs_name)

    if "alzheimer" in s or re.search(r"\bad\b", s):
        return "AD_dementia"
    if "parkinson" in s:
        return "Parkinsons"
    if (
        "amyotrophic lateral sclerosis" in s
        or "motor neuron" in s
        or re.search(r"(^| )als($| )", s)
    ):
        return "ALS_MND"
    if "heart failure" in s:
        return "Heart_failure"
    if "coronary atherosclerosis" in s or "coronary artery disease" in s:
        return "CAD_coronary"
    if "type 2 diabetes" in s or re.search(r"(^| )t2d($| )", s):
        return "T2D"
    if (
        "coronary vascular disease" in s
        or "cardiovascular disease" in s
        or re.search(r"(^| )cvd($| )", s)
    ):
        return "CVD_composite"
    if "timmers" in s or "longevity" in s or "healthspan" in s or "mvage" in s:
        return "Longevity_healthspan"

    return "unmapped"


# -----------------------------------------------------------------------------
# Chemistry feature specification and aliases
# -----------------------------------------------------------------------------


def load_chemistry_feature_spec(path: Path) -> pd.DataFrame:
    spec = pd.read_csv(path, sep="\t", dtype="string", keep_default_na=True)
    required = {
        "canonical_feature",
        "display_label",
        "aliases",
        "source_type",
        "source_feature_1",
        "component_1",
        "derivation",
        "search_terms",
    }
    missing = sorted(required - set(spec.columns))
    if missing:
        raise ValueError(f"Chemistry feature specification is missing: {missing}")

    spec = spec.copy()
    spec["canonical_feature"] = spec["canonical_feature"].map(slugify)
    if spec["canonical_feature"].duplicated().any():
        dupes = spec.loc[
            spec["canonical_feature"].duplicated(keep=False), "canonical_feature"
        ].tolist()
        raise ValueError(f"Duplicate canonical chemistry features: {dupes}")
    return spec


def chemistry_alias_lookup(spec: pd.DataFrame) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for row in spec.itertuples(index=False):
        canonical = slugify(row.canonical_feature)
        candidates = [canonical, getattr(row, "display_label", "")]
        for field in ["aliases", "legacy_feature"]:
            value = getattr(row, field, "")
            if not _is_missing(value):
                candidates.extend(str(value).split("|"))
        for candidate in candidates:
            key = normalize_text(candidate)
            if not key:
                continue
            existing = lookup.get(key)
            if existing is not None and existing != canonical:
                raise ValueError(
                    f"Chemistry alias {candidate!r} maps to both {existing!r} "
                    f"and {canonical!r}."
                )
            lookup[key] = canonical
    return lookup


def canonicalize_chemistry_feature(
    value: Any,
    alias_lookup: Mapping[str, str],
) -> str:
    key = normalize_text(value)
    if not key:
        return ""
    return alias_lookup.get(key, slugify(value))


# -----------------------------------------------------------------------------
# Discovery-target harmonization
# -----------------------------------------------------------------------------


def canonical_discovery_model(value: Any, framework: str) -> str:
    s = normalize_text(value)
    if framework == "VAE":
        if "resid" in s or "no age" in s:
            return "VAE_resid"
        if "age" in s or "included" in s:
            return "VAE_age"
        # VAE source tables usually have only these two models; keep an explicit
        # label rather than silently assigning an unknown row.
        return "VAE_unknown"

    if "spline" in s:
        return "Spline"
    if "resid" in s or "no age" in s:
        return "OLS_resid"
    if "age" in s or "included" in s:
        return "OLS_age"
    return "Univariate_unknown"


def _empty_like_index(index: pd.Index, dtype: str = "string") -> pd.Series:
    return pd.Series(pd.NA, index=index, dtype=dtype)


def harmonize_linear_targets(
    path: Path,
    omic_type: str,
    chemistry_spec: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Convert the supplied non-VAE OLS/spline target table to common schema."""

    raw = pd.read_csv(path, sep="\t", keep_default_na=True)
    required = {
        "model",
        "omic_name",
        "prs_name",
        "beta",
        "effect_direction",
        "fdr_q",
        "p_val",
        "omic_type",
    }
    if omic_type == "proteomics":
        required |= {"uniprot", "gene_name"}
    missing = sorted(required - set(raw.columns))
    if missing:
        raise ValueError(f"{path.name} is missing required columns: {missing}")

    out = pd.DataFrame(index=raw.index)
    out["discovery_framework"] = "Univariate"
    out["discovery_method"] = raw["model"].map(
        lambda x: "Spline" if "spline" in normalize_text(x) else "OLS"
    )
    out["discovery_model"] = raw["model"].map(
        lambda x: canonical_discovery_model(x, "Univariate")
    )
    out["discovery_model_label"] = out["discovery_model"].map(
        DISCOVERY_MODEL_LABELS
    ).fillna(out["discovery_model"])
    out["omic_type"] = omic_type
    out["omic_label"] = OMIC_LABELS[omic_type]
    out["feature_name_original"] = raw["omic_name"].astype("string")
    out["feature_display"] = raw["omic_name"].astype("string")
    out["prs_name"] = raw["prs_name"].astype("string")
    out["analysis_domain"] = out["prs_name"].map(map_prs_to_domain)
    out["analysis_domain_label"] = out["analysis_domain"].map(DOMAIN_LABELS)
    out["discovery_effect"] = numeric(raw["beta"])
    supplied_direction = raw["effect_direction"].map(canonicalize_direction)
    derived_direction = out["discovery_effect"].map(direction_from_effect)
    out["discovery_direction"] = supplied_direction.where(
        supplied_direction.isin(["positive", "negative"]), derived_direction
    )
    # Omnibus spline tests do not have one signed effect in this input.
    out.loc[out["discovery_model"].eq("Spline"), "discovery_direction"] = (
        "not_defined"
    )
    out["direction_evaluable"] = out["discovery_direction"].isin(
        ["positive", "negative"]
    )
    out["discovery_q"] = numeric(raw["fdr_q"])
    out["discovery_p"] = numeric(raw["p_val"])
    out["n_runs"] = np.nan
    out["percent_reproducibility"] = np.nan
    out["source_file"] = str(Path(path).resolve())
    out["source_candidate_id"] = (
        raw["candidate_id"].astype("string")
        if "candidate_id" in raw.columns
        else _empty_like_index(raw.index)
    )

    if omic_type == "proteomics":
        out["uniprot"] = raw["uniprot"].map(normalize_uniprot).replace("", pd.NA)
        out["gene_name"] = raw["gene_name"].map(normalize_gene).replace("", pd.NA)
        out["feature_key"] = np.where(
            out["uniprot"].notna(),
            "uniprot:" + out["uniprot"].astype(str),
            "gene:" + out["gene_name"].astype(str),
        )
        out["protein_identifier_crosswalk"] = np.where(
            out["uniprot"].notna(), "supplied_uniprot", "supplied_gene_only"
        )
        out["chemistry_canonical_feature"] = pd.NA
    else:
        if chemistry_spec is None:
            raise ValueError("chemistry_spec is required for chemistry targets")
        aliases = chemistry_alias_lookup(chemistry_spec)
        out["chemistry_canonical_feature"] = raw["omic_name"].map(
            lambda x: canonicalize_chemistry_feature(x, aliases)
        )
        display_map = chemistry_spec.set_index("canonical_feature")[
            "display_label"
        ].to_dict()
        out["feature_display"] = out["chemistry_canonical_feature"].map(
            display_map
        ).fillna(out["feature_name_original"])
        out["feature_key"] = "chem:" + out[
            "chemistry_canonical_feature"
        ].astype(str)
        out["uniprot"] = pd.NA
        out["gene_name"] = pd.NA
        out["protein_identifier_crosswalk"] = pd.NA

    return finalize_target_manifest(out)


def harmonize_vae_targets(
    path: Path,
    omic_type: str,
    chemistry_spec: Optional[pd.DataFrame] = None,
    min_n_runs: int = 6,
) -> pd.DataFrame:
    """Convert a legacy VAE target file to the common target schema.

    The adapter accepts the column variants observed in the supplied notebooks.
    """

    raw = pd.read_csv(path, sep="\t", keep_default_na=True)

    model_col = first_existing_column(raw, ["model", "model_code"], required=True)
    prs_col = first_existing_column(
        raw,
        ["disease_or_prs", "prs_name", "source_prs", "prs_base"],
        required=True,
    )
    feature_col = first_existing_column(
        raw,
        ["feature", "feature_lower", "omic_name", "omic_label", "omic_raw"],
        required=True,
    )
    effect_col = first_existing_column(
        raw,
        ["avg_bayes_k", "mean_bayes_k", "bayes_k", "effect_value", "beta"],
        required=True,
    )
    direction_col = first_existing_column(
        raw,
        ["direction", "vae_direction_from_bayes_k", "effect_direction"],
    )
    q_col = first_existing_column(raw, ["fdr", "fdr_q", "q_value", "q_val"])
    p_col = first_existing_column(raw, ["p_val", "p_value", "pvalue"])
    domain_col = first_existing_column(raw, ["analysis_domain", "domain"])
    n_runs_col = first_existing_column(raw, ["n_runs"])
    reproducibility_col = first_existing_column(
        raw, ["percent_reproducibility", "reproducibility"]
    )
    gene_col = first_existing_column(raw, ["gene_name", "gene", "assay"])
    uniprot_col = first_existing_column(raw, ["uniprot", "uniprot_id"])

    if n_runs_col is not None:
        raw[n_runs_col] = numeric(raw[n_runs_col])
        raw = raw[raw[n_runs_col].isna() | raw[n_runs_col].ge(min_n_runs)].copy()

    # When a broad fallback table is used, isolate proteins/non-proteins here.
    if omic_type == "proteomics":
        if gene_col is None and uniprot_col is None:
            raise ValueError(
                f"{path.name} does not contain gene_name or uniprot for proteins."
            )
        protein_mask = pd.Series(False, index=raw.index)
        if gene_col is not None:
            protein_mask |= raw[gene_col].notna() & raw[gene_col].astype(str).str.strip().ne("")
        if uniprot_col is not None:
            protein_mask |= raw[uniprot_col].notna() & raw[uniprot_col].astype(str).str.strip().ne("")
        raw = raw[protein_mask].copy()
    elif gene_col is not None or uniprot_col is not None:
        # Chemistry input files are normally already filtered.  If not, remove
        # clear protein rows while retaining rows with missing protein metadata.
        protein_mask = pd.Series(False, index=raw.index)
        if gene_col is not None:
            protein_mask |= raw[gene_col].notna() & raw[gene_col].astype(str).str.strip().ne("")
        if uniprot_col is not None:
            protein_mask |= raw[uniprot_col].notna() & raw[uniprot_col].astype(str).str.strip().ne("")
        raw = raw[~protein_mask].copy()

    out = pd.DataFrame(index=raw.index)
    out["discovery_framework"] = "VAE"
    out["discovery_method"] = "VAE"
    out["discovery_model"] = raw[model_col].map(
        lambda x: canonical_discovery_model(x, "VAE")
    )
    out["discovery_model_label"] = out["discovery_model"].map(
        DISCOVERY_MODEL_LABELS
    ).fillna(out["discovery_model"])
    out["omic_type"] = omic_type
    out["omic_label"] = OMIC_LABELS[omic_type]
    out["feature_name_original"] = raw[feature_col].astype("string")
    out["feature_display"] = raw[feature_col].astype("string")
    out["prs_name"] = raw[prs_col].astype("string")
    if domain_col is not None:
        supplied_domain = raw[domain_col].astype("string")
        out["analysis_domain"] = supplied_domain.where(
            supplied_domain.isin(list(DOMAIN_LABELS)),
            out["prs_name"].map(map_prs_to_domain),
        )
    else:
        out["analysis_domain"] = out["prs_name"].map(map_prs_to_domain)
    out["analysis_domain_label"] = out["analysis_domain"].map(DOMAIN_LABELS)
    out["discovery_effect"] = numeric(raw[effect_col])
    if direction_col is not None:
        supplied_direction = raw[direction_col].map(canonicalize_direction)
    else:
        supplied_direction = pd.Series("not_defined", index=raw.index)
    derived_direction = out["discovery_effect"].map(direction_from_effect)
    out["discovery_direction"] = supplied_direction.where(
        supplied_direction.isin(["positive", "negative"]), derived_direction
    )
    out["direction_evaluable"] = out["discovery_direction"].isin(
        ["positive", "negative"]
    )
    out["discovery_q"] = numeric(raw[q_col]) if q_col else np.nan
    out["discovery_p"] = numeric(raw[p_col]) if p_col else np.nan
    out["n_runs"] = numeric(raw[n_runs_col]) if n_runs_col else np.nan
    out["percent_reproducibility"] = (
        numeric(raw[reproducibility_col]) if reproducibility_col else np.nan
    )
    out["source_file"] = str(Path(path).resolve())
    out["source_candidate_id"] = _empty_like_index(raw.index)

    if omic_type == "proteomics":
        out["gene_name"] = (
            raw[gene_col].map(normalize_gene).replace("", pd.NA)
            if gene_col
            else pd.NA
        )
        out["uniprot"] = (
            raw[uniprot_col].map(normalize_uniprot).replace("", pd.NA)
            if uniprot_col
            else pd.NA
        )
        out["feature_key"] = np.where(
            out["uniprot"].notna(),
            "uniprot:" + out["uniprot"].astype(str),
            "gene:" + out["gene_name"].astype(str),
        )
        out["protein_identifier_crosswalk"] = np.where(
            out["uniprot"].notna(), "supplied_uniprot", "supplied_gene_only"
        )
        out["chemistry_canonical_feature"] = pd.NA
    else:
        if chemistry_spec is None:
            raise ValueError("chemistry_spec is required for chemistry targets")
        aliases = chemistry_alias_lookup(chemistry_spec)
        # Prefer feature_lower when available because it was curated in the VAE
        # workflow, but preserve the human-readable feature for display.
        canonical_source_col = (
            "feature_lower" if "feature_lower" in raw.columns else feature_col
        )
        out["chemistry_canonical_feature"] = raw[canonical_source_col].map(
            lambda x: canonicalize_chemistry_feature(x, aliases)
        )
        display_map = chemistry_spec.set_index("canonical_feature")[
            "display_label"
        ].to_dict()
        out["feature_display"] = out["chemistry_canonical_feature"].map(
            display_map
        ).fillna(out["feature_name_original"])
        out["feature_key"] = "chem:" + out[
            "chemistry_canonical_feature"
        ].astype(str)
        out["uniprot"] = pd.NA
        out["gene_name"] = pd.NA
        out["protein_identifier_crosswalk"] = pd.NA

    return finalize_target_manifest(out)


def apply_protein_gene_uniprot_crosswalk(
    targets: pd.DataFrame,
    crosswalk_source: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fill missing protein UniProt accessions from a one-to-one gene crosswalk.

    This is primarily useful when a legacy VAE target table identifies proteins
    only by gene symbol while a newer univariate table contains both gene and
    UniProt identifiers.  Only unambiguous gene-to-UniProt mappings are applied;
    genes mapping to more than one accession are retained unchanged and reported
    in the audit table.
    """

    out = targets.copy()
    source = crosswalk_source.copy()
    for table, label in [(out, "targets"), (source, "crosswalk_source")]:
        missing = {"gene_name", "uniprot", "omic_type"} - set(table.columns)
        if missing:
            raise ValueError(f"{label} is missing required columns: {sorted(missing)}")

    source = source[source["omic_type"].eq("proteomics")].copy()
    source["gene_key"] = source["gene_name"].map(normalize_gene)
    source["uniprot_key"] = source["uniprot"].map(normalize_uniprot)
    source = source[source["gene_key"].ne("") & source["uniprot_key"].ne("")].copy()

    mapping_rows: list[dict[str, Any]] = []
    for gene, group in source.groupby("gene_key", sort=True):
        accessions = sorted(set(group["uniprot_key"]))
        mapping_rows.append(
            {
                "gene_name": gene,
                "n_uniprot_accessions": len(accessions),
                "uniprot_accessions": "; ".join(accessions),
                "uniprot_to_apply": accessions[0] if len(accessions) == 1 else pd.NA,
                "crosswalk_status": (
                    "one_to_one" if len(accessions) == 1 else "ambiguous_gene_to_uniprot"
                ),
            }
        )
    audit = pd.DataFrame(mapping_rows)

    if audit.empty:
        out["protein_identifier_crosswalk"] = "not_applied"
        return finalize_target_manifest(out), audit

    one_to_one = audit[audit["crosswalk_status"].eq("one_to_one")].set_index(
        "gene_name"
    )["uniprot_to_apply"].to_dict()

    out["protein_identifier_crosswalk"] = "not_needed"
    protein = out["omic_type"].eq("proteomics")
    missing_uniprot = protein & out["uniprot"].map(normalize_uniprot).eq("")
    gene_key = out["gene_name"].map(normalize_gene)
    mapped = gene_key.map(one_to_one)
    apply_mask = missing_uniprot & mapped.notna()
    ambiguous_genes = set(
        audit.loc[
            audit["crosswalk_status"].eq("ambiguous_gene_to_uniprot"), "gene_name"
        ]
    )
    ambiguous_mask = missing_uniprot & gene_key.isin(ambiguous_genes)

    out.loc[apply_mask, "uniprot"] = mapped.loc[apply_mask].values
    out.loc[apply_mask, "feature_key"] = (
        "uniprot:" + mapped.loc[apply_mask].astype(str)
    ).values
    out.loc[apply_mask, "protein_identifier_crosswalk"] = (
        "filled_from_univariate_gene_crosswalk"
    )
    out.loc[ambiguous_mask, "protein_identifier_crosswalk"] = (
        "not_applied_ambiguous_gene"
    )
    out.loc[missing_uniprot & ~apply_mask & ~ambiguous_mask, "protein_identifier_crosswalk"] = (
        "no_crosswalk_match"
    )

    target_missing_counts = gene_key.loc[missing_uniprot].value_counts()
    target_updated_counts = gene_key.loc[apply_mask].value_counts()
    audit["n_target_rows_missing_uniprot"] = (
        audit["gene_name"].map(target_missing_counts).fillna(0).astype(int)
    )
    audit["n_target_rows_updated"] = (
        audit["gene_name"].map(target_updated_counts).fillna(0).astype(int)
    )

    # Recompute stable target IDs because feature_key may have changed.  The
    # standard output schema intentionally does not retain the temporary audit
    # column; the returned audit table records the crosswalk itself.
    return finalize_target_manifest(out), audit


def finalize_target_manifest(out: pd.DataFrame) -> pd.DataFrame:
    out = out.copy()
    out["analysis_domain"] = out["analysis_domain"].fillna("unmapped")
    out["analysis_domain_label"] = out["analysis_domain"].map(DOMAIN_LABELS).fillna(
        out["analysis_domain"]
    )

    reasons = pd.Series("", index=out.index, dtype="string")
    reasons.loc[out["analysis_domain"].eq("Longevity_healthspan")] = (
        "No prespecified AoU disease case-control phenotype for longevity/healthspan"
    )
    reasons.loc[out["analysis_domain"].eq("unmapped")] = "PRS could not be mapped to an AoU disease domain"
    reasons.loc[~out["analysis_domain"].isin(CASE_CONTROL_DOMAINS) & reasons.eq("")] = (
        "Analysis domain is not in the case-control registry"
    )
    reasons.loc[out["feature_key"].astype(str).isin(["", "nan", "gene:<NA>", "chem:"])] = (
        "Feature identifier is missing"
    )

    out["not_testable_reason"] = reasons.replace("", pd.NA)
    out["clinically_testable"] = out["not_testable_reason"].isna()

    ids = []
    for row in out.itertuples(index=False):
        ids.append(
            stable_id(
                "target",
                [
                    row.discovery_framework,
                    row.discovery_model,
                    row.omic_type,
                    row.analysis_domain,
                    row.feature_key,
                    row.prs_name,
                    row.source_candidate_id,
                ],
            )
        )
    out["target_id"] = ids

    ordered = [
        "target_id",
        "discovery_framework",
        "discovery_method",
        "discovery_model",
        "discovery_model_label",
        "omic_type",
        "omic_label",
        "analysis_domain",
        "analysis_domain_label",
        "prs_name",
        "feature_name_original",
        "feature_display",
        "feature_key",
        "chemistry_canonical_feature",
        "uniprot",
        "gene_name",
        "protein_identifier_crosswalk",
        "discovery_effect",
        "discovery_direction",
        "direction_evaluable",
        "discovery_p",
        "discovery_q",
        "n_runs",
        "percent_reproducibility",
        "clinically_testable",
        "not_testable_reason",
        "source_candidate_id",
        "source_file",
    ]
    for column in ordered:
        if column not in out.columns:
            out[column] = pd.NA
    out = out[ordered].drop_duplicates().reset_index(drop=True)

    # A hash collision is extraordinarily unlikely, but duplicate target IDs more
    # often signal truly duplicated source rows and should be made explicit.
    if out["target_id"].duplicated().any():
        duplicate_mask = out["target_id"].duplicated(keep=False)
        # Preserve exact duplicates only once; fail if the payload differs.
        distinct = out.loc[duplicate_mask].drop(columns="target_id").drop_duplicates()
        if len(distinct) != out.loc[duplicate_mask, "target_id"].nunique():
            raise ValueError("Non-identical target rows produced duplicate target_id values")
        out = out.drop_duplicates(subset="target_id").reset_index(drop=True)

    return out


def group_targets_for_analysis(targets: pd.DataFrame) -> pd.DataFrame:
    """Collapse repeated PRS rows that imply the same external AoU test.

    Grouping retains discovery framework/model so VAE and univariate evidence are
    compared without duplicating identical rows within a model.
    """

    group_cols = [
        "discovery_framework",
        "discovery_method",
        "discovery_model",
        "discovery_model_label",
        "omic_type",
        "omic_label",
        "analysis_domain",
        "analysis_domain_label",
        "feature_key",
        "feature_display",
        "chemistry_canonical_feature",
        "uniprot",
        "gene_name",
        "clinically_testable",
        "not_testable_reason",
    ]

    rows: list[dict[str, Any]] = []
    for keys, group in targets.groupby(group_cols, dropna=False, sort=False):
        record = dict(zip(group_cols, keys))
        record.update(
            {
                "analysis_target_id": stable_id(
                    "analysis_target",
                    [
                        record["discovery_framework"],
                        record["discovery_model"],
                        record["omic_type"],
                        record["analysis_domain"],
                        record["feature_key"],
                    ],
                ),
                "source_target_ids": join_unique(group["target_id"]),
                "prs_names": join_unique(group["prs_name"]),
                "n_discovery_rows": int(len(group)),
                "discovery_direction": consensus_direction(
                    group["discovery_direction"]
                ),
                "direction_evaluable": consensus_direction(
                    group["discovery_direction"]
                )
                in {"positive", "negative"},
                "min_discovery_p": numeric(group["discovery_p"]).min(skipna=True),
                "min_discovery_q": numeric(group["discovery_q"]).min(skipna=True),
                "mean_discovery_effect": numeric(group["discovery_effect"]).mean(
                    skipna=True
                ),
                "max_abs_discovery_effect": numeric(group["discovery_effect"])
                .abs()
                .max(skipna=True),
                "min_n_runs": numeric(group["n_runs"]).min(skipna=True),
                "mean_percent_reproducibility": numeric(
                    group["percent_reproducibility"]
                ).mean(skipna=True),
                "source_files": join_unique(group["source_file"]),
            }
        )
        rows.append(record)

    grouped = pd.DataFrame(rows)
    if grouped.empty:
        return grouped
    return grouped.sort_values(
        [
            "clinically_testable",
            "omic_type",
            "analysis_domain",
            "discovery_framework",
            "discovery_model",
            "min_discovery_q",
        ],
        ascending=[False, True, True, True, True, True],
        na_position="last",
    ).reset_index(drop=True)


def build_target_overlap_audit(targets: pd.DataFrame) -> pd.DataFrame:
    """Summarize whether a domain-feature pair occurs in VAE and univariate sets."""

    key_cols = ["omic_type", "analysis_domain", "feature_key", "feature_display"]
    rows = []
    for keys, group in targets.groupby(key_cols, dropna=False):
        frameworks = set(group["discovery_framework"].dropna().astype(str))
        rows.append(
            {
                **dict(zip(key_cols, keys)),
                "in_vae": "VAE" in frameworks,
                "in_univariate": "Univariate" in frameworks,
                "overlap_class": (
                    "VAE and univariate"
                    if frameworks == {"VAE", "Univariate"}
                    else "VAE only"
                    if "VAE" in frameworks
                    else "Univariate only"
                ),
                "discovery_models": join_unique(group["discovery_model_label"]),
                "prs_names": join_unique(group["prs_name"]),
                "n_target_rows": len(group),
            }
        )
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Disease definition and BigQuery SQL helpers
# -----------------------------------------------------------------------------


def load_disease_definitions(path: Path) -> list[dict[str, str]]:
    with open(path, "r", encoding="utf-8") as handle:
        definitions = json.load(handle)
    required = {"domain", "label", "case_regex", "control_exclude_regex", "seed_exclude_regex"}
    for i, definition in enumerate(definitions):
        missing = required - set(definition)
        if missing:
            raise ValueError(f"Disease definition {i} is missing {sorted(missing)}")
    return definitions


def _bq_string(value: str) -> str:
    return "'" + str(value).replace("\\", "\\\\").replace("'", "\\'") + "'"


def build_disease_patterns_cte(definitions: Sequence[Mapping[str, str]]) -> str:
    rows = []
    for definition in definitions:
        for set_type, regex_key in [
            ("case", "case_regex"),
            ("control_exclude", "control_exclude_regex"),
        ]:
            rows.append(
                "SELECT "
                f"{_bq_string(definition['domain'])} AS disease_domain, "
                f"{_bq_string(set_type)} AS set_type, "
                f"r{_bq_string(definition[regex_key])} AS include_regex, "
                f"r{_bq_string(definition['seed_exclude_regex'])} AS exclude_regex"
            )
    return "\nUNION ALL\n".join(rows)


def disease_concept_registry_query(cdr: str, patterns_cte: str) -> str:
    return f"""
WITH disease_patterns AS (
  {patterns_cte}
),
seed_concepts AS (
  SELECT DISTINCT
    p.disease_domain,
    p.set_type,
    c.concept_id AS seed_concept_id,
    c.concept_name AS seed_concept_name,
    p.exclude_regex
  FROM `{cdr}.concept` c
  CROSS JOIN disease_patterns p
  WHERE c.domain_id = 'Condition'
    AND c.standard_concept = 'S'
    AND c.invalid_reason IS NULL
    AND REGEXP_CONTAINS(LOWER(c.concept_name), p.include_regex)
    AND NOT REGEXP_CONTAINS(LOWER(c.concept_name), p.exclude_regex)
),
expanded_concepts AS (
  SELECT DISTINCT
    s.disease_domain,
    s.set_type,
    s.seed_concept_id,
    s.seed_concept_name,
    c.concept_id,
    c.concept_name,
    c.concept_code,
    c.vocabulary_id,
    c.concept_class_id
  FROM seed_concepts s
  JOIN `{cdr}.concept_ancestor` ca
    ON s.seed_concept_id = ca.ancestor_concept_id
  JOIN `{cdr}.concept` c
    ON ca.descendant_concept_id = c.concept_id
  WHERE c.domain_id = 'Condition'
    AND c.standard_concept = 'S'
    AND c.invalid_reason IS NULL
    AND NOT REGEXP_CONTAINS(LOWER(c.concept_name), s.exclude_regex)

  UNION DISTINCT

  SELECT DISTINCT
    s.disease_domain,
    s.set_type,
    s.seed_concept_id,
    s.seed_concept_name,
    c.concept_id,
    c.concept_name,
    c.concept_code,
    c.vocabulary_id,
    c.concept_class_id
  FROM seed_concepts s
  JOIN `{cdr}.concept` c
    ON s.seed_concept_id = c.concept_id
)
SELECT
  ec.*,
  COUNT(co.condition_occurrence_id) AS n_records_all_aou,
  COUNT(DISTINCT co.person_id) AS n_people_all_aou
FROM expanded_concepts ec
LEFT JOIN `{cdr}.condition_occurrence` co
  ON co.condition_concept_id = ec.concept_id
GROUP BY
  disease_domain, set_type, seed_concept_id, seed_concept_name,
  concept_id, concept_name, concept_code, vocabulary_id, concept_class_id
ORDER BY disease_domain, set_type, n_people_all_aou DESC
""".strip()


def disease_flag_hits_query(
    cdr: str,
    patterns_cte: str,
    definitions: Sequence[Mapping[str, str]],
) -> str:
    domains = [d["domain"] for d in definitions]
    expressions = []
    for domain in domains:
        expressions.extend(
            [
                f"MAX(IF(disease_domain = '{domain}' AND set_type = 'case', 1, 0)) AS {domain}_case",
                f"MAX(IF(disease_domain = '{domain}' AND set_type = 'control_exclude', 1, 0)) AS {domain}_control_exclude",
                f"MIN(IF(disease_domain = '{domain}' AND set_type = 'case', first_condition_date, NULL)) AS {domain}_first_case_date",
            ]
        )
    wide_sql = ",\n  ".join(expressions)

    return f"""
WITH disease_patterns AS (
  {patterns_cte}
),
seed_concepts AS (
  SELECT DISTINCT
    p.disease_domain,
    p.set_type,
    c.concept_id AS seed_concept_id,
    p.exclude_regex
  FROM `{cdr}.concept` c
  CROSS JOIN disease_patterns p
  WHERE c.domain_id = 'Condition'
    AND c.standard_concept = 'S'
    AND c.invalid_reason IS NULL
    AND REGEXP_CONTAINS(LOWER(c.concept_name), p.include_regex)
    AND NOT REGEXP_CONTAINS(LOWER(c.concept_name), p.exclude_regex)
),
expanded_concepts AS (
  SELECT DISTINCT
    s.disease_domain,
    s.set_type,
    c.concept_id
  FROM seed_concepts s
  JOIN `{cdr}.concept_ancestor` ca
    ON s.seed_concept_id = ca.ancestor_concept_id
  JOIN `{cdr}.concept` c
    ON ca.descendant_concept_id = c.concept_id
  WHERE c.domain_id = 'Condition'
    AND c.standard_concept = 'S'
    AND c.invalid_reason IS NULL
    AND NOT REGEXP_CONTAINS(LOWER(c.concept_name), s.exclude_regex)

  UNION DISTINCT

  SELECT disease_domain, set_type, seed_concept_id AS concept_id
  FROM seed_concepts
),
condition_hits AS (
  SELECT
    co.person_id,
    ec.disease_domain,
    ec.set_type,
    MIN(co.condition_start_date) AS first_condition_date
  FROM `{cdr}.condition_occurrence` co
  JOIN expanded_concepts ec
    ON co.condition_concept_id = ec.concept_id
  GROUP BY co.person_id, ec.disease_domain, ec.set_type
)
SELECT
  person_id,
  {wide_sql}
FROM condition_hits
GROUP BY person_id
""".strip()


def add_case_control_flags(
    cohort: pd.DataFrame,
    flag_hits: pd.DataFrame,
    base_domains: Sequence[str],
) -> pd.DataFrame:
    """Merge AoU-wide condition hits and derive domain/clean control flags."""

    if cohort["person_id"].duplicated().any():
        raise ValueError("cohort must contain one row per person_id")
    flags = cohort[["person_id"]].merge(
        flag_hits,
        on="person_id",
        how="left",
        validate="one_to_one",
    )
    base_flag_cols = []
    for domain in base_domains:
        base_flag_cols.extend([f"{domain}_case", f"{domain}_control_exclude"])
        for col in [f"{domain}_case", f"{domain}_control_exclude"]:
            if col not in flags.columns:
                flags[col] = 0
    flags[base_flag_cols] = flags[base_flag_cols].fillna(0).astype(int)

    for domain in base_domains:
        flags[f"{domain}_domain_control_candidate"] = (
            flags[f"{domain}_case"].eq(0)
            & flags[f"{domain}_control_exclude"].eq(0)
        ).astype(int)

    exclude_cols = [f"{d}_control_exclude" for d in base_domains]
    flags["clean_control_candidate"] = flags[exclude_cols].sum(axis=1).eq(0).astype(int)

    # The broad CVD phenotype is intentionally derived from CAD/coronary and HF.
    flags["CVD_composite_case"] = (
        flags["CAD_coronary_case"].eq(1) | flags["Heart_failure_case"].eq(1)
    ).astype(int)
    flags["CVD_composite_control_exclude"] = (
        flags["CAD_coronary_control_exclude"].eq(1)
        | flags["Heart_failure_control_exclude"].eq(1)
    ).astype(int)
    flags["CVD_composite_domain_control_candidate"] = (
        flags["CVD_composite_case"].eq(0)
        & flags["CVD_composite_control_exclude"].eq(0)
    ).astype(int)

    # Retain first diagnosis dates when available. They are not used in the
    # primary cross-sectional model, but keeping them makes later temporal
    # sensitivity analyses possible without repeating the disease extraction.
    date_columns = []
    for domain in base_domains:
        column = f"{domain}_first_case_date"
        if column in flags.columns:
            flags[column] = pd.to_datetime(flags[column], errors="coerce")
            date_columns.append(column)
    cad_date = "CAD_coronary_first_case_date"
    hf_date = "Heart_failure_first_case_date"
    if cad_date in flags.columns or hf_date in flags.columns:
        date_frame = pd.DataFrame(index=flags.index)
        date_frame[cad_date] = (
            flags[cad_date] if cad_date in flags.columns else pd.NaT
        )
        date_frame[hf_date] = (
            flags[hf_date] if hf_date in flags.columns else pd.NaT
        )
        flags["CVD_composite_first_case_date"] = date_frame.min(axis=1)

    return flags


# -----------------------------------------------------------------------------
# Genetic ancestry PCs
# -----------------------------------------------------------------------------


def parse_ancestry_pcs(
    path: Path,
    n_pcs: int = 8,
) -> pd.DataFrame:
    """Parse AoU ancestry predictions and the stringified PCA vector."""

    raw = pd.read_csv(path, sep="\t", low_memory=False)
    person_col = first_existing_column(raw, ["person_id", "research_id"], required=True)
    if "pca_features" not in raw.columns:
        # Also accept a previously parsed table.
        expected = [f"genetic_pc{i}" for i in range(1, n_pcs + 1)]
        missing = [c for c in expected if c not in raw.columns]
        if missing:
            raise ValueError(f"Ancestry file lacks pca_features and {missing}")
        out_cols = [person_col]
        for optional in ["ancestry_pred", "probabilities", "ancestry_pred_other"]:
            if optional in raw.columns:
                out_cols.append(optional)
        out_cols += expected
        out = raw[out_cols].copy().rename(columns={person_col: "person_id"})
        out["person_id"] = pd.to_numeric(out["person_id"], errors="raise").astype("int64")
        return out.drop_duplicates("person_id")

    def parse_vector(value: Any) -> list[float]:
        if _is_missing(value):
            return []
        if isinstance(value, (list, tuple, np.ndarray)):
            return [float(x) for x in value]
        parsed = ast.literal_eval(str(value))
        return [float(x) for x in parsed]

    vectors = raw["pca_features"].map(parse_vector)
    lengths = vectors.map(len)
    if lengths.max() < n_pcs:
        raise ValueError(
            f"PCA vectors contain at most {lengths.max()} values; {n_pcs} required"
        )

    pc_data = pd.DataFrame(
        [v[:n_pcs] for v in vectors],
        columns=[f"genetic_pc{i}" for i in range(1, n_pcs + 1)],
        index=raw.index,
    )
    metadata_cols = [person_col]
    for optional in ["ancestry_pred", "probabilities", "ancestry_pred_other"]:
        if optional in raw.columns:
            metadata_cols.append(optional)
    out = pd.concat([raw[metadata_cols], pc_data], axis=1)
    out = out.rename(columns={person_col: "person_id"})
    out["person_id"] = pd.to_numeric(out["person_id"], errors="raise").astype("int64")
    if out["person_id"].duplicated().any():
        out = out.sort_values("person_id").drop_duplicates("person_id", keep="first")
    return out.reset_index(drop=True)


# -----------------------------------------------------------------------------
# Proteomics helpers
# -----------------------------------------------------------------------------


PROTEIN_GENE_ALIASES: dict[str, list[str]] = {
    "FIGF": ["VEGFD"],
    "VEGFD": ["FIGF"],
    "IL27,EBI3": ["IL27_EBI3", "IL27", "EBI3"],
    "IL27_EBI3": ["IL27", "EBI3"],
}


def gene_assay_candidates(gene_name: Any) -> list[str]:
    gene = normalize_gene(gene_name)
    if not gene:
        return []
    candidates = [gene, *PROTEIN_GENE_ALIASES.get(gene, [])]
    candidates += [gene.replace(",", "_"), gene.replace("-", "_")]
    normalized = []
    seen = set()
    for candidate in candidates:
        key = normalize_assay_symbol(candidate)
        if key and key not in seen:
            normalized.append(key)
            seen.add(key)
    return normalized


def scan_target_proteomics_rows(
    proteomics_file: Path,
    target_proteins: pd.DataFrame,
    value_preference: Sequence[str] = ("PCNormalizedNPX", "NPX"),
    chunksize: int = 1_000_000,
    progress_every: int = 10,
) -> tuple[pd.DataFrame, str]:
    """Read only target-relevant Olink rows in one pass through the large file."""

    header = pd.read_csv(proteomics_file, sep="\t", nrows=0).columns.tolist()
    value_col = next((c for c in value_preference if c in header), None)
    if value_col is None:
        raise ValueError(
            f"None of the requested Olink value columns {value_preference} exists. "
            f"Available: {header}"
        )

    required = ["ResearchID", "Assay", "UniProt", "AssayQC", "SampleQC", value_col]
    missing = [c for c in required if c not in header]
    if missing:
        raise ValueError(f"Proteomics file is missing columns: {missing}")

    target_uniprots = {
        normalize_uniprot(v)
        for v in target_proteins.get("uniprot", pd.Series(dtype=object))
        if normalize_uniprot(v)
    }
    target_assays = set()
    for gene in target_proteins.get("gene_name", pd.Series(dtype=object)):
        target_assays.update(gene_assay_candidates(gene))

    pieces: list[pd.DataFrame] = []
    reader = pd.read_csv(
        proteomics_file,
        sep="\t",
        usecols=required,
        chunksize=chunksize,
        low_memory=False,
        dtype={
            "ResearchID": "string",
            "Assay": "string",
            "UniProt": "string",
            "AssayQC": "string",
            "SampleQC": "string",
        },
    )

    for i, chunk in enumerate(reader, start=1):
        uniprot_key = chunk["UniProt"].map(normalize_uniprot)
        assay_key = chunk["Assay"].map(normalize_assay_symbol)
        mask = uniprot_key.isin(target_uniprots) | assay_key.isin(target_assays)
        if not mask.any():
            if progress_every and i % progress_every == 0:
                print(f"Processed {i} proteomics chunks; no new target rows in this chunk")
            continue

        sub = chunk.loc[mask].copy()
        sub["uniprot_key"] = uniprot_key.loc[mask].values
        sub["assay_key"] = assay_key.loc[mask].values
        sub = sub[
            sub["AssayQC"].fillna("").str.upper().eq("PASS")
            & sub["SampleQC"].fillna("").str.upper().eq("PASS")
        ].copy()
        sub["person_id"] = pd.to_numeric(sub["ResearchID"], errors="coerce")
        sub["protein_value"] = pd.to_numeric(sub[value_col], errors="coerce")
        sub = sub.dropna(subset=["person_id", "protein_value", "Assay"])
        sub["person_id"] = sub["person_id"].astype("int64")
        pieces.append(
            sub[
                [
                    "person_id",
                    "Assay",
                    "assay_key",
                    "UniProt",
                    "uniprot_key",
                    "protein_value",
                ]
            ]
        )
        if progress_every and i % progress_every == 0:
            print(f"Processed {i} proteomics chunks; retained {sum(len(x) for x in pieces):,} rows")

    if not pieces:
        return pd.DataFrame(
            columns=[
                "person_id",
                "Assay",
                "assay_key",
                "UniProt",
                "uniprot_key",
                "protein_value",
            ]
        ), value_col
    return pd.concat(pieces, ignore_index=True), value_col


def select_protein_assays(
    target_proteins: pd.DataFrame,
    olink_rows: pd.DataFrame,
) -> pd.DataFrame:
    """Choose one AoU Olink assay per target using exact UniProt first."""

    if olink_rows.empty:
        coverage = pd.DataFrame(
            columns=["Assay", "assay_key", "uniprot_key", "n_people_qc"]
        )
    else:
        coverage = (
            olink_rows.groupby(
                ["Assay", "assay_key", "uniprot_key"],
                dropna=False,
                as_index=False,
            )["person_id"]
            .nunique()
            .rename(columns={"person_id": "n_people_qc"})
        )

    unique_targets = target_proteins[
        ["feature_key", "feature_display", "uniprot", "gene_name"]
    ].drop_duplicates()
    records = []
    for target in unique_targets.itertuples(index=False):
        target_uniprot = normalize_uniprot(target.uniprot)
        gene_candidates = set(gene_assay_candidates(target.gene_name))

        exact_uniprot = coverage[
            coverage["uniprot_key"].astype(str).eq(target_uniprot)
        ] if target_uniprot else coverage.iloc[0:0]
        if not exact_uniprot.empty:
            candidates = exact_uniprot.copy()
            method = "exact_uniprot"
        else:
            candidates = coverage[coverage["assay_key"].isin(gene_candidates)].copy()
            method = "gene_or_alias" if not candidates.empty else "unmapped"

        if candidates.empty:
            records.append(
                {
                    "feature_key": target.feature_key,
                    "feature_display": target.feature_display,
                    "uniprot": target_uniprot or pd.NA,
                    "gene_name": normalize_gene(target.gene_name) or pd.NA,
                    "aou_assay": pd.NA,
                    "aou_assay_key": pd.NA,
                    "aou_uniprot": pd.NA,
                    "aou_n_people_qc": 0,
                    "protein_mapping_method": method,
                    "n_mapping_candidates": 0,
                    "alternative_assays": pd.NA,
                }
            )
            continue

        candidates = candidates.sort_values(
            ["n_people_qc", "Assay"], ascending=[False, True]
        )
        best = candidates.iloc[0]
        records.append(
            {
                "feature_key": target.feature_key,
                "feature_display": target.feature_display,
                "uniprot": target_uniprot or pd.NA,
                "gene_name": normalize_gene(target.gene_name) or pd.NA,
                "aou_assay": best["Assay"],
                "aou_assay_key": best["assay_key"],
                "aou_uniprot": best["uniprot_key"],
                "aou_n_people_qc": int(best["n_people_qc"]),
                "protein_mapping_method": method,
                "n_mapping_candidates": int(len(candidates)),
                "alternative_assays": join_unique(candidates["Assay"].iloc[1:]),
            }
        )
    return pd.DataFrame(records)


def build_protein_z_matrix(
    olink_rows: pd.DataFrame,
    mapping: pd.DataFrame,
) -> pd.DataFrame:
    selected_assays = mapping["aou_assay"].dropna().astype(str).unique().tolist()
    sub = olink_rows[olink_rows["Assay"].isin(selected_assays)].copy()
    if sub.empty:
        return pd.DataFrame(columns=["person_id"])
    wide = (
        sub.groupby(["person_id", "Assay"])["protein_value"]
        .mean()
        .unstack()
    )
    means = wide.mean(axis=0)
    sds = wide.std(axis=0, ddof=0).replace(0, np.nan)
    z = (wide - means) / sds
    return z.reset_index()


def parse_manifest_assay_dates(path: Path) -> pd.DataFrame:
    manifest = pd.read_csv(path, sep="\t", low_memory=False)
    person_col = first_existing_column(
        manifest, ["research_id", "ResearchID", "sample_name", "person_id"], required=True
    )
    date = pd.Series(pd.NaT, index=manifest.index, dtype="datetime64[ns]")

    for candidate in [
        "SampleCollectionDate",
        "sample_collection_date",
        "CollectionDate",
        "collection_date",
        "AssayDate",
        "assay_date",
    ]:
        if candidate in manifest.columns:
            parsed = pd.to_datetime(manifest[candidate], errors="coerce")
            date = date.fillna(parsed)

    if "PlateID" in manifest.columns:
        parsed_plate = pd.to_datetime(
            manifest["PlateID"].astype(str).str.extract(r"_(\d{2}-\d{2}-\d{4})_")[0],
            format="%m-%d-%Y",
            errors="coerce",
        )
        date = date.fillna(parsed_plate)

    out = pd.DataFrame(
        {
            "person_id": pd.to_numeric(manifest[person_col], errors="coerce"),
            "assay_date": date,
        }
    ).dropna(subset=["person_id"])
    out["person_id"] = out["person_id"].astype("int64")
    return (
        out.dropna(subset=["assay_date"])
        .groupby("person_id", as_index=False)["assay_date"]
        .min()
    )


def calculate_age_at_date(
    people: pd.DataFrame,
    date_col: str,
    output_col: str,
) -> pd.DataFrame:
    """Calculate age at an event date using birth datetime or year fallback."""

    if date_col not in people.columns:
        raise ValueError(f"Missing event-date column: {date_col}")

    out = people.copy()
    date = pd.to_datetime(out[date_col], errors="coerce")

    if "birth_datetime" in out.columns:
        birth = pd.to_datetime(out["birth_datetime"], errors="coerce", utc=True)
        birth = birth.dt.tz_convert(None)
        age = (date - birth).dt.days / 365.25
    else:
        age = pd.Series(np.nan, index=out.index, dtype=float)

    if "year_of_birth" in out.columns:
        year = pd.to_numeric(out["year_of_birth"], errors="coerce")
        age = age.fillna(date.dt.year - year)

    out[output_col] = pd.to_numeric(age, errors="coerce")
    return out


# -----------------------------------------------------------------------------
# Chemistry extraction helpers
# -----------------------------------------------------------------------------


def bq_quote(value: Any) -> str:
    if _is_missing(value):
        value = ""
    escaped = str(value).replace("\\", "\\\\").replace("'", "\\'")
    return f"'{escaped}'"


def build_selected_strata_cte(selected: pd.DataFrame) -> str:
    required = [
        "source_feature",
        "component",
        "measurement_concept_id",
        "unit",
        "target_unit",
        "conversion_factor_to_target",
    ]
    missing = [c for c in required if c not in selected.columns]
    if missing:
        raise ValueError(f"Selected chemistry strata missing columns: {missing}")

    table = selected.copy()
    if "conversion_offset_to_target" not in table.columns:
        table["conversion_offset_to_target"] = 0.0
    table["conversion_offset_to_target"] = pd.to_numeric(
        table["conversion_offset_to_target"], errors="coerce"
    ).fillna(0.0)

    cte_columns = [*required, "conversion_offset_to_target"]
    rows = []
    for record in table[cte_columns].drop_duplicates().itertuples(index=False):
        rows.append(
            "STRUCT("
            f"{bq_quote(record.source_feature)} AS source_feature, "
            f"{bq_quote(record.component)} AS component, "
            f"CAST({int(record.measurement_concept_id)} AS INT64) AS measurement_concept_id, "
            f"{bq_quote(record.unit)} AS unit, "
            f"{bq_quote(record.target_unit)} AS target_unit, "
            f"CAST({float(record.conversion_factor_to_target):.17g} AS FLOAT64) AS conversion_factor_to_target, "
            f"CAST({float(record.conversion_offset_to_target):.17g} AS FLOAT64) AS conversion_offset_to_target"
            ")"
        )
    if not rows:
        raise ValueError("No complete included chemistry strata are available")
    return "selected_strata AS (\n  SELECT * FROM UNNEST([\n    " + ",\n    ".join(rows) + "\n  ])\n)"


def robust_zscore(series: pd.Series, q_low: float = 0.01, q_high: float = 0.99) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").astype(float)
    observed = values.dropna()
    if len(observed) < 3:
        return pd.Series(np.nan, index=values.index)
    lo, hi = observed.quantile([q_low, q_high])
    clipped = values.clip(lo, hi) if np.isfinite(lo) and np.isfinite(hi) and lo < hi else values
    sd = clipped.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=values.index)
    return (clipped - clipped.mean()) / sd


def normalize_selected_strata(
    legacy: pd.DataFrame,
    chemistry_spec: pd.DataFrame,
) -> pd.DataFrame:
    """Convert the old ``vae_feature`` strata table to generic source features."""

    aliases = chemistry_alias_lookup(chemistry_spec)
    table = legacy.copy()
    if "final_action" in table.columns:
        table = table[
            table["final_action"].astype(str).eq("include_for_normalized_combining")
        ].copy()
    source_col = first_existing_column(
        table, ["source_feature", "vae_feature", "feature"], required=True
    )
    table["source_feature"] = table[source_col].map(
        lambda x: canonicalize_chemistry_feature(x, aliases)
    )
    table["measurement_concept_id"] = pd.to_numeric(
        table["measurement_concept_id"], errors="coerce"
    )
    table["conversion_factor_to_target"] = pd.to_numeric(
        table["conversion_factor_to_target"], errors="coerce"
    )
    if "conversion_offset_to_target" not in table.columns:
        table["conversion_offset_to_target"] = 0.0
    table["conversion_offset_to_target"] = pd.to_numeric(
        table["conversion_offset_to_target"], errors="coerce"
    ).fillna(0.0)
    table["unit"] = table["unit"].fillna("Missing unit").astype(str)
    required_nonmissing = [
        "source_feature",
        "component",
        "measurement_concept_id",
        "target_unit",
        "conversion_factor_to_target",
    ]
    table = table.dropna(subset=required_nonmissing).copy()
    table["measurement_concept_id"] = table["measurement_concept_id"].astype("int64")
    return table[
        [
            "source_feature",
            "component",
            "measurement_concept_id",
            "unit",
            "target_unit",
            "conversion_factor_to_target",
            "conversion_offset_to_target",
        ]
    ].drop_duplicates()


def required_chemistry_components(
    target_features: Iterable[str],
    chemistry_spec: pd.DataFrame,
) -> pd.DataFrame:
    spec = chemistry_spec.set_index("canonical_feature", drop=False)
    rows = []
    for feature in sorted(set(target_features)):
        if feature not in spec.index:
            rows.append(
                {
                    "canonical_feature": feature,
                    "source_feature": feature,
                    "component": pd.NA,
                    "component_role": "unknown",
                }
            )
            continue
        row = spec.loc[feature]
        for number in [1, 2, 3]:
            sf = row.get(f"source_feature_{number}", pd.NA)
            comp = row.get(f"component_{number}", pd.NA)
            if not _is_missing(sf) and str(sf).strip() and not _is_missing(comp) and str(comp).strip():
                rows.append(
                    {
                        "canonical_feature": feature,
                        "source_feature": slugify(sf),
                        "component": str(comp).strip(),
                        "component_role": f"component_{number}",
                    }
                )
    return pd.DataFrame(rows).drop_duplicates()


def chemistry_component_availability(
    target_features: Iterable[str],
    chemistry_spec: pd.DataFrame,
    selected_strata: pd.DataFrame,
) -> pd.DataFrame:
    """Report whether every required source component has curated strata."""

    required = required_chemistry_components(target_features, chemistry_spec)
    if required.empty:
        return pd.DataFrame(
            columns=[
                "canonical_feature",
                "n_required_components",
                "n_available_components",
                "missing_components",
                "all_components_available",
            ]
        )

    available = selected_strata[["source_feature", "component"]].drop_duplicates()
    check = required.merge(
        available.assign(component_available=True),
        on=["source_feature", "component"],
        how="left",
    )
    check["component_available"] = check["component_available"].eq(True)

    rows = []
    for feature, group in check.groupby("canonical_feature", sort=False):
        missing = group.loc[
            ~group["component_available"], ["source_feature", "component"]
        ]
        rows.append(
            {
                "canonical_feature": feature,
                "n_required_components": int(len(group)),
                "n_available_components": int(group["component_available"].sum()),
                "missing_components": join_unique(
                    missing.apply(
                        lambda r: f"{r.source_feature}:{r.component}", axis=1
                    )
                ),
                "all_components_available": bool(group["component_available"].all()),
            }
        )
    return pd.DataFrame(rows)


def assemble_chemistry_feature_matrix(
    component_long: pd.DataFrame,
    target_features: Iterable[str],
    chemistry_spec: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create one z-scored person-level column per canonical chemistry feature."""

    raw = component_long.pivot_table(
        index="person_id",
        columns=["source_feature", "component"],
        values="value_raw_median",
        aggfunc="first",
    )
    z = component_long.pivot_table(
        index="person_id",
        columns=["source_feature", "component"],
        values="value_z_median",
        aggfunc="first",
    )
    people = raw.index.union(z.index).sort_values()
    raw = raw.reindex(people)
    z = z.reindex(people)

    def get(df: pd.DataFrame, source: Any, component: Any) -> pd.Series:
        key = (slugify(source), str(component))
        if key in df.columns:
            return pd.to_numeric(df[key], errors="coerce")
        return pd.Series(np.nan, index=people)

    spec = chemistry_spec.set_index("canonical_feature", drop=False)
    matrix = pd.DataFrame(index=people)
    metadata = []

    for feature in sorted(set(target_features)):
        if feature not in spec.index:
            metadata.append(
                {
                    "canonical_feature": feature,
                    "matrix_column": f"chem_{feature}_z",
                    "status": "not_in_feature_spec",
                    "n_people": 0,
                }
            )
            continue

        row = spec.loc[feature]
        derivation = str(row.get("derivation", "direct") or "direct").strip().lower()
        source1, comp1 = row.get("source_feature_1"), row.get("component_1")
        source2, comp2 = row.get("source_feature_2"), row.get("component_2")

        if derivation == "direct":
            # Direct biomarkers use the within-stratum z value, followed by one
            # robust standardization after combining strata/person medians.
            values = get(z, source1, comp1)
            feature_z = robust_zscore(values)
        elif derivation == "pulse_pressure":
            feature_z = robust_zscore(get(raw, source1, comp1) - get(raw, source2, comp2))
        elif derivation == "mean_arterial_pressure":
            feature_z = robust_zscore(
                (get(raw, source1, comp1) + 2 * get(raw, source2, comp2)) / 3
            )
        elif derivation == "ratio":
            denominator = get(raw, source2, comp2).replace(0, np.nan)
            feature_z = robust_zscore(get(raw, source1, comp1) / denominator)
        elif derivation == "homa_ir":
            # Formula assumes glucose in mg/dL and insulin in microIU/mL; the
            # curated strata target units must enforce those units.
            feature_z = robust_zscore(
                get(raw, source1, comp1) * get(raw, source2, comp2) / 405.0
            )
        else:
            feature_z = pd.Series(np.nan, index=people)

        column = f"chem_{feature}_z"
        matrix[column] = feature_z
        metadata.append(
            {
                "canonical_feature": feature,
                "display_label": row.get("display_label", feature),
                "matrix_column": column,
                "derivation": derivation,
                "source_components": "; ".join(
                    x
                    for x in [
                        f"{slugify(source1)}:{comp1}" if not _is_missing(source1) else "",
                        f"{slugify(source2)}:{comp2}" if not _is_missing(source2) else "",
                    ]
                    if x
                ),
                "status": "ok" if feature_z.notna().any() else "no_usable_values",
                "n_people": int(feature_z.notna().sum()),
                "mean": feature_z.mean(skipna=True),
                "sd": feature_z.std(skipna=True, ddof=0),
            }
        )

    matrix = matrix.reset_index().rename(columns={"index": "person_id"})
    return matrix, pd.DataFrame(metadata)


# -----------------------------------------------------------------------------
# Common regression engine
# -----------------------------------------------------------------------------


def zscore_series(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").astype(float)
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=values.index)
    return (values - values.mean()) / sd


def bh_fdr(p_values: Iterable[Any]) -> np.ndarray:
    p = pd.to_numeric(pd.Series(list(p_values)), errors="coerce").to_numpy(float)
    q = np.full(p.shape, np.nan, dtype=float)
    mask = np.isfinite(p)
    if not mask.any():
        return q
    valid = np.clip(p[mask], np.nextafter(0.0, 1.0), 1.0)
    order = np.argsort(valid)
    ranked = valid[order]
    m = len(ranked)
    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    q[mask] = restored
    return q


def add_bh_by_group(
    df: pd.DataFrame,
    group_cols: Sequence[str],
    p_col: str,
    q_col: str,
) -> pd.DataFrame:
    out = df.copy()
    out[q_col] = np.nan
    if p_col not in out.columns or out.empty:
        return out
    missing_groups = [column for column in group_cols if column not in out.columns]
    if missing_groups:
        raise ValueError(f"Missing BH grouping columns: {missing_groups}")
    for _, index in out.groupby(list(group_cols), dropna=False).groups.items():
        idx = list(index)
        out.loc[idx, q_col] = bh_fdr(out.loc[idx, p_col])
    return out


def load_model_plan(path: Path) -> pd.DataFrame:
    plan = pd.read_csv(path, sep="\t", dtype="string")
    required = {
        "analysis_name",
        "analysis_label",
        "is_primary",
        "control_set",
        "age_form",
        "age_subset",
        "include_sex",
        "n_pcs",
    }
    missing = sorted(required - set(plan.columns))
    if missing:
        raise ValueError(f"Model plan missing columns: {missing}")
    plan["is_primary"] = plan["is_primary"].str.lower().isin(["true", "1", "yes"])
    plan["include_sex"] = plan["include_sex"].str.lower().isin(["true", "1", "yes"])
    plan["n_pcs"] = pd.to_numeric(plan["n_pcs"], errors="raise").astype(int)
    return plan


def build_design_matrix(
    dat: pd.DataFrame,
    age_col: str,
    spec: Mapping[str, Any],
    pc_columns: Sequence[str],
) -> pd.DataFrame:
    """Create the common disease-case design matrix."""

    from patsy import dmatrix

    X = pd.DataFrame(index=dat.index)
    X["intercept"] = 1.0
    X["case"] = dat["case"].astype(float)

    age_form = str(spec["age_form"])
    if age_form == "spline":
        if dat[age_col].nunique(dropna=True) >= 6:
            basis = dmatrix(
                f"0 + bs({age_col}, df=4, degree=3, include_intercept=False, lower_bound=18, upper_bound=110)",
                data=dat,
                return_type="dataframe",
            )
            basis.index = dat.index
            basis.columns = [f"age_spline_{i+1}" for i in range(basis.shape[1])]
            X = pd.concat([X, basis], axis=1)
        else:
            X["age_z"] = zscore_series(dat[age_col])
    elif age_form == "linear":
        X["age_z"] = zscore_series(dat[age_col])
    elif age_form == "none":
        pass
    else:
        raise ValueError(f"Unknown age_form: {age_form}")

    if bool(spec["include_sex"]):
        sex_col = "sex_adjustment" if "sex_adjustment" in dat.columns else "sex"
        sex = dat[sex_col].fillna("Other/Unknown").astype(str)
        dummies = pd.get_dummies(sex, prefix="sex", drop_first=True, dtype=float)
        X = pd.concat([X, dummies], axis=1)

    n_pcs = int(spec["n_pcs"])
    for pc in list(pc_columns)[:n_pcs]:
        X[pc] = zscore_series(dat[pc])

    # Remove nuisance terms with no variation.  Case is checked separately.
    keep = ["intercept", "case"]
    keep += [
        c
        for c in X.columns
        if c not in {"intercept", "case"} and X[c].nunique(dropna=False) > 1
    ]
    return X[keep].astype(float)


def fit_case_control_ols(
    matrix: pd.DataFrame,
    outcome_col: str,
    domain: str,
    age_col: str,
    spec: Mapping[str, Any],
    pc_columns: Sequence[str],
    min_cases: int,
    min_controls: int,
    age_cutoff: float = 45.0,
) -> dict[str, Any]:
    """Fit one standardized omic outcome on disease case status using HC3 SEs."""

    import statsmodels.api as sm
    from scipy.stats import norm

    case_col = f"{domain}_case"
    control_set = str(spec["control_set"])
    control_col = (
        f"{domain}_domain_control_candidate"
        if control_set == "domain"
        else "clean_control_candidate"
        if control_set == "clean"
        else None
    )
    if control_col is None:
        raise ValueError(f"Unknown control_set: {control_set}")

    base = {
        "analysis_name": spec["analysis_name"],
        "analysis_label": spec["analysis_label"],
        "is_primary": bool(spec["is_primary"]),
        "control_set": control_set,
        "age_form": spec["age_form"],
        "age_subset": spec["age_subset"],
        "n_pcs": int(spec["n_pcs"]),
        "include_sex": bool(spec["include_sex"]),
        "domain": domain,
        "outcome_col": outcome_col,
        "status": "not_run",
        "n_people": 0,
        "n_cases": 0,
        "n_controls": 0,
        "beta_case": np.nan,
        "se_case": np.nan,
        "z_case": np.nan,
        "p_case": np.nan,
        "minus_log10_p": np.nan,
        "ci_low": np.nan,
        "ci_high": np.nan,
        "model_r2": np.nan,
        "condition_number": np.nan,
        "n_model_columns": np.nan,
    }

    needed = [outcome_col, age_col, case_col, control_col]
    if bool(spec["include_sex"]):
        needed.append("sex_adjustment" if "sex_adjustment" in matrix.columns else "sex")
    pcs = list(pc_columns)[: int(spec["n_pcs"])]
    needed += pcs
    missing = [c for c in needed if c not in matrix.columns]
    if missing:
        base["status"] = "missing_columns"
        base["error_detail"] = ";".join(missing)
        return base

    mask = matrix[case_col].eq(1) | matrix[control_col].eq(1)
    dat = matrix.loc[mask, list(dict.fromkeys(needed))].copy()
    dat["case"] = matrix.loc[mask, case_col].astype(int).values

    age_subset = str(spec["age_subset"])
    if age_subset == "lt45":
        dat = dat[dat[age_col] < age_cutoff]
    elif age_subset == "ge45":
        dat = dat[dat[age_col] >= age_cutoff]
    elif age_subset != "all":
        raise ValueError(f"Unknown age_subset: {age_subset}")

    dat = dat.dropna(subset=[outcome_col, age_col, "case", *pcs]).copy()
    n_cases = int(dat["case"].sum())
    n_controls = int(dat["case"].eq(0).sum())
    base.update(n_people=len(dat), n_cases=n_cases, n_controls=n_controls)

    if n_cases < min_cases or n_controls < min_controls:
        base["status"] = "skipped_low_n"
        return base
    if dat["case"].nunique() < 2:
        base["status"] = "no_case_control_variation"
        return base

    try:
        X = build_design_matrix(dat, age_col, spec, pc_columns)
        valid = X.notna().all(axis=1) & pd.to_numeric(dat[outcome_col], errors="coerce").notna()
        X = X.loc[valid]
        y = pd.to_numeric(dat.loc[valid, outcome_col], errors="coerce").astype(float)
        n_cases = int(dat.loc[valid, "case"].sum())
        n_controls = int(dat.loc[valid, "case"].eq(0).sum())
        base.update(n_people=len(y), n_cases=n_cases, n_controls=n_controls)
        if n_cases < min_cases or n_controls < min_controls:
            base["status"] = "skipped_low_n_after_covariates"
            return base
        if len(y) <= X.shape[1] + 1:
            base["status"] = "insufficient_residual_df"
            return base

        fit = sm.OLS(y, X).fit(cov_type="HC3")
        beta = float(fit.params["case"])
        se = float(fit.bse["case"])
        if not np.isfinite(se) or se <= 0:
            base["status"] = "invalid_case_se"
            return base
        z = beta / se
        log_p = math.log(2.0) + float(norm.logsf(abs(z)))
        tiny_log = math.log(np.nextafter(0.0, 1.0))
        p = float(math.exp(max(log_p, tiny_log)))
        minus_log10_p = float(-log_p / math.log(10.0))
        base.update(
            {
                "status": "ok",
                "beta_case": beta,
                "se_case": se,
                "z_case": z,
                "p_case": p,
                "minus_log10_p": minus_log10_p,
                "ci_low": beta - 1.96 * se,
                "ci_high": beta + 1.96 * se,
                "model_r2": float(fit.rsquared),
                "condition_number": float(np.linalg.cond(X.to_numpy())),
                "n_model_columns": int(X.shape[1]),
            }
        )
    except Exception as exc:  # status is retained in output for troubleshooting
        base["status"] = "model_error"
        base["error_detail"] = repr(exc)
    return base


def fit_case_spline_interaction_ols(
    matrix: pd.DataFrame,
    outcome_col: str,
    domain: str,
    age_col: str,
    spec: Mapping[str, Any],
    pc_columns: Sequence[str],
    min_cases: int,
    min_controls: int,
    reference_ages: Sequence[float] = (40.0, 55.0, 70.0),
) -> tuple[dict[str, Any], pd.DataFrame]:
    """Fit a disease-by-spline(age) interaction and age-specific contrasts.

    The joint p-value tests all case-by-spline terms simultaneously.  The
    reference-age contrasts estimate the adjusted standardized case-control
    difference at fixed ages and are intended for plotting, not for selecting
    the most favorable age.
    """

    import statsmodels.api as sm
    from patsy import build_design_matrices, dmatrix
    from scipy.stats import norm

    control_set = str(spec["control_set"])
    case_col = f"{domain}_case"
    control_col = (
        f"{domain}_domain_control_candidate"
        if control_set == "domain"
        else "clean_control_candidate"
        if control_set == "clean"
        else None
    )
    if control_col is None:
        raise ValueError(f"Unknown control_set: {control_set}")

    base: dict[str, Any] = {
        "analysis_name": "primary_spline_age_interaction",
        "analysis_label": (
            "Discovery-concordant spline test: disease x spline(age); "
            "sex + PC adjustment"
        ),
        "control_set": control_set,
        "n_pcs": int(spec["n_pcs"]),
        "include_sex": bool(spec["include_sex"]),
        "domain": domain,
        "outcome_col": outcome_col,
        "status": "not_run",
        "n_people": 0,
        "n_cases": 0,
        "n_controls": 0,
        "interaction_statistic": np.nan,
        "interaction_df": np.nan,
        "p_interaction": np.nan,
        "minus_log10_p_interaction": np.nan,
        "condition_number": np.nan,
        "n_model_columns": np.nan,
        "error_detail": pd.NA,
    }

    needed = [outcome_col, age_col, case_col, control_col]
    if bool(spec["include_sex"]):
        needed.append("sex_adjustment" if "sex_adjustment" in matrix.columns else "sex")
    pcs = list(pc_columns)[: int(spec["n_pcs"])]
    needed += pcs
    missing = [column for column in needed if column not in matrix.columns]
    if missing:
        base["status"] = "missing_columns"
        base["error_detail"] = ";".join(missing)
        return base, pd.DataFrame()

    mask = matrix[case_col].eq(1) | matrix[control_col].eq(1)
    dat = matrix.loc[mask, list(dict.fromkeys(needed))].copy()
    dat["case"] = matrix.loc[mask, case_col].astype(int).values
    dat = dat.dropna(subset=[outcome_col, age_col, "case", *pcs]).copy()

    n_cases = int(dat["case"].sum())
    n_controls = int(dat["case"].eq(0).sum())
    base.update(n_people=len(dat), n_cases=n_cases, n_controls=n_controls)
    if n_cases < min_cases or n_controls < min_controls:
        base["status"] = "skipped_low_n"
        return base, pd.DataFrame()
    if dat["case"].nunique() < 2:
        base["status"] = "no_case_control_variation"
        return base, pd.DataFrame()
    if dat[age_col].nunique(dropna=True) < 6:
        base["status"] = "insufficient_age_variation"
        return base, pd.DataFrame()

    try:
        basis = dmatrix(
            f"0 + bs({age_col}, df=4, degree=3, include_intercept=False, "
            "lower_bound=18, upper_bound=110)",
            data=dat,
            return_type="dataframe",
        )
        basis.index = dat.index
        design_info = basis.design_info
        basis.columns = [f"age_spline_{i + 1}" for i in range(basis.shape[1])]

        X = pd.DataFrame(index=dat.index)
        X["intercept"] = 1.0
        X["case"] = dat["case"].astype(float)
        X = pd.concat([X, basis], axis=1)

        interaction_columns = []
        for column in basis.columns:
            interaction = f"case_x_{column}"
            X[interaction] = X["case"] * X[column]
            interaction_columns.append(interaction)

        if bool(spec["include_sex"]):
            sex_col = "sex_adjustment" if "sex_adjustment" in dat.columns else "sex"
            sex = dat[sex_col].fillna("Other/Unknown").astype(str)
            dummies = pd.get_dummies(sex, prefix="sex", drop_first=True, dtype=float)
            X = pd.concat([X, dummies], axis=1)

        for pc in pcs:
            X[pc] = zscore_series(dat[pc])

        # Drop only nuisance columns with no variation; all interaction terms
        # must remain for the prespecified joint test.
        keep = ["intercept", "case", *basis.columns, *interaction_columns]
        keep += [
            column
            for column in X.columns
            if column not in keep and X[column].nunique(dropna=False) > 1
        ]
        X = X[keep].astype(float)
        valid = X.notna().all(axis=1) & pd.to_numeric(
            dat[outcome_col], errors="coerce"
        ).notna()
        X = X.loc[valid]
        y = pd.to_numeric(dat.loc[valid, outcome_col], errors="coerce").astype(float)
        n_cases = int(dat.loc[valid, "case"].sum())
        n_controls = int(dat.loc[valid, "case"].eq(0).sum())
        base.update(n_people=len(y), n_cases=n_cases, n_controls=n_controls)
        if n_cases < min_cases or n_controls < min_controls:
            base["status"] = "skipped_low_n_after_covariates"
            return base, pd.DataFrame()
        if len(y) <= X.shape[1] + 1:
            base["status"] = "insufficient_residual_df"
            return base, pd.DataFrame()

        fit = sm.OLS(y, X).fit(cov_type="HC3")
        interaction_positions = [X.columns.get_loc(column) for column in interaction_columns]
        restriction = np.zeros((len(interaction_positions), X.shape[1]))
        for row_index, column_index in enumerate(interaction_positions):
            restriction[row_index, column_index] = 1.0
        wald = fit.wald_test(restriction, use_f=False, scalar=True)
        statistic = float(np.asarray(wald.statistic).squeeze())
        p_joint = float(wald.pvalue)
        if not np.isfinite(p_joint):
            raise ValueError("Joint interaction test returned a non-finite p-value")
        p_joint = float(np.clip(p_joint, np.nextafter(0.0, 1.0), 1.0))

        base.update(
            {
                "status": "ok",
                "interaction_statistic": statistic,
                "interaction_df": len(interaction_positions),
                "p_interaction": p_joint,
                "minus_log10_p_interaction": -math.log10(p_joint),
                "condition_number": float(np.linalg.cond(X.to_numpy())),
                "n_model_columns": int(X.shape[1]),
                "observed_age_min": float(dat.loc[valid, age_col].min()),
                "observed_age_max": float(dat.loc[valid, age_col].max()),
            }
        )

        params = fit.params.to_numpy(float)
        covariance = fit.cov_params().to_numpy(float)
        contrast_rows = []
        for reference_age in reference_ages:
            new_basis = np.asarray(
                build_design_matrices(
                    [design_info],
                    {age_col: [float(reference_age)]},
                    return_type="dataframe",
                )[0]
            ).reshape(-1)
            contrast = np.zeros(X.shape[1], dtype=float)
            contrast[X.columns.get_loc("case")] = 1.0
            for value, column in zip(new_basis, interaction_columns):
                contrast[X.columns.get_loc(column)] = float(value)
            estimate = float(contrast @ params)
            variance = float(contrast @ covariance @ contrast)
            se = math.sqrt(max(variance, 0.0))
            if se > 0 and np.isfinite(se):
                z_value = estimate / se
                log_p = math.log(2.0) + float(norm.logsf(abs(z_value)))
                p_value = float(
                    math.exp(max(log_p, math.log(np.nextafter(0.0, 1.0))))
                )
            else:
                z_value = np.nan
                p_value = np.nan
            contrast_rows.append(
                {
                    "reference_age": float(reference_age),
                    "beta_case_at_age": estimate,
                    "se_case_at_age": se,
                    "z_case_at_age": z_value,
                    "p_case_at_age": p_value,
                    "ci_low_at_age": estimate - 1.96 * se,
                    "ci_high_at_age": estimate + 1.96 * se,
                    "age_within_observed_range": bool(
                        base["observed_age_min"] <= reference_age <= base["observed_age_max"]
                    ),
                }
            )
        return base, pd.DataFrame(contrast_rows)
    except Exception as exc:
        base["status"] = "model_error"
        base["error_detail"] = repr(exc)
        return base, pd.DataFrame()


def run_unique_spline_interaction_fits(
    matrix: pd.DataFrame,
    analysis_targets: pd.DataFrame,
    feature_column_map: Mapping[str, str],
    omic_type: str,
    age_col: str,
    primary_spec: Mapping[str, Any],
    pc_columns: Sequence[str],
    min_cases: int,
    min_controls: int,
    reference_ages: Sequence[float] = (40.0, 55.0, 70.0),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run discovery-concordant disease-by-age tests for spline candidates."""

    required_target_columns = {
        "clinically_testable",
        "omic_type",
        "analysis_domain",
        "feature_key",
        "discovery_model",
    }
    missing_target_columns = sorted(
        required_target_columns - set(analysis_targets.columns)
    )
    if missing_target_columns:
        raise ValueError(
            "analysis_targets is missing required spline columns: "
            f"{missing_target_columns}"
        )

    spline_targets = analysis_targets[
        analysis_targets["clinically_testable"].fillna(False)
        & analysis_targets["omic_type"].eq(omic_type)
        & analysis_targets["discovery_model"].eq("Spline")
    ].copy()
    pairs = spline_targets[["analysis_domain", "feature_key"]].drop_duplicates()
    fit_rows: list[dict[str, Any]] = []
    contrast_rows: list[pd.DataFrame] = []
    contrast_test_ids: set[str] = set()
    cache: dict[tuple[str, str], tuple[dict[str, Any], pd.DataFrame]] = {}

    for pair in pairs.itertuples(index=False):
        outcome_col = feature_column_map.get(pair.feature_key)
        link_id = stable_id(
            "spline_fit_link",
            [omic_type, pair.analysis_domain, pair.feature_key],
        )
        if not outcome_col:
            unique_test_id = stable_id(
                "spline_test_unavailable",
                [omic_type, pair.analysis_domain, pair.feature_key],
            )
            result = {
                "spline_fit_link_id": link_id,
                "unique_spline_test_id": unique_test_id,
                "omic_type": omic_type,
                "analysis_domain": pair.analysis_domain,
                "feature_key": pair.feature_key,
                "outcome_col": pd.NA,
                "analysis_name": "primary_spline_age_interaction",
                "analysis_label": "Discovery-concordant spline interaction",
                "status_interaction": "feature_not_available",
                "p_interaction": np.nan,
            }
            fit_rows.append(result)
            continue

        unique_test_id = stable_id(
            "spline_test",
            [omic_type, pair.analysis_domain, outcome_col],
        )
        cache_key = (str(pair.analysis_domain), str(outcome_col))
        if cache_key not in cache:
            cache[cache_key] = fit_case_spline_interaction_ols(
                matrix=matrix,
                outcome_col=outcome_col,
                domain=pair.analysis_domain,
                age_col=age_col,
                spec=primary_spec,
                pc_columns=pc_columns,
                min_cases=min_cases,
                min_controls=min_controls,
                reference_ages=reference_ages,
            )
        fit_result, contrasts = cache[cache_key]
        result = dict(fit_result)
        result["status_interaction"] = result.pop("status")
        result.update(
            {
                "spline_fit_link_id": link_id,
                "unique_spline_test_id": unique_test_id,
                "omic_type": omic_type,
                "analysis_domain": pair.analysis_domain,
                "feature_key": pair.feature_key,
            }
        )
        fit_rows.append(result)

        if not contrasts.empty and unique_test_id not in contrast_test_ids:
            contrasts = contrasts.copy()
            contrasts["unique_spline_test_id"] = unique_test_id
            contrasts["omic_type"] = omic_type
            contrasts["analysis_domain"] = pair.analysis_domain
            contrasts["outcome_col"] = outcome_col
            contrast_rows.append(contrasts)
            contrast_test_ids.add(unique_test_id)

    fit_schema = [
        "spline_fit_link_id",
        "unique_spline_test_id",
        "omic_type",
        "analysis_domain",
        "feature_key",
        "outcome_col",
        "analysis_name",
        "analysis_label",
        "status_interaction",
        "n_people",
        "n_cases",
        "n_controls",
        "interaction_statistic",
        "interaction_df",
        "p_interaction",
        "minus_log10_p_interaction",
        "condition_number",
        "n_model_columns",
        "observed_age_min",
        "observed_age_max",
        "error_detail",
    ]
    contrast_schema = [
        "unique_spline_test_id",
        "omic_type",
        "analysis_domain",
        "outcome_col",
        "reference_age",
        "beta_case_at_age",
        "se_case_at_age",
        "z_case_at_age",
        "p_case_at_age",
        "ci_low_at_age",
        "ci_high_at_age",
        "age_within_observed_range",
    ]
    fits = pd.DataFrame(fit_rows)
    contrasts = (
        pd.concat(contrast_rows, ignore_index=True)
        if contrast_rows
        else pd.DataFrame(columns=contrast_schema)
    )
    for column in fit_schema:
        if column not in fits.columns:
            fits[column] = np.nan
    fits = fits[[*fit_schema, *[c for c in fits.columns if c not in fit_schema]]]
    for column in contrast_schema:
        if column not in contrasts.columns:
            contrasts[column] = np.nan
    contrasts = contrasts[
        [*contrast_schema, *[c for c in contrasts.columns if c not in contrast_schema]]
    ]
    if fits.empty:
        for column in [
            "q_interaction_within_omic_domain",
            "q_interaction_within_omic",
            "q_interaction_global",
        ]:
            fits[column] = np.nan
            contrasts[column] = np.nan
        return fits, contrasts

    if "p_interaction" not in fits.columns:
        fits["p_interaction"] = np.nan
    q_columns = [
        "q_interaction_within_omic_domain",
        "q_interaction_within_omic",
        "q_interaction_global",
    ]
    unique_tests = fits.drop_duplicates("unique_spline_test_id").copy()
    unique_tests = add_bh_by_group(
        unique_tests,
        ["omic_type", "analysis_domain"],
        "p_interaction",
        q_columns[0],
    )
    unique_tests = add_bh_by_group(
        unique_tests,
        ["omic_type"],
        "p_interaction",
        q_columns[1],
    )
    unique_tests[q_columns[2]] = np.nan
    ok = unique_tests["status_interaction"].eq("ok") & pd.to_numeric(
        unique_tests["p_interaction"], errors="coerce"
    ).notna()
    unique_tests.loc[ok, q_columns[2]] = bh_fdr(
        unique_tests.loc[ok, "p_interaction"]
    )
    q_map = unique_tests[["unique_spline_test_id", *q_columns]]
    fits = fits.merge(
        q_map,
        on="unique_spline_test_id",
        how="left",
        validate="many_to_one",
    )
    if not contrasts.empty:
        contrasts = contrasts.merge(
            q_map,
            on="unique_spline_test_id",
            how="left",
            validate="many_to_one",
        )
    return fits, contrasts



def attach_discovery_concordant_results(
    primary_main_effect_results: pd.DataFrame,
    spline_fit_links: pd.DataFrame,
    nominal_threshold: float = 0.05,
    fdr_threshold: float = 0.10,
) -> pd.DataFrame:
    """Attach the discovery-concordant replication test to each primary target.

    VAE and OLS targets use the common adjusted disease main-effect model.
    Spline discovery targets instead use the prespecified joint disease-by-
    spline(age) interaction test.  The common main-effect result remains in the
    table for cross-method comparability; these additional columns identify the
    test that most closely matches each discovery hypothesis.
    """

    required_primary = {
        "analysis_target_id",
        "discovery_model",
        "omic_type",
        "analysis_domain",
        "feature_key",
        "status",
        "p_case",
        "q_within_analysis_omic_domain",
        "replication_category",
        "unique_aou_test_id",
    }
    missing_primary = sorted(required_primary - set(primary_main_effect_results.columns))
    if missing_primary:
        raise ValueError(
            "Primary results are missing required columns: " + ", ".join(missing_primary)
        )

    out = primary_main_effect_results.copy()
    spline = spline_fit_links.copy()
    spline_keys = ["omic_type", "analysis_domain", "feature_key"]

    if not spline.empty:
        required_spline = {
            *spline_keys,
            "unique_spline_test_id",
            "status_interaction",
            "p_interaction",
            "q_interaction_within_omic_domain",
        }
        missing_spline = sorted(required_spline - set(spline.columns))
        if missing_spline:
            raise ValueError(
                "Spline fit links are missing required columns: "
                + ", ".join(missing_spline)
            )
        if spline.duplicated(spline_keys).any():
            duplicates = spline.loc[
                spline.duplicated(spline_keys, keep=False), spline_keys
            ].drop_duplicates()
            raise ValueError(
                "Spline fit links are not unique by omic/domain/feature. "
                f"First duplicate keys: {duplicates.head().to_dict(orient='records')}"
            )

        retain = spline_keys + [
            column
            for column in [
                "spline_fit_link_id",
                "unique_spline_test_id",
                "outcome_col",
                "status_interaction",
                "interaction_statistic",
                "interaction_df",
                "p_interaction",
                "minus_log10_p_interaction",
                "q_interaction_within_omic_domain",
                "q_interaction_within_omic",
                "q_interaction_global",
                "n_people",
                "n_cases",
                "n_controls",
                "condition_number",
                "n_model_columns",
                "observed_age_min",
                "observed_age_max",
                "error_detail",
            ]
            if column in spline.columns
        ]
        spline = spline[retain].copy()
        spline["discovery_model"] = "Spline"
        rename = {
            "outcome_col": "interaction_outcome_col",
            "n_people": "interaction_n_people",
            "n_cases": "interaction_n_cases",
            "n_controls": "interaction_n_controls",
            "condition_number": "interaction_condition_number",
            "n_model_columns": "interaction_n_model_columns",
            "error_detail": "interaction_error_detail",
        }
        spline = spline.rename(columns=rename)
        out = out.merge(
            spline,
            on=[*spline_keys, "discovery_model"],
            how="left",
            validate="many_to_one",
        )

    # Ensure columns exist even when there were no spline targets/fits.
    interaction_defaults: dict[str, Any] = {
        "spline_fit_link_id": pd.NA,
        "unique_spline_test_id": pd.NA,
        "interaction_outcome_col": pd.NA,
        "status_interaction": pd.NA,
        "interaction_statistic": np.nan,
        "interaction_df": np.nan,
        "p_interaction": np.nan,
        "minus_log10_p_interaction": np.nan,
        "q_interaction_within_omic_domain": np.nan,
        "q_interaction_within_omic": np.nan,
        "q_interaction_global": np.nan,
        "interaction_n_people": np.nan,
        "interaction_n_cases": np.nan,
        "interaction_n_controls": np.nan,
        "interaction_condition_number": np.nan,
        "interaction_n_model_columns": np.nan,
        "observed_age_min": np.nan,
        "observed_age_max": np.nan,
        "interaction_error_detail": pd.NA,
    }
    for column, default in interaction_defaults.items():
        if column not in out.columns:
            out[column] = default

    is_spline = out["discovery_model"].eq("Spline")
    main_p = pd.to_numeric(out["p_case"], errors="coerce")
    main_q = pd.to_numeric(out["q_within_analysis_omic_domain"], errors="coerce")
    interaction_p = pd.to_numeric(out["p_interaction"], errors="coerce")
    interaction_q = pd.to_numeric(
        out["q_interaction_within_omic_domain"], errors="coerce"
    )

    out["replication_test_family"] = "disease main effect"
    out.loc[is_spline, "replication_test_family"] = (
        "disease x spline(age) interaction"
    )
    out["discovery_concordant_status"] = out["status"]
    out.loc[is_spline, "discovery_concordant_status"] = out.loc[
        is_spline, "status_interaction"
    ]
    out["discovery_concordant_p"] = main_p
    out.loc[is_spline, "discovery_concordant_p"] = interaction_p.loc[is_spline]
    out["discovery_concordant_q"] = main_q
    out.loc[is_spline, "discovery_concordant_q"] = interaction_q.loc[is_spline]
    out["discovery_concordant_q_scope"] = (
        "within primary main-effect tests x omic x disease domain"
    )
    out.loc[is_spline, "discovery_concordant_q_scope"] = (
        "within spline-interaction tests x omic x disease domain"
    )
    out["discovery_concordant_unique_test_id"] = out["unique_aou_test_id"]
    out.loc[is_spline, "discovery_concordant_unique_test_id"] = out.loc[
        is_spline, "unique_spline_test_id"
    ]
    out["discovery_concordant_test_statistic"] = pd.to_numeric(
        out.get("z_case"), errors="coerce"
    )
    out.loc[is_spline, "discovery_concordant_test_statistic"] = pd.to_numeric(
        out.loc[is_spline, "interaction_statistic"], errors="coerce"
    )
    out["discovery_concordant_test_df"] = np.where(
        out["status"].eq("ok"), 1.0, np.nan
    )
    out.loc[is_spline, "discovery_concordant_test_df"] = pd.to_numeric(
        out.loc[is_spline, "interaction_df"], errors="coerce"
    )
    out["discovery_concordant_n_cases"] = pd.to_numeric(
        out.get("n_cases"), errors="coerce"
    )
    out["discovery_concordant_n_controls"] = pd.to_numeric(
        out.get("n_controls"), errors="coerce"
    )
    out.loc[is_spline, "discovery_concordant_n_cases"] = pd.to_numeric(
        out.loc[is_spline, "interaction_n_cases"], errors="coerce"
    )
    out.loc[is_spline, "discovery_concordant_n_controls"] = pd.to_numeric(
        out.loc[is_spline, "interaction_n_controls"], errors="coerce"
    )
    out["discovery_concordant_direction_relation"] = out.get(
        "direction_relation", pd.Series("not_evaluable", index=out.index)
    )
    out.loc[is_spline, "discovery_concordant_direction_relation"] = "not_evaluable"

    category = out["replication_category"].fillna("not_tested").astype(str)
    category = category.where(~is_spline, "not_tested")
    spline_ok = is_spline & out["discovery_concordant_status"].eq("ok")
    category.loc[spline_ok] = "not_significant"
    category.loc[spline_ok & interaction_p.lt(nominal_threshold)] = (
        "nominal_nondirectional_age_interaction"
    )
    category.loc[spline_ok & interaction_q.lt(fdr_threshold)] = (
        "fdr10_nondirectional_age_interaction"
    )
    out["discovery_concordant_replication_category"] = category

    return out


def discovery_concordant_method_summary(primary: pd.DataFrame) -> pd.DataFrame:
    """Summarize replication using each discovery model's concordant AoU test."""

    required = {
        "discovery_framework",
        "discovery_model",
        "omic_type",
        "analysis_target_id",
        "discovery_concordant_status",
        "discovery_concordant_p",
        "discovery_concordant_q",
    }
    missing = sorted(required - set(primary.columns))
    if missing:
        raise ValueError(
            "Discovery-concordant summary is missing columns: " + ", ".join(missing)
        )

    rows: list[dict[str, Any]] = []
    group_cols = ["discovery_framework", "discovery_model", "omic_type"]
    for keys, group in primary.groupby(group_cols, dropna=False):
        status = group["discovery_concordant_status"].astype("string")
        ok = status.eq("ok").fillna(False)
        p = pd.to_numeric(group["discovery_concordant_p"], errors="coerce")
        q = pd.to_numeric(group["discovery_concordant_q"], errors="coerce")
        directional = (
            group.get("direction_evaluable", pd.Series(False, index=group.index))
            .fillna(False)
            .astype(bool)
            & ~group["discovery_model"].eq("Spline")
        )
        same = group.get(
            "discovery_concordant_direction_relation",
            pd.Series("not_evaluable", index=group.index),
        ).eq("same_direction")
        rows.append(
            {
                **dict(zip(group_cols, keys)),
                "replication_test_family": join_unique(
                    group.get(
                        "replication_test_family",
                        pd.Series(pd.NA, index=group.index),
                    )
                ),
                "n_analysis_targets": int(group["analysis_target_id"].nunique()),
                "n_modelled": int(ok.sum()),
                "n_nominal": int((ok & p.lt(0.05)).sum()),
                "n_fdr10": int((ok & q.lt(0.10)).sum()),
                "n_direction_evaluable": int((ok & directional).sum()),
                "n_nominal_directional": int(
                    (ok & directional & same & p.lt(0.05)).sum()
                ),
                "n_fdr10_directional": int(
                    (ok & directional & same & q.lt(0.10)).sum()
                ),
                "median_n_cases": pd.to_numeric(
                    group.get("discovery_concordant_n_cases"), errors="coerce"
                ).median(),
                "median_n_controls": pd.to_numeric(
                    group.get("discovery_concordant_n_controls"), errors="coerce"
                ).median(),
            }
        )
    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary["nominal_fraction_of_modelled"] = (
            summary["n_nominal"] / summary["n_modelled"].replace(0, np.nan)
        )
        summary["fdr10_fraction_of_modelled"] = (
            summary["n_fdr10"] / summary["n_modelled"].replace(0, np.nan)
        )
        summary["directional_nominal_fraction"] = (
            summary["n_nominal_directional"]
            / summary["n_direction_evaluable"].replace(0, np.nan)
        )
    return summary


def build_discovery_concordant_forest_rows(
    primary: pd.DataFrame,
    spline_contrasts: pd.DataFrame,
    p_threshold: float = 0.05,
    require_fdr10: bool = False,
    max_targets: Optional[int] = None,
) -> pd.DataFrame:
    """Construct plot rows without selecting a favorable spline reference age.

    Main-effect targets contribute one adjusted case coefficient. Spline targets
    are selected only by their prespecified joint interaction p/q value and then
    contribute *all* fixed, in-range reference-age contrasts.  Age-specific
    contrast p-values are retained for context but are never used for selection.
    """

    required = {
        "analysis_target_id",
        "discovery_model",
        "discovery_concordant_status",
        "discovery_concordant_p",
        "discovery_concordant_q",
        "beta_case",
        "ci_low",
        "ci_high",
    }
    missing = sorted(required - set(primary.columns))
    if missing:
        raise ValueError("Primary plot input is missing columns: " + ", ".join(missing))

    selected = primary[primary["discovery_concordant_status"].eq("ok")].copy()
    selected["discovery_concordant_p"] = pd.to_numeric(
        selected["discovery_concordant_p"], errors="coerce"
    )
    selected["discovery_concordant_q"] = pd.to_numeric(
        selected["discovery_concordant_q"], errors="coerce"
    )
    selected = selected[
        selected["discovery_concordant_p"].lt(float(p_threshold))
    ].copy()
    if require_fdr10:
        selected = selected[selected["discovery_concordant_q"].lt(0.10)].copy()
    if max_targets is not None:
        selected = selected.nsmallest(int(max_targets), "discovery_concordant_p")

    main = selected[~selected["discovery_model"].eq("Spline")].copy()
    main["plot_estimate"] = pd.to_numeric(main["beta_case"], errors="coerce")
    main["plot_ci_low"] = pd.to_numeric(main["ci_low"], errors="coerce")
    main["plot_ci_high"] = pd.to_numeric(main["ci_high"], errors="coerce")
    main["plot_reference_age"] = np.nan
    main["plot_reference_age_label"] = pd.NA
    main["plot_estimate_family"] = "adjusted disease main effect"
    main["age_specific_contrast_p"] = np.nan
    main["age_within_observed_range"] = pd.NA

    spline = selected[selected["discovery_model"].eq("Spline")].copy()
    if spline.empty:
        spline_rows = spline
    else:
        if spline_contrasts.empty:
            raise ValueError(
                "Spline targets passed the plot threshold but no spline contrast "
                "table was supplied"
            )
        contrast_required = {
            "unique_spline_test_id",
            "reference_age",
            "beta_case_at_age",
            "ci_low_at_age",
            "ci_high_at_age",
            "p_case_at_age",
            "age_within_observed_range",
        }
        missing_contrast = sorted(contrast_required - set(spline_contrasts.columns))
        if missing_contrast:
            raise ValueError(
                "Spline contrast table is missing columns: "
                + ", ".join(missing_contrast)
            )
        contrast_cols = list(contrast_required)
        contrasts = spline_contrasts[contrast_cols].drop_duplicates().copy()
        spline_rows = spline.merge(
            contrasts,
            on="unique_spline_test_id",
            how="left",
            validate="many_to_many",
        )
        spline_rows = spline_rows[
            spline_rows["age_within_observed_range"].fillna(False)
        ].copy()
        spline_rows["plot_estimate"] = pd.to_numeric(
            spline_rows["beta_case_at_age"], errors="coerce"
        )
        spline_rows["plot_ci_low"] = pd.to_numeric(
            spline_rows["ci_low_at_age"], errors="coerce"
        )
        spline_rows["plot_ci_high"] = pd.to_numeric(
            spline_rows["ci_high_at_age"], errors="coerce"
        )
        spline_rows["plot_reference_age"] = pd.to_numeric(
            spline_rows["reference_age"], errors="coerce"
        )
        spline_rows["plot_reference_age_label"] = spline_rows[
            "plot_reference_age"
        ].map(lambda value: f"age {value:g}" if pd.notna(value) else pd.NA)
        spline_rows["plot_estimate_family"] = (
            "adjusted case-control contrast from joint age interaction"
        )
        spline_rows["age_specific_contrast_p"] = pd.to_numeric(
            spline_rows["p_case_at_age"], errors="coerce"
        )

    plot_rows = pd.concat([main, spline_rows], ignore_index=True, sort=False)
    if plot_rows.empty:
        return plot_rows
    plot_rows["plot_p"] = pd.to_numeric(
        plot_rows["discovery_concordant_p"], errors="coerce"
    )
    plot_rows["plot_q"] = pd.to_numeric(
        plot_rows["discovery_concordant_q"], errors="coerce"
    )
    plot_rows["plot_selection_test"] = plot_rows["replication_test_family"]
    plot_rows["plot_direction_relation"] = plot_rows.get(
        "discovery_concordant_direction_relation",
        pd.Series("not_evaluable", index=plot_rows.index),
    )
    return plot_rows[
        plot_rows[["plot_estimate", "plot_ci_low", "plot_ci_high"]]
        .notna()
        .all(axis=1)
    ].reset_index(drop=True)

def run_unique_replication_fits(
    matrix: pd.DataFrame,
    analysis_targets: pd.DataFrame,
    feature_column_map: Mapping[str, str],
    omic_type: str,
    age_col: str,
    model_plan: pd.DataFrame,
    pc_columns: Sequence[str],
    min_cases: int,
    min_controls: int,
    age_cutoff: float = 45.0,
) -> pd.DataFrame:
    """Fit each unique AoU outcome/domain/model combination only once.

    Different discovery target identifiers can map to the same AoU outcome
    column (for example, a legacy gene-symbol key and a newer UniProt key).
    The regression is cached by domain + outcome column + model, then linked
    back to every discovery feature key.  FDR is calculated on unique AoU tests,
    not on the expanded target-link rows.
    """

    required_target_columns = {
        "clinically_testable",
        "omic_type",
        "analysis_domain",
        "feature_key",
    }
    missing_target_columns = sorted(required_target_columns - set(analysis_targets.columns))
    if missing_target_columns:
        raise ValueError(
            f"analysis_targets is missing required columns: {missing_target_columns}"
        )

    testable = analysis_targets[
        analysis_targets["clinically_testable"].fillna(False)
        & analysis_targets["omic_type"].eq(omic_type)
    ].copy()
    pairs = testable[["analysis_domain", "feature_key"]].drop_duplicates()

    rows: list[dict[str, Any]] = []
    fit_cache: dict[tuple[str, str, str], dict[str, Any]] = {}
    total = len(pairs) * len(model_plan)
    counter = 0

    for pair in pairs.itertuples(index=False):
        outcome_col = feature_column_map.get(pair.feature_key)
        for spec_row in model_plan.to_dict(orient="records"):
            counter += 1
            if counter == 1 or counter % 25 == 0 or counter == total:
                print(f"{omic_type}: target link {counter}/{total}")

            link_fit_id = stable_id(
                "fit_link",
                [
                    omic_type,
                    pair.analysis_domain,
                    pair.feature_key,
                    spec_row["analysis_name"],
                ],
            )

            if not outcome_col:
                unique_test_id = stable_id(
                    "aou_test_unavailable",
                    [
                        omic_type,
                        pair.analysis_domain,
                        pair.feature_key,
                        spec_row["analysis_name"],
                    ],
                )
                result = {
                    "fit_id": link_fit_id,
                    "unique_aou_test_id": unique_test_id,
                    "omic_type": omic_type,
                    "analysis_domain": pair.analysis_domain,
                    "feature_key": pair.feature_key,
                    "outcome_col": pd.NA,
                    "analysis_name": spec_row["analysis_name"],
                    "analysis_label": spec_row["analysis_label"],
                    "is_primary": bool(spec_row["is_primary"]),
                    "control_set": spec_row["control_set"],
                    "age_form": spec_row["age_form"],
                    "age_subset": spec_row["age_subset"],
                    "n_pcs": int(spec_row["n_pcs"]),
                    "include_sex": bool(spec_row["include_sex"]),
                    "domain": pair.analysis_domain,
                    "status": "feature_not_available",
                }
            else:
                unique_test_id = stable_id(
                    "aou_test",
                    [
                        omic_type,
                        pair.analysis_domain,
                        outcome_col,
                        spec_row["analysis_name"],
                    ],
                )
                cache_key = (
                    str(pair.analysis_domain),
                    str(outcome_col),
                    str(spec_row["analysis_name"]),
                )
                if cache_key not in fit_cache:
                    fit_cache[cache_key] = fit_case_control_ols(
                        matrix=matrix,
                        outcome_col=outcome_col,
                        domain=pair.analysis_domain,
                        age_col=age_col,
                        spec=spec_row,
                        pc_columns=pc_columns,
                        min_cases=min_cases,
                        min_controls=min_controls,
                        age_cutoff=age_cutoff,
                    )
                result = dict(fit_cache[cache_key])
                result.update(
                    {
                        "fit_id": link_fit_id,
                        "unique_aou_test_id": unique_test_id,
                        "omic_type": omic_type,
                        "analysis_domain": pair.analysis_domain,
                        "feature_key": pair.feature_key,
                    }
                )
            rows.append(result)

    fits = pd.DataFrame(rows)
    if fits.empty:
        return fits

    numeric_result_columns = [
        "n_people",
        "n_cases",
        "n_controls",
        "beta_case",
        "se_case",
        "z_case",
        "p_case",
        "minus_log10_p",
        "ci_low",
        "ci_high",
        "model_r2",
        "condition_number",
        "n_model_columns",
    ]
    for column in numeric_result_columns:
        if column not in fits.columns:
            fits[column] = np.nan

    q_columns = [
        "q_within_analysis_omic_domain",
        "q_within_analysis_omic",
        "q_global_all_unique_fits",
    ]

    # One row per actual AoU test is used as the FDR denominator.  The resulting
    # q-values are then copied back to all discovery feature keys linked to that
    # same outcome.
    unique_tests = fits.drop_duplicates("unique_aou_test_id").copy()
    unique_tests = add_bh_by_group(
        unique_tests,
        ["analysis_name", "omic_type", "analysis_domain"],
        "p_case",
        q_columns[0],
    )
    unique_tests = add_bh_by_group(
        unique_tests,
        ["analysis_name", "omic_type"],
        "p_case",
        q_columns[1],
    )
    unique_tests[q_columns[2]] = np.nan
    ok_unique = (
        unique_tests["status"].eq("ok")
        & pd.to_numeric(unique_tests["p_case"], errors="coerce").notna()
    )
    unique_tests.loc[ok_unique, q_columns[2]] = bh_fdr(
        unique_tests.loc[ok_unique, "p_case"]
    )

    q_map = unique_tests[["unique_aou_test_id", *q_columns]]
    fits = fits.drop(columns=[c for c in q_columns if c in fits.columns]).merge(
        q_map,
        on="unique_aou_test_id",
        how="left",
        validate="many_to_one",
    )
    return fits


def attach_fits_to_targets(
    analysis_targets: pd.DataFrame,
    fits: pd.DataFrame,
) -> pd.DataFrame:
    merged = analysis_targets.merge(
        fits,
        on=["omic_type", "analysis_domain", "feature_key"],
        how="left",
        validate="many_to_many",
        suffixes=("_discovery", "_aou"),
    )
    merged["aou_direction"] = np.where(
        pd.to_numeric(merged.get("beta_case"), errors="coerce") > 0,
        "positive",
        np.where(
            pd.to_numeric(merged.get("beta_case"), errors="coerce") < 0,
            "negative",
            "not_defined",
        ),
    )
    merged["direction_relation"] = "not_evaluable"
    evaluable = (
        merged["direction_evaluable"].fillna(False)
        & merged["status"].eq("ok")
        & merged["aou_direction"].isin(["positive", "negative"])
    )
    merged.loc[evaluable, "direction_relation"] = np.where(
        merged.loc[evaluable, "discovery_direction"].eq(
            merged.loc[evaluable, "aou_direction"]
        ),
        "same_direction",
        "opposite_direction",
    )

    p = pd.to_numeric(merged.get("p_case"), errors="coerce")
    q = pd.to_numeric(
        merged.get("q_within_analysis_omic_domain"), errors="coerce"
    )
    merged["replication_category"] = "not_tested"
    ok = merged["status"].eq("ok")
    merged.loc[ok, "replication_category"] = "not_significant"
    merged.loc[ok & p.lt(0.05) & merged["direction_relation"].eq("same_direction"), "replication_category"] = "nominal_directional"
    merged.loc[ok & q.lt(0.10) & merged["direction_relation"].eq("same_direction"), "replication_category"] = "fdr10_directional"
    merged.loc[ok & p.lt(0.05) & merged["direction_relation"].eq("opposite_direction"), "replication_category"] = "nominal_opposite_direction"
    merged.loc[ok & p.lt(0.05) & merged["direction_relation"].eq("not_evaluable"), "replication_category"] = "nominal_nondirectional"
    merged.loc[ok & q.lt(0.10) & merged["direction_relation"].eq("not_evaluable"), "replication_category"] = "fdr10_nondirectional"
    return merged


def primary_result_table(
    combined_results: pd.DataFrame,
    primary_analysis_name: str,
) -> pd.DataFrame:
    primary = combined_results[
        combined_results["analysis_name"].eq(primary_analysis_name)
    ].copy()
    return primary.sort_values(
        ["omic_type", "analysis_domain", "p_case", "discovery_framework", "discovery_model"],
        na_position="last",
    ).reset_index(drop=True)


def method_summary(primary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_cols = ["discovery_framework", "discovery_model", "omic_type"]
    for keys, group in primary.groupby(group_cols, dropna=False):
        ok = group["status"].eq("ok")
        p = pd.to_numeric(group["p_case"], errors="coerce")
        q = pd.to_numeric(group["q_within_analysis_omic_domain"], errors="coerce")
        directional = group["direction_evaluable"].fillna(False)
        same = group["direction_relation"].eq("same_direction")
        rows.append(
            {
                **dict(zip(group_cols, keys)),
                "n_analysis_targets": len(group),
                "n_modelled": int(ok.sum()),
                "n_nominal_any_direction": int((ok & p.lt(0.05)).sum()),
                "n_fdr10_any_direction": int((ok & q.lt(0.10)).sum()),
                "n_direction_evaluable": int((ok & directional).sum()),
                "n_nominal_directional": int((ok & directional & same & p.lt(0.05)).sum()),
                "n_fdr10_directional": int((ok & directional & same & q.lt(0.10)).sum()),
                "median_n_cases": pd.to_numeric(group.get("n_cases"), errors="coerce").median(),
                "median_n_controls": pd.to_numeric(group.get("n_controls"), errors="coerce").median(),
            }
        )
    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary["nominal_fraction_of_modelled"] = summary["n_nominal_any_direction"] / summary["n_modelled"].replace(0, np.nan)
        summary["directional_nominal_fraction"] = summary["n_nominal_directional"] / summary["n_direction_evaluable"].replace(0, np.nan)
    return summary
