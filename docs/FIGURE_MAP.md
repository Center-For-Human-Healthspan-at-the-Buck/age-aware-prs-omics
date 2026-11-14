# Figure-to-script and output map

The scripts export individual panels rather than assembled multipanel figures.
PDF and PNG versions are written with the same filename stem.

## Main figures — `scripts/02_plot_main_figures.R`

| Manuscript item | Output stem |
|---|---|
| Figure 1 | `fig01_workflow_schematic` |
| Figure 2a | `fig02a_age_residualized_chord` |
| Figure 2b | `fig02b_age_included_chord` |
| Figure 2c | `fig02c_age_residualized_top_features` |
| Figure 2d | `fig02d_age_included_top_features` |
| Figure 3a | `fig03a_neurodegeneration_fingerprint` |
| Figure 3b | `fig03b_pd_feature_by_model_matrix` |
| Figure 3c | `fig03c_als_feature_by_model_matrix` |
| Figure 3b–c legend | `fig03bc_shared_legend` |
| Figure 4a | `fig04a_ad_lollipop` |
| Figure 4b | `fig04b_ad_age_context` |
| Figure 5a | `fig05a_t2d_lollipop` |
| Figure 5b | `fig05b_cvd_lollipop` |
| Figure 5c | `fig05c_cvd_recurring_associations` |
| Figure 6a workflow | `fig06a1_disease_association_workflow` |
| Figure 6a chemistry ancestry | `fig06a2_chemistry_ancestry` |
| Figure 6a proteomics ancestry | `fig06a3_proteomics_ancestry` |
| Figure 6b clinical features | `fig06b1_clinical_directional` |
| Figure 6b proteins | `fig06b2_protein_directional` |
| Figure 6c clinical features | `fig06c1_clinical_discordant` |
| Figure 6c proteins | `fig06c2_protein_discordant` |
| Figure 6 legend | `fig06_shared_legend` |

## Supplementary figures — `scripts/03_plot_supplementary_figures.R`

| Manuscript item | Output stem |
|---|---|
| Supplementary Figure 1 | `supp_fig01_prs_cosine_age_residualized` |
| Supplementary Figure 2 | `supp_fig02_prs_cosine_age_included` |
| Supplementary Figure 3 | `supp_fig03_mantel_genetic_omic_dissimilarity` |
| Supplementary Figure 4 | `supp_fig04_ad_forest` |
| Supplementary Figure 5 | `supp_fig05_pd_forest` |
| Supplementary Figure 6 | `supp_fig06_als_forest` |
| Supplementary Figure 7 | `supp_fig07_t2d_forest` |
| Supplementary Figure 8 | `supp_fig08_cvd_forest` |
| Supplementary Figure 9 | `supp_fig09_longevity_forest` |

## Extended Data figures — `scripts/04_plot_extended_data_figures.R`

| Manuscript item | Output stem |
|---|---|
| Extended Data Figure 1a | `extended_data_fig01a_domain_model_counts` |
| Extended Data Figure 1b | `extended_data_fig01b_pairwise_overlap` |
| Extended Data Figure 1c | `extended_data_fig01c_age_context` |
| Extended Data Figure 2a | `extended_data_fig02a_positive_controls_age_residualized` |
| Extended Data Figure 2b | `extended_data_fig02b_positive_controls_age_included` |
| Extended Data Figure 3, longevity | `extended_data_fig03_longevity_lollipop` |
| Extended Data Figure 4a, AoU concordance | `extended_data_fig04a_aou_coverage_concordance` |
| Extended Data Figure 4b, proteins | `extended_data_fig04b_protein_domains` |
| Extended Data Figure 4c, chemistries | `extended_data_fig04c_chemistry_domains` |

The assembled Extended Data Figure 4 PDF and PNG are generated directly from
the aggregate AoU input by `scripts/05_assemble_extended_data_figure_4.py`;
the underlying calculated coverage and domain-count tables are exported beside
the figure in `figures/`.

Panel-level TSV files are written beside the relevant plots where the script
constructs a transformed plotting table.
