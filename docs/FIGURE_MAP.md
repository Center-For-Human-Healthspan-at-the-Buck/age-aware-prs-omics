# Figure-to-script and output map

The R scripts export individual vector PDF panels and 300-dpi PNG previews.
`scripts/05_assemble_multipanel_figures.py` places the panel PDFs into final
multi-panel PDFs without rasterization. Extended Data Figure 4 is rendered by
`scripts/05_assemble_extended_data_figure_4.py`.

## Final assembled PDFs — `figures/`

| Manuscript item | Output file |
|---|---|
| Figure 1 | `Fig_1.pdf` |
| Figure 2 | `Fig_2.pdf` |
| Figure 3 | `Fig_3.pdf` |
| Figure 4 | `Fig_4.pdf` |
| Figure 5 | `Fig_5.pdf` |
| Figure 6 | `Fig_6.pdf` |
| Extended Data Figure 1 | `Extended_Data_Fig1_Domain_level_model_comparison.pdf` |
| Extended Data Figure 2 | `Extended_Data_Fig2_Vae_Positive_Control_PRS_Associations.pdf` |
| Extended Data Figure 3 | `Extended_Data_Fig3_VAE_longevity_associations.pdf` |
| Extended Data Figure 4 | `Extended_Data_Fig4_AoU_directional_concordance.pdf` |

Assembly order follows the manuscript artwork: Figure 2 uses a/b above c and
d; Figure 3 uses a above b/c and their shared legend; Figure 4 stacks a above
b and its separate legend; Figure 5 uses a/b above c; Figure 6 uses the workflow and ancestry panels
above the concordant and discordant forest-plot rows. Extended Data Figure 1
uses a/b above c, Extended Data Figure 2 uses a/b side by side, and Extended
Data Figure 4 uses a above b/c and the shared disease-domain legend.

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
| Figure 4b legend | `fig04b_legend` |
| Figure 5a | `fig05a_t2d_lollipop` |
| Figure 5b | `fig05b_cvd_lollipop` |
| Figure 5c | `fig05c_cvd_recurring_associations` |
| Figure 6a workflow | `fig06a1_disease_association_workflow` |
| Figure 6a chemistry ancestry | `fig06a2_chemistry_ancestry` |
| Figure 6a proteomics ancestry | `fig06a3_proteomics_ancestry` |
| Figure 6a shared ancestry legend | `fig06a_ancestry_legend` |
| Figure 6b clinical features | `fig06b1_clinical_directional` |
| Figure 6b proteins | `fig06b2_protein_directional` |
| Figure 6c clinical features | `fig06c1_clinical_discordant` |
| Figure 6c proteins | `fig06c2_protein_discordant` |
| Figure 6 reference legend (separate; not placed in the supplied composite layout) | `fig06_shared_legend` |

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

Extended Data Figure 4 is generated directly from the aggregate AoU input; the
underlying coverage and domain-count tables are exported beside the figure.

Panel-level TSV files are written beside the relevant plots where the script
constructs a transformed plotting table.
