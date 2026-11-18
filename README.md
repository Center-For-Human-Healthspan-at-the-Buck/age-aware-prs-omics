# Age-informed PRS–omics analysis

Analysis and figure-generation code for:

> Nogal et al. *Age-informed modeling of multi-omic manifestations of polygenic risk in a perturbable latent space.*

The repository contains:

- figure-generation code for six main figures, nine Supplementary Figures and four Extended Data Figures;
- aggregate VAE, OLS and All of Us results used by the figures;
- an analyte-wise implementation of the OLS main-effect, linear PRS × age and PRS × spline(age) models;
- the All of Us Controlled Tier analysis workflow.

The Arivale discovery cohort contained 1,491 participants. VAE training is not
rerun; figures are generated from frozen aggregate association results. The VAE
analyses used MOVE v1.5.0.

## Requirements

- R 4.x
- Python 3.9 or later

Install the R dependencies:

```bash
Rscript scripts/install_r_dependencies.R
```

Install the Python dependencies:

```bash
python3 -m pip install -r requirements.txt
```

## Reproduce the figures

Validate the aggregate inputs:

```bash
Rscript scripts/00_check_inputs.R
```

Generate the panel files and assembled figures:

```bash
Rscript scripts/02_plot_main_figures.R
Rscript scripts/03_plot_supplementary_figures.R
Rscript scripts/04_plot_extended_data_figures.R
python3 scripts/05_assemble_multipanel_figures.py
python3 scripts/05_assemble_extended_data_figure_4.py
```

The equivalent command is:

```bash
make all
```

Outputs are written to:

- `outputs/main/`
- `outputs/supplementary/`
- `outputs/extended_data/`
- `figures/` (assembled Main and Extended Data PDFs)

Panel-level PDF and 300-dpi PNG files are produced. Final multi-panel PDFs are
assembled directly from the panel PDFs without flattening them to raster
images. Main figures are limited to 183 mm wide by 170 mm tall; Extended Data
figures are limited to 180 mm wide by 170 mm tall. See `docs/FIGURE_MAP.md` for
the complete panel-to-figure mapping.

## Linear-model workflow

`scripts/01_fit_linear_models.R` fits:

1. age-residualized OLS main effects;
2. age-included OLS main effects;
3. linear PRS × age interactions;
4. joint PRS × natural-spline(age) interactions.

A synthetic example is included:

```bash
make linear-example
```

The generic workflow uses pairwise complete observations, Benjamini–Hochberg
correction and a nested-model F test for spline interactions. Study-specific
covariates and model settings are passed as command-line arguments.

## Data organization

- `data/vae/`: aggregate VAE associations for age-residualized and age-included models
- `data/manuscript/`: compact OLS and model-comparison inputs
- `data/metadata/`: omics metadata and the PRS correlation matrix
- `data/aou/`: aggregate All of Us disease-association results and target lists
- `data/templates/`: synthetic inputs for the linear-model example
- `config/`: MOVE configuration files

The upstream All of Us filename `replication_summary.tsv` is retained for
compatibility. The analysis compares directional concordance between Arivale
discovery associations and All of Us disease–control contrasts; it is not a
same-estimand replication.

## Access restrictions

Participant-level Arivale data are not distributed because of consent and
data-use restrictions. Participant-level All of Us data are available only in
an authorized Researcher Workbench. The aggregate inputs included here are
sufficient to generate the figures.

The All of Us workflow is documented in `scripts/aou_workbench/README.md`.

## Repository structure

```text
R/                 plotting and data-processing functions
scripts/           validation, analysis and figure-generation scripts
config/            MOVE configuration files
data/              aggregate and synthetic inputs
assets/            static workflow schematics
outputs/           generated panel files
figures/           assembled Main and Extended Data figures
docs/              figure map
```

## Citation

Citation metadata are provided in `CITATION.cff`.
