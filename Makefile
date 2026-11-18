.PHONY: install-r check linear-example main-panels main supplementary \
	extended-panels extended extended-assembled assembled validate all clean

install-r:
	Rscript scripts/install_r_dependencies.R

check:
	Rscript scripts/00_check_inputs.R

linear-example:
	Rscript scripts/01_fit_linear_models.R \
		--matrix data/templates/example_analysis_matrix.tsv \
		--features data/templates/feature_manifest.tsv \
		--prs data/templates/prs_manifest.tsv \
		--age-column age \
		--covariates sex,PC1,PC2 \
		--min-n 100 \
		--cores 1 \
		--output-dir outputs/example_linear_models

main-panels:
	Rscript scripts/02_plot_main_figures.R

main: main-panels
	python3 scripts/05_assemble_multipanel_figures.py --main-only

supplementary:
	Rscript scripts/03_plot_supplementary_figures.R

extended-panels:
	Rscript scripts/04_plot_extended_data_figures.R

extended: extended-panels
	python3 scripts/05_assemble_multipanel_figures.py --extended-only
	python3 scripts/05_assemble_extended_data_figure_4.py

extended-assembled: extended-panels
	python3 scripts/05_assemble_multipanel_figures.py --extended-only
	python3 scripts/05_assemble_extended_data_figure_4.py

assembled: main-panels extended-panels
	python3 scripts/05_assemble_multipanel_figures.py
	python3 scripts/05_assemble_extended_data_figure_4.py

validate:
	python3 scripts/06_validate_package.py

all: check main supplementary extended

clean:
	rm -f outputs/main/* outputs/supplementary/* outputs/extended_data/* outputs/tables/*
	rm -f figures/Fig_*.pdf figures/Extended_Data_*.pdf figures/Extended_Data_*.png
	rm -f figures/Extended_Data_*_coverage.tsv figures/Extended_Data_*_domains.tsv
	rm -rf outputs/example_linear_models
