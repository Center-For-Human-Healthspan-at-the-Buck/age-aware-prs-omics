.PHONY: install-r check linear-example main supplementary extended extended-assembled all clean

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

main:
	Rscript scripts/02_plot_main_figures.R

supplementary:
	Rscript scripts/03_plot_supplementary_figures.R

extended:
	Rscript scripts/04_plot_extended_data_figures.R

extended-assembled:
	python scripts/05_assemble_extended_data_figure_4.py

all: check main supplementary extended extended-assembled

clean:
	rm -f outputs/main/* outputs/supplementary/* outputs/extended_data/* outputs/tables/*
	rm -rf outputs/example_linear_models
