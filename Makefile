# Makefile — Madagascar Influenza Analysis
# Activate the conda environment first:
#   conda activate madagascar_influenza

PYTHON   = conda run -n madagascar_influenza python
QUARTO   = conda run -n madagascar_influenza quarto
DATA_DIR = data

.PHONY: all analyse report clean

all: analyse report

analyse:
	$(PYTHON) scripts/analyse_gisaid.py

report: $(DATA_DIR)/combined_metadata.tsv
	cd report && $(QUARTO) render report.qmd --output-dir ../docs

clean:
	rm -f $(DATA_DIR)/madagascar_clean.tsv $(DATA_DIR)/avian_africa_clean.tsv
	rm -f docs/index.html
	rm -rf report/.quarto report/report_files
