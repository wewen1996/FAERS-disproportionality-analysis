# Source this file to plot a completed analysis without reloading FAERS tables.
CODE_DIR <- "F:/GLP1-RA/new"
RUN_DIR <- file.path(CODE_DIR, "05_\u7ed3\u679c", "new_reanalysis", "run_20261019_180238_21536")
source(file.path(CODE_DIR, "07_additional_figures.R"), encoding = "UTF-8")
additional_figures_output <- plot_additional_results(RUN_DIR)
