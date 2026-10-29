PROJECT_ROOT <- "F:/GLP1-RA/new"
CODE_DIR <- PROJECT_ROOT

CONFIG <- list(
  project_root = PROJECT_ROOT,
  data_dir = file.path(PROJECT_ROOT, "raw", "2004Q1-2026Q2\u5df2\u6e05\u6d17\u6570\u636e", "Rdata\u96c6"),
  mapping_dir = file.path(PROJECT_ROOT, "Drug_name_mapping"),
  output_root = file.path(PROJECT_ROOT, "05_\u7ed3\u679c", "new_reanalysis"),
  target_pt = "COLITIS ISCHAEMIC",
  drugs = c("albiglutide", "dulaglutide", "exenatide", "liraglutide", "lixisenatide", "semaglutide", "tirzepatide"),
  first_quarter_code = 41L,
  last_quarter_code = 262L,
  age_storage = "already_years",
  maximum_age = 120,
  minimum_weibull_n = 10L,
  minimum_comed_reports = 5L,
  make_figures = TRUE,
  active_comparator_targets = c("semaglutide", "tirzepatide"),
  notoriety_cutoffs = as.Date(c("2023-07-01", "2021-06-04", "2023-11-08"))
)

required_packages <- c("data.table", "survival", "openxlsx")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install missing packages in RStudio: install.packages(c(",
  paste(sprintf('"%s"', missing_packages), collapse = ", "), "))")
suppressPackageStartupMessages(library(data.table))
for (script in c("01_helpers.R", "11_additional_definitions.R", "02_build_cohort.R", "03_descriptive_and_signals.R",
                  "04_time_to_onset.R", "06_additional_analyses.R", "05_figures_and_flow.R", "07_additional_figures.R")) {
  source(file.path(CODE_DIR, script), encoding = "UTF-8")
}

run_reanalysis <- function(config) {
  stopifnot(identical(config$age_storage, "already_years"), length(config$drugs) > 0,
    !anyDuplicated(config$drugs), config$target_pt == "COLITIS ISCHAEMIC")
  config$first_quarter <- quarter_index(config$first_quarter_code)
  config$last_quarter <- quarter_index(config$last_quarter_code)
  stopifnot(!is.na(config$first_quarter), !is.na(config$last_quarter), config$first_quarter <= config$last_quarter)
  last_start <- as.Date(sprintf("%04d-%02d-01", config$last_quarter %/% 4L, (config$last_quarter %% 4L) * 3L + 1L))
  config$study_end <- seq(last_start, by = "3 months", length.out = 2L)[2L] - 1
  input_paths <- c(file.path(config$data_dir, paste0(c("rawDEMO", "DRUG", "REAC", "OUTC", "THER", "INDI"), ".RData")),
    file.path(config$mapping_dir, c("drug_name_mapping.rds", "drug_lists_for_analysis.rds")))
  if (any(!file.exists(input_paths))) stop("Missing inputs: ", paste(input_paths[!file.exists(input_paths)], collapse = "; "))
  context <- new.env(parent = emptyenv())
  context$config <- config
  context$tables <- list()
  context$flow <- data.table::data.table()
  context$schema <- data.table::data.table()
  run_name <- paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S", tz = "Asia/Shanghai"), "_", Sys.getpid())
  context$out <- file.path(config$output_root, run_name)
  if (dir.exists(context$out)) stop("Output directory already exists; run again with a new timestamp.")
  for (directory in c("tables", "data", "figures", "code")) dir.create(file.path(context$out, directory), recursive = TRUE)
  status_path <- file.path(context$out, "RUN_STATUS.txt")
  writeLines("RUNNING", status_path)
  completed <- FALSE
  on.exit(if (!completed) writeLines("FAILED_OR_INTERRUPTED: inspect run_log.txt; do not use as complete results.", status_path), add = TRUE)
  log_file <- file(file.path(context$out, "run_log.txt"), open = "wt", encoding = "UTF-8")
  sink(log_file, split = TRUE)
  on.exit({sink(); close(log_file)}, add = TRUE)
  saveRDS(config, file.path(context$out, "analysis_config.rds"))
  dput(config, file = file.path(context$out, "analysis_config.txt"))
  info <- file.info(input_paths)
  write_result(context, data.table::data.table(path = input_paths, size_bytes = info$size,
    modified_at = format(info$mtime, tz = "Asia/Shanghai", usetz = TRUE)), "input_manifest")
  code_files <- file.path(CODE_DIR, c("00_run_reanalysis.R", "01_helpers.R", "11_additional_definitions.R",
    "02_build_cohort.R", "03_descriptive_and_signals.R", "04_time_to_onset.R", "06_additional_analyses.R",
    "05_figures_and_flow.R", "07_additional_figures.R", "08_plot_existing_results.R", "99_self_check.R", "README.md"))
  if (!all(file.exists(code_files))) stop("The code folder is incomplete; retain all supplied scripts and README.md.")
  if (!all(file.copy(code_files, file.path(context$out, "code")))) stop("Code snapshot failed.")
  small_inputs <- c(code_files, tail(input_paths, 2L))
  write_result(context, data.table::data.table(path = small_inputs,
    md5 = unname(tools::md5sum(small_inputs))), "code_mapping_checksums")
  cat("Output:", context$out, "\nTarget:", config$target_pt, "\n")
  tryCatch({
    cohort <- build_cohort(context)
    run_descriptive_and_signals(context, cohort)
    tto <- run_time_to_onset(context, cohort)
    run_comedication(context, cohort)
    run_additional_analyses(context, cohort)
    write_flow_and_figures(context, cohort, tto)
    write_additional_figures(context)
    workbook_tables <- context$tables[vapply(context$tables, nrow, integer(1)) <= 1048575L]
    names(workbook_tables) <- make.unique(substr(names(workbook_tables), 1L, 28L), sep = "_")
    openxlsx::write.xlsx(workbook_tables, file.path(context$out, "analysis_summary.xlsx"), overwrite = FALSE)
    capture.output(sessionInfo(), file = file.path(context$out, "sessionInfo.txt"))
  }, error = function(condition) {
    cat("FAILED:", conditionMessage(condition), "\n")
    stop(condition)
  })
  completed <- TRUE
  writeLines("COMPLETED", status_path)
  cat("Completed. Results:", context$out, "\n")
  invisible(context$out)
}

reanalysis_output <- run_reanalysis(CONFIG)
