# Run with Rscript, or source this file and call run_conflict_primaryid_sensitivity().
# Standalone, frozen-cohort sensitivity: no fallback to an older report version.
# Historical result directories and raw inputs are never changed.

conflict_source_quarter <- function(value) {
  text <- toupper(as.character(value))
  tokens <- regmatches(text, regexec("([0-9]{2,4})Q([1-4])", text))
  vapply(tokens, function(x) {
    if (length(x) != 3L) return(NA_integer_)
    year <- as.integer(x[2L]); if (year < 100L) year <- year + 2000L
    if (year < 2004L || year > 2099L) return(NA_integer_)
    year * 4L + as.integer(x[3L]) - 1L
  }, integer(1))
}

conflict_compare <- function(before, after, keys, scenario) {
  fields <- c("a", "b", "c", "d", "total_reports", "ROR", "ROR_low", "ROR_high",
              "PRR", "PRR_low", "PRR_high", "ROR_signal", "PRR_signal", "fisher_p", "fisher_q_BH", "status")
  result <- merge(before[, c(keys, fields), with = FALSE], after[, c(keys, fields), with = FALSE],
    by = keys, suffixes = c("_before", "_after"), all = TRUE)
  stopifnot(nrow(result) == nrow(before), nrow(result) == nrow(after), !anyNA(result$a_after))
  result[, `:=`(scenario = scenario, removed_a = a_before - a_after, removed_b = b_before - b_after,
    removed_c = c_before - c_after, removed_d = d_before - d_after,
    ROR_percent_change = 100 * (ROR_after / ROR_before - 1),
    ROR_signal_changed = ROR_signal_before != ROR_signal_after,
    PRR_signal_changed = PRR_signal_before != PRR_signal_after)]
  stopifnot(all(result$removed_a >= 0), all(result$removed_b >= 0),
    all(result$removed_c >= 0), all(result$removed_d >= 0))
  result
}

run_conflict_primaryid_sensitivity <- function(baseline_dir = NULL, project_root = "F:/GLP1-RA/new", output_root = NULL) {
  suppressPackageStartupMessages(library(data.table))
  old_threads <- getDTthreads(); setDTthreads(4L)
  on.exit(setDTthreads(old_threads), add = TRUE)
  result_root <- file.path(project_root, "05_\u7ed3\u679c")
  if (is.null(baseline_dir)) {
    candidates <- sort(list.dirs(file.path(result_root, "new_reanalysis"), recursive = FALSE), decreasing = TRUE)
    complete <- vapply(candidates, function(p) {
      status <- file.path(p, "RUN_STATUS.txt")
      file.exists(status) && identical(trimws(readLines(status, warn = FALSE)[1L]), "COMPLETED")
    }, logical(1))
    if (!any(complete)) stop("No COMPLETED baseline run is available.")
    baseline_dir <- candidates[complete][1L]
  }
  stopifnot(identical(trimws(readLines(file.path(baseline_dir, "RUN_STATUS.txt"), warn = FALSE)[1L]), "COMPLETED"))
  # Freeze baseline numerical and drug-matching definitions. The corrected current
  # source table is not needed here; baseline definition metadata are not cited.
  for (script in c("01_helpers.R", "06_additional_analyses.R"))
    sys.source(file.path(baseline_dir, "code", script), envir = environment())
  sys.source(file.path(project_root, "11_additional_definitions.R"), envir = environment())
  config <- readRDS(file.path(baseline_dir, "analysis_config.rds"))
  definitions <- additional_definitions()
  if (is.null(output_root)) output_root <- file.path(result_root, "conflict_primaryid_sensitivity")
  out <- file.path(output_root,
    paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S", tz = "Asia/Shanghai"), "_", Sys.getpid()))
  stopifnot(!dir.exists(out))
  for (folder in c("tables", "audit", "code")) dir.create(file.path(out, folder), recursive = TRUE)
  completed <- FALSE
  writeLines("RUNNING", file.path(out, "RUN_STATUS.txt"))
  on.exit(if (!completed) writeLines("FAILED_OR_INTERRUPTED", file.path(out, "RUN_STATUS.txt")), add = TRUE)
  log <- file(file.path(out, "run_log.txt"), "wt")
  sink(log, split = TRUE); on.exit({sink(); close(log)}, add = TRUE)
  emit <- function(x, name, detail = FALSE) {
    fwrite(x, file.path(out, if (detail) "audit" else "tables", paste0(name, ".csv")), bom = TRUE, na = "")
    invisible(x)
  }
  read_detail <- function(name, select = NULL) fread(file.path(baseline_dir, "data", paste0(name, ".csv")),
    select = select, colClasses = c(PRIMARYID = "character"), showProgress = FALSE, encoding = "UTF-8")
  cat("Baseline:", baseline_dir, "\nOutput:", out, "\n")
  writeLines(c(paste("baseline_dir =", baseline_dir), "Primary sensitivity: remove every conflicted PRIMARYID from the frozen final universe.",
    "Conservative sensitivity: remove all final reports whose CASEID appears in any conflicted DEMO record.",
    "No older report is substituted; all other definitions remain those of the baseline.",
    "Active comparators: all configured targets/comparators and all four existing indication strata.",
    "Indication strata remain report-level, not drug-linked; this analysis does not fix indication confounding.",
    "ROR signal: a >= 3 and lower 95% CI > 1. PRR signal: a >= 3, PRR >= 2 and Pearson chi-square >= 4.",
    "Fisher BH: eight main rows per scenario; active comparisons within target drug and stratum per scenario."),
    file.path(out, "ANALYSIS_PLAN.txt"))
  snapshot <- c(file.path(project_root, "09_conflict_primaryid_sensitivity.R"),
    file.path(project_root, "11_additional_definitions.R"),
    file.path(baseline_dir, "code", c("01_helpers.R", "06_additional_analyses.R")))
  stopifnot(all(file.copy(snapshot, file.path(out, "code"))))
  reference_files <- file.path(baseline_dir, "tables", c("signals_main.csv", "active_comparator_signals.csv", "dedup_primaryid_case_conflicts.csv"))
  emit(data.table(path = c(snapshot, reference_files), md5 = unname(tools::md5sum(c(snapshot, reference_files)))), "code_reference_checksums")
  manifest <- fread(file.path(baseline_dir, "tables/input_manifest.csv"), encoding = "UTF-8")
  raw_paths <- file.path(config$data_dir, paste0(c("rawDEMO", "DRUG", "REAC", "INDI", "THER", "OUTC"), ".RData"))
  stopifnot(all(file.exists(raw_paths)))
  info <- file.info(raw_paths)
  source_manifest <- data.table(path = raw_paths, size_bytes = info$size,
    modified_at = format(info$mtime, tz = "Asia/Shanghai", usetz = TRUE))
  prior <- manifest[match(raw_paths, path)]
  stopifnot(!anyNA(prior$size_bytes), identical(as.numeric(prior$size_bytes), as.numeric(info$size)),
    identical(prior$modified_at, source_manifest$modified_at))
  emit(source_manifest, "raw_input_manifest")
  conflicts <- fread(reference_files[3L], colClasses = c(PRIMARYID = "character"), encoding = "UTF-8")
  ids <- conflicts$PRIMARYID
  stopifnot(length(ids) > 0L, !anyDuplicated(ids))
  emit(conflicts, "conflict_id_list")

  # Read each raw table once. Preserve FILE and drug sequence in the audit extracts.
  raw_subset <- function(name) {
    cat("Reading raw", name, "at", format(Sys.time()), "\n"); flush.console()
    e <- new.env(parent = baseenv())
    loaded <- load(file.path(config$data_dir, paste0(name, ".RData")), envir = e)
    choices <- loaded[vapply(loaded, function(n) is.data.frame(e[[n]]), logical(1))]
    selected <- if (name %in% choices) name else choices
    stopifnot(length(selected) == 1L)
    d <- as.data.table(e[[selected]]); setnames(d, toupper(names(d)))
    stopifnot("PRIMARYID" %in% names(d), !anyDuplicated(names(d)))
    emit(data.table(column = names(d), storage_class = vapply(d, function(x) class(x)[1L], character(1)), source_rows = nrow(d)),
      paste0("schema_", name), TRUE)
    if (is.numeric(d$PRIMARYID)) {
      numeric_ids <- suppressWarnings(as.numeric(ids))
      sub <- d[PRIMARYID %in% numeric_ids[!is.na(numeric_ids)]]
    } else sub <- d[PRIMARYID %chin% ids]
    sub[, PRIMARYID := identifier_text(PRIMARYID)]
    if ("CASEID" %in% names(sub)) sub[, CASEID := identifier_text(CASEID)]
    rm(d, e); gc(verbose = FALSE)
    sub
  }
  demo <- raw_subset("rawDEMO")
  demo[, SOURCE_QUARTER := quarter_index(GETDATAYEAR)]
  demo <- demo[!is.na(SOURCE_QUARTER) & SOURCE_QUARTER >= config$first_quarter & SOURCE_QUARTER <= config$last_quarter]
  discovered <- demo[, .(case_count = uniqueN(CASEID)), by = PRIMARYID][case_count > 1L]
  stopifnot(setequal(discovered$PRIMARYID, ids))
  implicated_cases <- unique(demo$CASEID)
  emit(demo, "conflicting_DEMO_rows", TRUE)
  identity <- unique(demo[, .(PRIMARYID, SOURCE_QUARTER, CASEID)])
  identity_summary <- identity[, .(source_CASEID_count = .N, source_CASEIDs = paste(sort(CASEID), collapse = " | "),
    unique_source_CASEID = if (.N == 1L) CASEID else NA_character_), by = .(PRIMARYID, SOURCE_QUARTER)]
  emit(identity_summary, "quarter_case_identity", TRUE)
  cat("Reading frozen analysis core at", format(Sys.time()), "\n")
  core <- as.data.table(readRDS(file.path(baseline_dir, "data/analysis_core.rds")))
  core <- core[, .(PRIMARYID, CASEID, QUARTER, TARGET)]
  stopifnot(!anyDuplicated(core$PRIMARYID), !anyDuplicated(core$CASEID), !anyNA(core$TARGET))
  baseline_conflict <- core[PRIMARYID %chin% ids]
  retained_identity <- baseline_conflict[, .(PRIMARYID, retained_CASEID = CASEID, retained_QUARTER = QUARTER)]
  emit(core[CASEID %chin% implicated_cases], "retained_implicated_case_reports", TRUE)
  raw <- list(rawDEMO = demo)
  source_audits <- list(); value_audits <- list()
  value_columns <- list(DRUG = c("DRUG_SEQ", "DRUGNAME", "PROD_AI", "ROLE_COD"), REAC = "PT",
    INDI = c("DRUG_SEQ", "INDI_PT"), THER = c("DRUG_SEQ", "START_DT", "END_DT"), OUTC = "OUTC_COD")
  for (name in c("DRUG", "REAC", "INDI", "THER", "OUTC")) {
    sub <- raw_subset(name)
    has_source <- "FILE" %in% names(sub)
    sub[, SOURCE_QUARTER := if (has_source) conflict_source_quarter(FILE) else NA_integer_]
    if ("DRUG_SEQ" %in% names(sub)) sub[, DRUG_SEQ := clean_text(DRUG_SEQ)]
    sub <- merge(sub, identity_summary, by = c("PRIMARYID", "SOURCE_QUARTER"), all.x = TRUE, sort = FALSE)
    sub <- merge(sub, retained_identity, by = "PRIMARYID", all.x = TRUE, sort = FALSE)
    sub[, `:=`(retained_in_baseline = !is.na(retained_CASEID),
      source_quarter_mismatch = !is.na(retained_QUARTER) & !is.na(SOURCE_QUARTER) & SOURCE_QUARTER != retained_QUARTER,
      source_CASEID_mismatch = !is.na(retained_CASEID) & !is.na(unique_source_CASEID) & unique_source_CASEID != retained_CASEID)]
    cols <- intersect(value_columns[[name]], names(sub))
    # Compare clinical field sets with the selected-quarter records. This is a
    # provenance discrepancy, not proof that the two source CASEIDs are two people.
    sub[, clinical_signature := do.call(paste, c(lapply(.SD, function(x) ifelse(is.na(x), "<NA>", as.character(x))), sep = "\034")), .SDcols = cols]
    signatures <- unique(sub[, .(PRIMARYID, SOURCE_QUARTER, clinical_signature)])
    variation <- signatures[, .(clinical_value_sets_differ = {
      sets <- lapply(split(clinical_signature, SOURCE_QUARTER), function(x) sort(unique(x)))
      length(sets) > 1L && !all(vapply(sets[-1L], identical, logical(1), sets[[1L]]))
    }), by = PRIMARYID]
    selected_values <- unique(sub[retained_in_baseline & !is.na(SOURCE_QUARTER) & SOURCE_QUARTER == retained_QUARTER,
      .(PRIMARYID, clinical_signature)])
    selected_values[, selected_quarter_value_present := TRUE]
    sub <- merge(sub, selected_values, by = c("PRIMARYID", "clinical_signature"), all.x = TRUE, sort = FALSE)
    sub[, other_source_only_value := source_quarter_mismatch & is.na(selected_quarter_value_present)]
    audit <- sub[, .(raw_rows = .N, source_available = has_source,
      source_quarters = paste(sort(unique(quarter_label(SOURCE_QUARTER[!is.na(SOURCE_QUARTER)]))), collapse = " | "),
      source_quarter_count = uniqueN(SOURCE_QUARTER[!is.na(SOURCE_QUARTER)]),
      unparsed_source_rows = sum(is.na(SOURCE_QUARTER)),
      rows_without_quarter_case_match = sum(is.na(source_CASEID_count)),
      retained_in_baseline = any(retained_in_baseline),
      source_quarter_mismatch_rows = sum(source_quarter_mismatch),
      source_CASEID_mismatch_rows = sum(source_CASEID_mismatch),
      other_source_only_value_rows = sum(other_source_only_value)), by = PRIMARYID]
    audit <- merge(audit, variation, by = "PRIMARYID", all.x = TRUE)
    audit[, table := name]; source_audits[[name]] <- audit
    sub[, clinical_signature := NULL]
    raw[[name]] <- sub
    emit(sub, paste0("conflict_raw_", name), TRUE)
  }
  audit <- rbindlist(source_audits, fill = TRUE)
  emit(audit, "source_linkage_by_id")
  emit(audit[, .(conflict_ids_with_rows = .N, raw_rows = sum(raw_rows),
    source_available = all(source_available), multi_quarter_ids = sum(source_quarter_count > 1L),
    clinical_value_sets_differ_ids = sum(clinical_value_sets_differ),
    retained_ids_with_rows = sum(retained_in_baseline),
    retained_ids_with_other_quarter_rows = sum(source_quarter_mismatch_rows > 0),
    retained_ids_with_other_CASEID_rows = sum(source_CASEID_mismatch_rows > 0),
    retained_ids_with_other_source_only_values = sum(other_source_only_value_rows > 0),
    source_quarter_mismatch_rows = sum(source_quarter_mismatch_rows),
    source_CASEID_mismatch_rows = sum(source_CASEID_mismatch_rows), unparsed_source_rows = sum(unparsed_source_rows)), by = table], "source_linkage_summary")

  cat("Reading frozen drug exposures and indication flags at", format(Sys.time()), "\n")
  pairs <- read_detail("primary_exposure_pairs")
  matched <- read_detail("matched_drug_records")
  hits <- read_detail("additional_drug_name_matches", c("PRIMARYID", "ingredient", "ROLE_COD", "combination_flag"))
  hits <- unique(hits[ROLE_COD == "PS" & !combination_flag, .(PRIMARYID, ingredient)])
  flags <- read_detail("additional_indication_flags", c("PRIMARYID", "diabetes_recorded", "weight_management_recorded", "metabolic_or_weight_recorded"))
  stopifnot(!anyDuplicated(flags$PRIMARYID), setequal(flags$PRIMARYID, core$PRIMARYID))
  mapping <- fread(file.path(baseline_dir, "data/mapping_snapshot.csv"), encoding = "UTF-8")
  raw_drug <- copy(raw$DRUG)
  raw_drug[, ROLE_COD := toupper(clean_text(ROLE_COD))]
  map_index <- match(raw_drug$DRUGNAME, mapping$original_drug_name)
  raw_drug[, `:=`(study_drug = mapping$standardized_ingredient[map_index], main = mapping$include_main[map_index])]
  raw_study <- raw_drug[!is.na(main) & main & !is.na(study_drug)]
  raw_target <- raw$REAC[toupper(clean_text(PT)) == config$target_pt]
  raw_groups <- c(list(study_drug_group = unique(raw_study$PRIMARYID)),
    setNames(lapply(config$drugs, function(g) unique(raw_study[study_drug == g, PRIMARYID])), config$drugs))
  group_impact <- rbindlist(lapply(names(raw_groups), function(g) {
    ps_ids <- if (g == "study_drug_group") unique(raw_study[ROLE_COD == "PS", PRIMARYID]) else unique(raw_study[study_drug == g & ROLE_COD == "PS", PRIMARYID])
    retained_ids <- if (g == "study_drug_group") unique(pairs$PRIMARYID) else unique(pairs[drug == g, PRIMARYID])
    data.table(drug = g, raw_conflict_any_role_ids = length(raw_groups[[g]]), raw_conflict_PS_ids = length(ps_ids),
      raw_conflict_PS_target_ids = length(intersect(ps_ids, raw_target$PRIMARYID)),
      baseline_conflict_PS_reports = sum(baseline_conflict$PRIMARYID %chin% retained_ids),
      baseline_conflict_PS_target_reports = sum(baseline_conflict$TARGET & baseline_conflict$PRIMARYID %chin% retained_ids))
  }))
  emit(group_impact, "conflict_study_drug_impact")
  case_removed <- core[CASEID %chin% implicated_cases, PRIMARYID]
  exclusions <- list(exclude_conflict_PRIMARYID = intersect(ids, core$PRIMARYID), exclude_implicated_CASEID = case_removed)
  impact <- rbindlist(lapply(names(exclusions), function(s) {
    removed <- core[PRIMARYID %chin% exclusions[[s]]]
    exposed <- removed$PRIMARYID %chin% pairs$PRIMARYID
    data.table(scenario = s, excluded_reports = nrow(removed), excluded_target_reports = sum(removed$TARGET),
      excluded_study_PS_reports = sum(exposed), excluded_study_PS_target_reports = sum(exposed & removed$TARGET),
      remaining_reports = nrow(core) - nrow(removed), remaining_target_reports = sum(core$TARGET) - sum(removed$TARGET))
  }))
  emit(impact, "exclusion_population_summary")
  emit(data.table(metric = c("conflict_PRIMARYIDs", "raw_conflicting_DEMO_rows", "implicated_CASEIDs",
    "raw_conflict_target_PRIMARYIDs", "raw_conflict_target_REAC_rows", "baseline_retained_conflict_PRIMARYIDs",
    "baseline_not_retained_conflict_PRIMARYIDs"), value = c(length(ids), nrow(demo), length(implicated_cases),
      uniqueN(raw_target$PRIMARYID), nrow(raw_target), nrow(baseline_conflict), length(setdiff(ids, core$PRIMARYID)))), "audit_population_counts")

  # Expose target-specific source origin without treating a missing source as absence.
  target_origin <- raw_target[, .(PRIMARYID, PT, FILE, SOURCE_QUARTER, unique_source_CASEID,
    retained_CASEID, retained_QUARTER, retained_in_baseline, source_CASEID_mismatch, other_source_only_value)]
  emit(target_origin, "target_event_source_audit", TRUE)
  emit(raw_study, "study_drug_source_audit", TRUE)
  # Therapy/indication joins also use DRUG_SEQ in the baseline. Count only
  # cross-quarter pairs for the same ID and sequence; source-less rows are unknown.
  seq_links <- list()
  for (name in c("INDI", "THER")) {
    drug_keys <- unique(raw$DRUG[!is.na(SOURCE_QUARTER) & nzchar(DRUG_SEQ), .(PRIMARYID, DRUG_SEQ, drug_quarter = SOURCE_QUARTER)])
    child_keys <- unique(raw[[name]][!is.na(SOURCE_QUARTER) & nzchar(DRUG_SEQ), .(PRIMARYID, DRUG_SEQ, child_quarter = SOURCE_QUARTER)])
    linked <- merge(drug_keys, child_keys, by = c("PRIMARYID", "DRUG_SEQ"), allow.cartesian = TRUE)
    linked[, `:=`(table = name, cross_quarter_pair = drug_quarter != child_quarter,
      retained_in_baseline = PRIMARYID %chin% core$PRIMARYID,
      baseline_study_PS = PRIMARYID %chin% pairs$PRIMARYID)]
    seq_links[[name]] <- linked
  }
  emit(rbindlist(seq_links), "drug_sequence_source_join_pairs", TRUE)

  comparator_ids <- lapply(definitions$comparators, function(ingredients) unique(hits[ingredient %chin% ingredients, PRIMARYID]))
  target_ids <- setNames(lapply(config$active_comparator_targets, function(g) unique(matched[main & ROLE_COD == "PS" & drug == g, PRIMARYID])), config$active_comparator_targets)
  all_exposure_ids <- unique(c(unlist(comparator_ids, use.names = FALSE), unlist(target_ids, use.names = FALSE)))
  active_results <- function(removed_ids) {
    eligible_core <- core[PRIMARYID %chin% all_exposure_ids & !PRIMARYID %chin% removed_ids]
    relevant_flags <- flags[PRIMARYID %chin% eligible_core$PRIMARYID]
    rows <- list()
    for (g in names(target_ids)) for (comparator in names(comparator_ids)) {
      for (stratum in c("all_reports", "diabetes_recorded", "weight_management_recorded", "metabolic_or_weight_recorded")) {
        row <- active_comparator_one(eligible_core, relevant_flags, target_ids[[g]], comparator_ids[[comparator]],
          g, comparator, definitions$comparators[[comparator]], stratum)
        row[, eligible_background_reports := if (stratum == "all_reports") nrow(core) - length(removed_ids) else
          flags[!PRIMARYID %chin% removed_ids & get(stratum), .N]]
        rows[[length(rows) + 1L]] <- row
      }
    }
    result <- rbindlist(rows)
    result[, fisher_q_BH := p.adjust(fisher_p, "BH"), by = .(target_drug, stratum)]
    result
  }
  main_results <- function(removed_ids) {
    result <- signals_for_groups(core[!PRIMARYID %chin% removed_ids], pairs, config$drugs, "main_PS")
    result[, fisher_q_BH := p.adjust(fisher_p, "BH")]
    result
  }
  check_baseline <- function(actual, file, keys) {
    expected <- fread(file, encoding = "UTF-8")
    metrics <- c("a", "b", "c", "d", "ROR", "ROR_low", "ROR_high", "PRR", "fisher_p", "fisher_q_BH", "ROR_signal", "PRR_signal")
    joined <- merge(actual, expected, by = keys, suffixes = c("_actual", "_expected"))
    stopifnot(nrow(joined) == nrow(actual), nrow(actual) == nrow(expected))
    check <- rbindlist(lapply(metrics, function(m) data.table(table = basename(file), metric = m,
      pass = isTRUE(all.equal(joined[[paste0(m, "_actual")]], joined[[paste0(m, "_expected")]], tolerance = 1e-10, check.attributes = FALSE)))))
    stopifnot(all(check$pass)); check
  }
  cat("Reconstructing baseline tables at", format(Sys.time()), "\n")
  before_main <- main_results(character()); before_active <- active_results(character())
  main_keys <- c("drug", "stratum"); active_keys <- c("target_drug", "comparator", "stratum")
  checks <- rbindlist(list(check_baseline(before_main, reference_files[1L], main_keys),
    check_baseline(before_active, reference_files[2L], active_keys)))
  emit(checks, "baseline_reproduction_checks")
  emit(before_main, "main_before"); emit(before_active, "active_comparator_before")
  main_comparisons <- list(); active_comparisons <- list()
  for (scenario in names(exclusions)) {
    cat("Calculating", scenario, "at", format(Sys.time()), "\n")
    removed_ids <- exclusions[[scenario]]
    after_main <- main_results(removed_ids); after_active <- active_results(removed_ids)
    main_comparisons[[scenario]] <- conflict_compare(before_main, after_main, main_keys, scenario)
    active_comparisons[[scenario]] <- conflict_compare(before_active, after_active, active_keys, scenario)
    emit(after_main, paste0("main_", scenario)); emit(after_active, paste0("active_", scenario))
    emit(core[PRIMARYID %chin% removed_ids], paste0("excluded_reports_", scenario), TRUE)
  }
  main_comparison <- rbindlist(main_comparisons); active_comparison <- rbindlist(active_comparisons)
  emit(main_comparison, "main_before_after_comparison")
  emit(active_comparison, "active_comparator_before_after_comparison")
  combined <- rbindlist(list(main_comparison[, .(scenario, family = "main", ROR_signal_changed, PRR_signal_changed, ROR_percent_change)],
    active_comparison[, .(scenario, family = "active_comparator", ROR_signal_changed, PRR_signal_changed, ROR_percent_change)]))
  emit(combined[, .(comparisons = .N, ROR_signal_changes = sum(ROR_signal_changed), PRR_signal_changes = sum(PRR_signal_changed),
    max_absolute_ROR_percent_change = max(abs(ROR_percent_change), na.rm = TRUE)), by = .(scenario, family)], "signal_change_summary")
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  completed <- TRUE; writeLines("COMPLETED", file.path(out, "RUN_STATUS.txt"))
  cat("COMPLETED:", out, "\n")
  invisible(out)
}

if (sys.nframe() == 0L) conflict_sensitivity_output <- run_conflict_primaryid_sensitivity()
