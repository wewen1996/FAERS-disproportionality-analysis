# Compare TTO-available and TTO-unavailable target reports.
# The comparison is report-level (pooled study_drug_group rows), while drug
# membership is tabulated separately and may overlap within a report.

run_tto_availability_comparison <- function(
    baseline_dir,
    output_root = "F:/GLP1-RA/new/05_结果/tto_availability_comparison") {
  suppressPackageStartupMessages(library(data.table))
  if (.Platform$OS.type == "windows") {
    try(suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936")), silent = TRUE)
  }
  baseline_dir <- normalizePath(baseline_dir, winslash = "/")
  stopifnot(file.exists(file.path(baseline_dir, "RUN_STATUS.txt")))
  stopifnot(identical(trimws(readLines(file.path(baseline_dir, "RUN_STATUS.txt"), warn = FALSE)[1L]), "COMPLETED"))
  if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required.")
  helper <- file.path(baseline_dir, "code", "01_helpers.R")
  if (file.exists(helper)) sys.source(helper, envir = environment())
  if (!exists("summarize_numeric", inherits = FALSE) || !exists("weibull_fit", inherits = FALSE)) {
    stop("Baseline helper functions summarize_numeric/weibull_fit are unavailable.")
  }
  old_threads <- data.table::getDTthreads(); data.table::setDTthreads(4L)
  on.exit(data.table::setDTthreads(old_threads), add = TRUE)

  out <- file.path(output_root, paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
  dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "audit"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "code"), recursive = TRUE, showWarnings = FALSE)
  writeLines("RUNNING", file.path(out, "RUN_STATUS.txt"))
  complete <- FALSE
  on.exit(if (!complete) writeLines("FAILED_OR_INTERRUPTED", file.path(out, "RUN_STATUS.txt")), add = TRUE)
  log_con <- file(file.path(out, "run_log.txt"), open = "wt")
  sink(log_con, split = TRUE)
  on.exit({ sink(); close(log_con) }, add = TRUE)
  emit <- function(x, name, folder = "tables") {
    data.table::fwrite(x, file.path(out, folder, paste0(name, ".csv")), bom = TRUE, na = "")
    invisible(x)
  }

  tto_path <- file.path(baseline_dir, "data", "tto_case_level.csv")
  tto_all <- data.table::fread(tto_path, encoding = "UTF-8",
    colClasses = list(character = c("PRIMARYID", "drug", "CASEID", "EVENT_DT", "FDA_DT")))
  date_cols <- intersect(c("START_DATE", "EVENT_DATE", "FDA_DATE"), names(tto_all))
  for (nm in date_cols) tto_all[, (nm) := as.IDate(get(nm))]
  tto <- tto_all[drug == "study_drug_group"]
  stopifnot(!anyDuplicated(tto$PRIMARYID), nrow(tto) == 241L)
  tto[, TTO_AVAILABLE := status == "valid"]
  tto[, `:=`(FDA_YEAR = as.integer(format(FDA_DATE, "%Y")),
    availability_group = fifelse(status == "valid", "available", "unavailable"))]
  stopifnot(sum(tto$TTO_AVAILABLE) == 65L, sum(!tto$TTO_AVAILABLE) == 176L)

  # Recover raw start-date information from the individual drug rows. The
  # pooled output intentionally blanks START_DATE for an invalid pooled row.
  individual <- tto_all[drug != "study_drug_group"]
  start_info <- individual[, .(
    n_index_drug_rows = .N,
    n_drugs_with_start_date = sum(!is.na(START_DATE)),
    any_start_date = any(!is.na(START_DATE)),
    all_drugs_start_date = all(!is.na(START_DATE)),
    any_missing_sequence = any(n_missing_sequence > 0L),
    any_unmatched_therapy = any(n_unmatched_therapy > 0L),
    any_invalid_start = any(n_invalid_start > 0L),
    any_multiple_start = any(n_distinct_start > 1L)), by = PRIMARYID]
  tto <- merge(tto, start_info, by = "PRIMARYID", all.x = TRUE, sort = FALSE)
  for (nm in c("any_start_date", "all_drugs_start_date", "any_missing_sequence",
               "any_unmatched_therapy", "any_invalid_start", "any_multiple_start")) {
    tto[is.na(get(nm)), (nm) := FALSE]
  }
  tto[, `:=`(
    event_date_present = !is.na(EVENT_DATE),
    fda_date_present = !is.na(FDA_DATE),
    all_required_date_fields_present = !is.na(EVENT_DATE) & !is.na(FDA_DATE) & all_drugs_start_date)]
  tto[, `:=`(
    date_order_valid = all_required_date_fields_present & !is.na(TTO_DAYS) & TTO_DAYS >= 0 & EVENT_DATE <= FDA_DATE,
    date_information_complete_but_rejected = all_required_date_fields_present & !TTO_AVAILABLE)]

  # Demographic variables are taken from the completed baseline export.
  demo <- data.table::fread(file.path(baseline_dir, "data", "exposed_case_demographics.csv"), encoding = "UTF-8",
    colClasses = list(character = c("PRIMARYID", "CASEID", "FDA_DT", "EVENT_DT", "SEX", "SEX_GROUP", "AGE_GROUP", "REPORTER")))
  demo <- unique(demo[PRIMARYID %chin% tto$PRIMARYID,
    .(PRIMARYID, AGE_YEARS, AGE_GROUP, SEX_GROUP, REPORTER, FDA_DATE, TARGET)])
  stopifnot(!anyDuplicated(demo$PRIMARYID))
  tto <- merge(tto, demo, by = "PRIMARYID", all.x = TRUE, sort = FALSE, suffixes = c("", "_demo"))
  stopifnot(nrow(tto) == 241L)
  tto[, `:=`(
    AGE_YEARS = suppressWarnings(as.numeric(AGE_YEARS)),
    AGE_GROUP = fifelse(is.na(AGE_GROUP) | !nzchar(AGE_GROUP), "Unknown", as.character(AGE_GROUP)),
    SEX_GROUP = fifelse(is.na(SEX_GROUP) | !nzchar(SEX_GROUP), "Unknown", as.character(SEX_GROUP)),
    REPORTER = fifelse(is.na(REPORTER) | !nzchar(REPORTER), "Unknown", as.character(REPORTER)),
    FDA_YEAR = fifelse(is.na(FDA_YEAR), as.integer(format(FDA_DATE, "%Y")), FDA_YEAR))]
  tto[, FDA_YEAR_BAND := cut(FDA_YEAR, breaks = c(2003, 2010, 2015, 2020, 2026),
    labels = c("2004-2010", "2011-2015", "2016-2020", "2021-2026"), include.lowest = TRUE)]
  tto[, FDA_YEAR_BAND := fifelse(is.na(FDA_YEAR_BAND), "Unknown", as.character(FDA_YEAR_BAND))]

  # Recorded serious outcomes. Absence of an OUTC row means no recorded
  # occurrence in this audit; it does not prove the outcome did not occur.
  outcome <- data.table::fread(file.path(baseline_dir, "data", "outcome_records.csv"), encoding = "UTF-8",
    colClasses = list(character = c("PRIMARYID", "OUTC_COD")))
  outcome <- unique(outcome[PRIMARYID %chin% tto$PRIMARYID])
  outcome_flags <- outcome[, .(
    death = any(OUTC_COD == "DE"), life_threatening = any(OUTC_COD == "LT"),
    hospitalization = any(OUTC_COD == "HO"), any_serious = any(OUTC_COD %chin% c("DE", "LT", "HO"))), by = PRIMARYID]
  tto <- merge(tto, outcome_flags, by = "PRIMARYID", all.x = TRUE, sort = FALSE)
  for (nm in c("death", "life_threatening", "hospitalization", "any_serious")) tto[is.na(get(nm)), (nm) := FALSE]

  target_drugs <- c("albiglutide", "dulaglutide", "exenatide", "liraglutide", "lixisenatide", "semaglutide", "tirzepatide")
  matched <- data.table::fread(file.path(baseline_dir, "data", "matched_drug_records.csv"), encoding = "UTF-8",
    select = c("PRIMARYID", "ROLE_COD", "drug", "main"),
    colClasses = list(character = c("PRIMARYID", "ROLE_COD", "drug")))
  matched <- unique(matched[PRIMARYID %chin% tto$PRIMARYID & main == TRUE & ROLE_COD == "PS" & drug %chin% target_drugs,
    .(PRIMARYID, drug)])
  drug_counts <- rbindlist(lapply(c("available", "unavailable"), function(g) {
    ids <- tto[availability_group == g, PRIMARYID]
    z <- matched[PRIMARYID %chin% ids, .(n = uniqueN(PRIMARYID)), by = drug]
    z[, `:=`(availability_group = g, denominator = length(ids), pct = 100 * n / length(ids))]
    z
  }), fill = TRUE)
  drug_counts <- merge(data.table::CJ(drug = target_drugs, availability_group = c("available", "unavailable")),
    drug_counts, by = c("drug", "availability_group"), all.x = TRUE)
  drug_counts[is.na(n), `:=`(n = 0L, denominator = ifelse(availability_group == "available", 65L, 176L), pct = 0)]
  drug_counts[, p_value_overall := {
    membership <- tto[, .(PRIMARYID, availability_group)]
    membership[, membership := PRIMARYID %chin% matched[drug == .BY$drug, PRIMARYID]]
    tab <- table(membership$availability_group, membership$membership)
    tryCatch(stats::fisher.test(tab)$p.value, error = function(e) NA_real_)
  }, by = drug]
  drug_counts[, smd := {
    pa <- n[availability_group == "available"] / denominator[availability_group == "available"]
    pu <- n[availability_group == "unavailable"] / denominator[availability_group == "unavailable"]
    p <- (pa + pu) / 2
    if (p > 0 && p < 1) (pa - pu) / sqrt(p * (1 - p)) else if (pa == pu) 0 else NA_real_
  }, by = drug]
  emit(drug_counts, "drug_distribution")

  # A descriptive comparison table for categorical variables. The p-value is
  # the omnibus association test and is repeated on each level row for export.
  smd_binary <- function(p1, p0) {
    p <- (p1 + p0) / 2
    if (!is.finite(p) || p <= 0 || p >= 1) return(ifelse(p1 == p0, 0, NA_real_))
    (p1 - p0) / sqrt(p * (1 - p))
  }
  categorical_compare <- function(variable, label) {
    x <- copy(tto)
    x[, level := as.character(get(variable))]
    x[is.na(level) | !nzchar(level), level := "Unknown"]
    tab <- table(x$availability_group, x$level)
    p <- tryCatch({
      if (nrow(tab) < 2L || ncol(tab) < 2L) NA_real_ else stats::fisher.test(tab, workspace = 2e7)$p.value
    }, error = function(e) tryCatch(stats::chisq.test(tab, correct = FALSE)$p.value, error = function(e) NA_real_))
    out <- x[, .(n = .N), by = .(availability_group, level)]
    den <- x[, .(denominator = .N), by = availability_group]
    out <- merge(out, den, by = "availability_group")
    out[, pct := 100 * n / denominator]
    wide <- dcast(out, level ~ availability_group, value.var = c("n", "pct"), fill = 0)
    rename_map <- c(n_available = "available_n", pct_available = "available_pct",
      n_unavailable = "unavailable_n", pct_unavailable = "unavailable_pct")
    for (old in names(rename_map)) if (old %in% names(wide)) data.table::setnames(wide, old, rename_map[[old]])
    for (g in c("available", "unavailable")) for (m in c("n", "pct")) {
      nm <- paste(g, m, sep = "_"); if (!nm %in% names(wide)) wide[, (nm) := 0]
    }
    wide[, `:=`(variable = label, p_value_overall = p,
      smd = smd_binary(available_pct / 100, unavailable_pct / 100))]
    setcolorder(wide, c("variable", "level", "available_n", "available_pct", "unavailable_n", "unavailable_pct", "p_value_overall", "smd"))
    wide[]
  }
  categorical <- rbindlist(list(
    categorical_compare("AGE_GROUP", "Age group"),
    categorical_compare("SEX_GROUP", "Sex"),
    categorical_compare("FDA_YEAR_BAND", "FDA year band"),
    categorical_compare("death", "Recorded death"),
    categorical_compare("life_threatening", "Recorded life-threatening outcome"),
    categorical_compare("hospitalization", "Recorded hospitalization"),
    categorical_compare("any_serious", "Any recorded serious outcome"),
    categorical_compare("event_date_present", "Event date present"),
    categorical_compare("fda_date_present", "FDA date present"),
    categorical_compare("any_start_date", "Any valid start date present"),
    categorical_compare("all_drugs_start_date", "All drug start dates present"),
    categorical_compare("all_required_date_fields_present", "All required dates present"),
    categorical_compare("date_order_valid", "Date order valid"),
    categorical_compare("date_information_complete_but_rejected", "Dates complete but TTO rejected")
  ), fill = TRUE)
  emit(categorical, "categorical_comparison")
  emit(categorical[variable %chin% c("Event date present", "FDA date present", "Any valid start date present",
    "All drug start dates present", "All required dates present", "Date order valid",
    "Dates complete but TTO rejected")], "date_completeness_comparison")

  numeric_compare <- function(variable, label) {
    z <- tto[is.finite(get(variable)), .(value = as.numeric(get(variable)), availability_group)]
    p <- if (length(unique(z$availability_group)) == 2L) tryCatch(stats::wilcox.test(value ~ availability_group, data = z, exact = FALSE)$p.value, error = function(e) NA_real_) else NA_real_
    smd <- if (all(c("available", "unavailable") %chin% z$availability_group)) {
      a <- z[availability_group == "available", value]; u <- z[availability_group == "unavailable", value]
      pooled <- sqrt(((length(a) - 1) * var(a) + (length(u) - 1) * var(u)) / max(1, length(a) + length(u) - 2))
      if (is.finite(pooled) && pooled > 0) (mean(a) - mean(u)) / pooled else NA_real_
    } else NA_real_
    z[, .(variable = label, available_n = sum(availability_group == "available"), available_median = median(value[availability_group == "available"]),
      available_q1 = quantile(value[availability_group == "available"], .25), available_q3 = quantile(value[availability_group == "available"], .75),
      unavailable_n = sum(availability_group == "unavailable"), unavailable_median = median(value[availability_group == "unavailable"]),
      unavailable_q1 = quantile(value[availability_group == "unavailable"], .25), unavailable_q3 = quantile(value[availability_group == "unavailable"], .75),
      p_value = p, smd = smd)]
  }
  numeric <- rbindlist(list(numeric_compare("AGE_YEARS", "Age in years"), numeric_compare("FDA_YEAR", "FDA receipt year")), fill = TRUE)
  emit(numeric, "numeric_comparison")

  # Full annual distribution is kept separately from the grouped-year test.
  annual <- tto[, .(n = .N), by = .(availability_group, FDA_YEAR)][order(FDA_YEAR, availability_group)]
  annual <- merge(annual, tto[, .(denominator = .N), by = availability_group], by = "availability_group")
  annual[, pct := 100 * n / denominator]
  emit(annual, "annual_distribution")

  unavailable_reason <- tto[TTO_AVAILABLE == FALSE, .(n = .N), by = status][order(-n)]
  unavailable_reason[, pct_of_unavailable := 100 * n / sum(n)]
  emit(unavailable_reason, "unavailable_reasons")

  valid <- tto[TTO_AVAILABLE == TRUE]
  zero_day <- rbindlist(lapply(c("study_drug_group", target_drugs), function(g) {
    z <- tto_all[drug == g & status == "valid" & is.finite(TTO_DAYS)]
    data.table(group = g, valid_n = nrow(z), zero_day_n = sum(z$TTO_DAYS == 0), zero_day_pct = if (nrow(z)) 100 * mean(z$TTO_DAYS == 0) else NA_real_, positive_n = sum(z$TTO_DAYS > 0))
  }), fill = TRUE)
  emit(zero_day, "zero_day_description")

  positive_fit <- rbindlist(lapply(c("study_drug_group", target_drugs), function(g) {
    days <- tto_all[drug == g & status == "valid" & is.finite(TTO_DAYS) & TTO_DAYS > 0, TTO_DAYS]
    fit <- weibull_fit(days, minimum_n = 10L)$summary
    cbind(data.table(group = g, zero_day_excluded = TRUE), summarize_numeric(days), fit)
  }), fill = TRUE)
  emit(positive_fit, "positive_tto_weibull_fit")

  emit(tto[, .(PRIMARYID, CASEID, drug = "study_drug_group", availability_group, status, TTO_DAYS,
    EVENT_DATE, FDA_DATE, START_DATE, AGE_YEARS, AGE_GROUP, SEX_GROUP, REPORTER, FDA_YEAR,
    death, life_threatening, hospitalization, any_serious, event_date_present, fda_date_present,
    any_start_date, all_drugs_start_date, all_required_date_fields_present, date_order_valid)], "report_level_audit", "audit")
  emit(individual, "individual_tto_rows", "audit")
  emit(matched, "target_drug_membership", "audit")
  emit(data.table(metric = c("tto_reports", "available_reports", "unavailable_reports", "available_fraction_pct", "valid_zero_day_reports", "valid_positive_reports"),
    value = c(nrow(tto), sum(tto$TTO_AVAILABLE), sum(!tto$TTO_AVAILABLE), 100 * mean(tto$TTO_AVAILABLE), sum(valid$TTO_DAYS == 0), sum(valid$TTO_DAYS > 0))), "population_summary")
  file.copy(c("13_tto_availability_comparison.R", helper), file.path(out, "code"), overwrite = TRUE)

  age_unknown_a <- categorical[variable == "Age group" & level == "Unknown", available_pct]
  age_unknown_u <- categorical[variable == "Age group" & level == "Unknown", unavailable_pct]
  sex_unknown_a <- categorical[variable == "Sex" & level == "Unknown", available_pct]
  sex_unknown_u <- categorical[variable == "Sex" & level == "Unknown", unavailable_pct]
  serious_a <- categorical[variable == "Any recorded serious outcome" & level == "TRUE", available_pct]
  serious_u <- categorical[variable == "Any recorded serious outcome" & level == "TRUE", unavailable_pct]
  serious_p <- categorical[variable == "Any recorded serious outcome", p_value_overall][1L]
  event_a <- categorical[variable == "Event date present" & level == "TRUE", available_pct]
  event_u <- categorical[variable == "Event date present" & level == "TRUE", unavailable_pct]
  all_date_a <- categorical[variable == "All required dates present" & level == "TRUE", available_pct]
  all_date_u <- categorical[variable == "All required dates present" & level == "TRUE", unavailable_pct]
  complete_rejected <- categorical[variable == "Dates complete but TTO rejected" & level == "TRUE", unavailable_n]
  report <- c(
    "# TTO availability comparison",
    "",
    sprintf("Report-level pooled study_drug_group analysis: %d target reports, %d TTO-available and %d TTO-unavailable.", nrow(tto), sum(tto$TTO_AVAILABLE), sum(!tto$TTO_AVAILABLE)),
    "TTO availability required parseable event, FDA receipt and all target-drug start dates, a non-negative interval, and the prespecified date-order rules.",
    "",
    "Drug membership is tabulated by report. A report with more than one primary-suspect study drug can contribute to more than one drug row, so percentages are within the availability group and do not sum to 100%.",
    sprintf("The largest drug groups were semaglutide (available %d/%d, %.1f%%; unavailable %d/%d, %.1f%%) and tirzepatide (available %d/%d, %.1f%%; unavailable %d/%d, %.1f%%). Drug-level p values and standardized differences are in tables/drug_distribution.csv.",
      drug_counts[drug == "semaglutide" & availability_group == "available", n], 65, drug_counts[drug == "semaglutide" & availability_group == "available", pct],
      drug_counts[drug == "semaglutide" & availability_group == "unavailable", n], 176, drug_counts[drug == "semaglutide" & availability_group == "unavailable", pct],
      drug_counts[drug == "tirzepatide" & availability_group == "available", n], 65, drug_counts[drug == "tirzepatide" & availability_group == "available", pct],
      drug_counts[drug == "tirzepatide" & availability_group == "unavailable", n], 176, drug_counts[drug == "tirzepatide" & availability_group == "unavailable", pct]),
    "",
    sprintf("Known-age distributions were similar (median 60 years in both groups; Wilcoxon p=%.3f), but age was unknown in %.1f%% of available and %.1f%% of unavailable reports. Sex was unknown in %.1f%% and %.1f%%, respectively. The grouped FDA-year distribution showed p=%.3f; numeric receipt year had p=%.3f.",
      numeric[variable == "Age in years", p_value], age_unknown_a, age_unknown_u, sex_unknown_a, sex_unknown_u,
      categorical[variable == "FDA year band", p_value_overall][1L], numeric[variable == "FDA receipt year", p_value][1L]),
    sprintf("Any recorded serious outcome occurred in %.1f%% of available and %.1f%% of unavailable reports (descriptive p=%.3f). Components are reported separately for death, life-threatening outcome and hospitalization.", serious_a, serious_u, serious_p),
    "",
    sprintf("Date completeness was strongly selective: an event date was present in %.1f%% of available versus %.1f%% of unavailable reports; all required dates were present in %.1f%% versus %.1f%%. %d unavailable reports had complete required dates but were rejected by a logical TTO rule or multi-drug consistency rule.", event_a, event_u, all_date_a, all_date_u, complete_rejected),
    "Unavailability reasons are in tables/unavailable_reasons.csv. The leading reasons were missing/invalid event date, no matching therapy for all index rows, and missing/invalid start date.",
    "",
    "Zero-day intervals were retained and described separately in tables/zero_day_description.csv. The pooled analysis had 7 zero-day reports and 58 positive intervals.",
    "The Weibull fit in tables/positive_tto_weibull_fit.csv uses only positive TTO values. For the pooled group, n=58, median=49.5 days (IQR 14.25-137.25), scale=101.6 days (95% CI 67.8-152.4), and shape=0.674 (95% CI 0.560-0.811). This describes selected reported intervals and is not a treated-population risk function.",
    "",
    "P values and standardized differences describe possible selection in date availability; they do not identify the missingness mechanism or establish causality. A missing OUTC record does not prove that a serious outcome did not occur."
  )
  writeLines(report, file.path(out, "TTO_AVAILABILITY_REPORT.md"), useBytes = TRUE)
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  writeLines("COMPLETED", file.path(out, "RUN_STATUS.txt"))
  complete <- TRUE
  cat("COMPLETED:", out, "\n")
  invisible(out)
}

if (sys.nframe() == 0L) {
  baseline <- Sys.getenv("BASELINE_DIR")
  if (!nzchar(baseline)) stop("Set BASELINE_DIR to a completed new_reanalysis run.")
  run_tto_availability_comparison(baseline_dir = baseline,
    output_root = Sys.getenv("OUT_ROOT", "F:/GLP1-RA/new/05_结果/tto_availability_comparison"))
}
