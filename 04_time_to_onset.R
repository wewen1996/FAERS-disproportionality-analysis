assemble_tto_records <- function(index_drugs, therapy, demographics, study_end) {
  index_drugs <- unique(index_drugs[, .(PRIMARYID, DRUG_SEQ, drug)])
  therapy <- unique(therapy[, .(PRIMARYID, DRUG_SEQ, START_DT, END_DT)])
  therapy[, therapy_present := TRUE]
  joined <- merge(index_drugs, therapy[nzchar(DRUG_SEQ)],
    by = c("PRIMARYID", "DRUG_SEQ"), all.x = TRUE, allow.cartesian = TRUE)
  joined[, START_DATE := exact_date(START_DT)]
  summaries <- joined[, .(n_index_rows = .N,
    n_missing_sequence = sum(!nzchar(DRUG_SEQ)),
    n_unmatched_therapy = sum(is.na(therapy_present)),
    n_invalid_start = sum(!is.na(therapy_present) & is.na(START_DATE)),
    n_distinct_start = data.table::uniqueN(START_DATE[!is.na(START_DATE)]),
    START_DATE = if (any(!is.na(START_DATE))) min(START_DATE, na.rm = TRUE) else as.Date(NA)),
    by = .(PRIMARYID, drug)]
  records <- merge(summaries, demographics[, .(PRIMARYID, CASEID, EVENT_DT, FDA_DT, FDA_DATE)], by = "PRIMARYID")
  records[, EVENT_DATE := exact_date(EVENT_DT)]
  records[, TTO_DAYS := as.numeric(EVENT_DATE - START_DATE)]
  records[, status := data.table::fcase(
    is.na(EVENT_DATE), "missing_partial_or_invalid_event_date",
    is.na(FDA_DATE), "missing_partial_or_invalid_fda_date",
    n_missing_sequence > 0, "missing_drug_sequence",
    n_unmatched_therapy > 0, "no_matching_therapy_for_all_index_rows",
    n_invalid_start > 0 | n_distinct_start == 0, "missing_partial_or_invalid_start_date",
    n_distinct_start > 1, "multiple_distinct_start_dates",
    EVENT_DATE > FDA_DATE, "event_after_fda_receipt",
    EVENT_DATE > study_end, "event_after_study_end",
    TTO_DAYS < 0, "event_before_drug_start", default = "valid")]
  list(records = records, joined = joined)
}

run_time_to_onset <- function(context, cohort) {
  targets <- cohort$core[TARGET & MAIN_EXPOSED]
  index_drugs <- cohort$matched[main & ROLE_COD == "PS" & PRIMARYID %chin% targets$PRIMARYID]
  therapy <- load_faers_table(context, "THER", c("PRIMARYID", "DRUG_SEQ", "START_DT"), "END_DT")
  therapy <- therapy[PRIMARYID %chin% targets$PRIMARYID]
  therapy[, DRUG_SEQ := clean_text(DRUG_SEQ)]
  assembled <- assemble_tto_records(index_drugs, therapy, targets, context$config$study_end)
  records <- assembled$records
  write_result(context, assembled$joined, "tto_drug_therapy_links", detail = TRUE)
  rm(therapy, assembled)
  gc(verbose = FALSE)
  pooled <- records[, {
    all_valid <- all(status == "valid")
    unique_time <- data.table::uniqueN(TTO_DAYS)
    pooled_status <- if (.N == 1L) status[1L] else if (!all_valid) {
      "incomplete_dates_across_study_drugs"
    } else if (unique_time > 1L) "different_start_dates_across_study_drugs" else "valid"
    .(CASEID = CASEID[1L], EVENT_DT = EVENT_DT[1L], FDA_DT = FDA_DT[1L],
      EVENT_DATE = EVENT_DATE[1L], FDA_DATE = FDA_DATE[1L],
      START_DATE = if (pooled_status == "valid") START_DATE[1L] else as.Date(NA),
      TTO_DAYS = if (pooled_status == "valid") TTO_DAYS[1L] else NA_real_, status = pooled_status)
  }, by = PRIMARYID]
  pooled[, drug := "study_drug_group"]
  records <- data.table::rbindlist(list(records, pooled), fill = TRUE)
  write_result(context, records, "tto_case_level", detail = TRUE)
  statuses <- c("missing_partial_or_invalid_event_date", "missing_partial_or_invalid_fda_date",
    "missing_drug_sequence", "no_matching_therapy_for_all_index_rows", "missing_partial_or_invalid_start_date",
    "multiple_distinct_start_dates", "event_after_fda_receipt", "event_after_study_end", "event_before_drug_start",
    "incomplete_dates_across_study_drugs", "different_start_dates_across_study_drugs", "valid")
  attrition <- data.table::CJ(drug = c("study_drug_group", context$config$drugs), status = statuses)
  attrition <- merge(attrition, records[, .(n = .N), by = .(drug, status)], by = c("drug", "status"), all.x = TRUE)
  attrition[is.na(n), n := 0L]
  write_result(context, attrition, "tto_attrition")
  summaries <- list()
  models <- list()
  for (group in c("study_drug_group", context$config$drugs)) {
    selected <- records[drug == group]
    valid_days <- selected[status == "valid", TTO_DAYS]
    fit <- weibull_fit(valid_days, context$config$minimum_weibull_n)
    summary <- cbind(data.table::data.table(drug = group, target_reports = nrow(selected),
      excluded_date_reports = sum(selected$status != "valid"), valid_nonnegative_TTO = length(valid_days),
      same_day_reports = sum(valid_days == 0), positive_TTO = sum(valid_days > 0)),
      summarize_numeric(valid_days), fit$summary)
    stopifnot(summary$target_reports == cohort$counts[drug == group, target_PT_reports],
      summary$target_reports == summary$excluded_date_reports + summary$valid_nonnegative_TTO,
      summary$valid_nonnegative_TTO == summary$same_day_reports + summary$positive_TTO)
    summaries[[group]] <- summary
    models[group] <- list(fit$model)
  }
  summaries <- data.table::rbindlist(summaries)
  write_result(context, summaries, "tto_summary_weibull")
  saveRDS(models, file.path(context$out, "data", "weibull_models.rds"))
  list(records = records, summary = summaries)
}
