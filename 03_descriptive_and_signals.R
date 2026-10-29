run_descriptive_and_signals <- function(context, cohort) {
  core <- cohort$core
  config <- context$config
  groups <- group_ids(cohort$pairs, config$drugs)
  categorical <- list()
  continuous <- list()
  annual <- list()
  for (group in names(groups)) {
    all_exposed <- core[PRIMARYID %chin% groups[[group]]]
    for (population in c("all_exposed", "target_PT")) {
      selected <- if (population == "target_PT") all_exposed[TARGET == TRUE] else all_exposed
      for (variable in c("SEX_GROUP", "AGE_GROUP", "REPORTER", "COUNTRY")) {
        distribution <- selected[, .(n = .N), by = .(level = get(variable))]
        distribution[, `:=`(drug = group, population = population, variable = variable,
          denominator = nrow(selected), percent = if (nrow(selected)) 100 * n / nrow(selected) else NA_real_)]
        categorical[[length(categorical) + 1L]] <- distribution
      }
      age_summary <- summarize_numeric(selected$AGE_YEARS)
      age_summary[, `:=`(drug = group, population = population, variable = "age_years",
        denominator = nrow(selected), n_missing = nrow(selected) - n_observed)]
      continuous[[length(continuous) + 1L]] <- age_summary
    }
    calendar <- data.table::data.table(year = seq(config$first_quarter %/% 4L, config$last_quarter %/% 4L))
    yearly <- all_exposed[, .(exposed_reports = .N, target_reports = sum(TARGET)), by = .(year = SOURCE_YEAR)]
    yearly <- merge(calendar, yearly, by = "year", all.x = TRUE)
    yearly[is.na(exposed_reports), `:=`(exposed_reports = 0L, target_reports = 0L)]
    coverage <- context$tables$quarter_coverage
    yearly[, supplied_quarters := vapply(year, function(selected_year) {
      sum(coverage$present & substr(coverage$quarter, 1L, 4L) == selected_year)
    }, integer(1))]
    yearly[, `:=`(drug = group, full_year_available = supplied_quarters == 4L,
      target_percent_of_reports = ifelse(exposed_reports > 0, 100 * target_reports / exposed_reports, NA_real_))]
    annual[[group]] <- yearly
  }
  write_result(context, data.table::rbindlist(categorical), "descriptive_categories")
  write_result(context, data.table::rbindlist(continuous), "descriptive_age")
  write_result(context, data.table::rbindlist(annual), "annual_report_counts")

  outcomes <- load_faers_table(context, "OUTC", c("PRIMARYID", "OUTC_COD"))
  outcomes <- unique(outcomes[PRIMARYID %chin% unique(cohort$pairs$PRIMARYID),
    .(PRIMARYID, OUTC_COD = toupper(clean_text(OUTC_COD)))])
  outcomes <- outcomes[nzchar(OUTC_COD)]
  outcome_rows <- list()
  serious_codes <- c("DE", "LT", "HO", "DS", "CA", "RI", "OT")
  for (group in names(groups)) {
    for (population in c("all_exposed", "target_PT")) {
      ids <- groups[[group]]
      if (population == "target_PT") ids <- intersect(ids, core[TARGET == TRUE, PRIMARYID])
      observed <- outcomes[PRIMARYID %chin% ids]
      for (category in c(serious_codes, "ANY_RECORDED_SERIOUS", "NO_RECORDED_OUTCOME")) {
        number <- if (category == "ANY_RECORDED_SERIOUS") {
          data.table::uniqueN(observed[OUTC_COD %chin% serious_codes, PRIMARYID])
        } else if (category == "NO_RECORDED_OUTCOME") {
          length(setdiff(ids, observed$PRIMARYID))
        } else data.table::uniqueN(observed[OUTC_COD == category, PRIMARYID])
        outcome_rows[[length(outcome_rows) + 1L]] <- data.table::data.table(drug = group,
          population = population, outcome = category, n = number, denominator = length(ids),
          percent = if (length(ids)) 100 * number / length(ids) else NA_real_)
      }
    }
  }
  write_result(context, outcomes, "outcome_records", detail = TRUE)
  write_result(context, data.table::rbindlist(outcome_rows), "descriptive_outcomes")
  rm(outcomes)
  gc(verbose = FALSE)

  indications <- load_faers_table(context, "INDI", c("PRIMARYID", "DRUG_SEQ", "INDI_PT"))
  indications <- unique(indications[PRIMARYID %chin% unique(cohort$pairs$PRIMARYID),
    .(PRIMARYID, DRUG_SEQ = clean_text(DRUG_SEQ), INDI_PT = clean_text(INDI_PT))])
  exposure_keys <- unique(cohort$matched[main & ROLE_COD == "PS" & nzchar(DRUG_SEQ), .(PRIMARYID, DRUG_SEQ, drug)])
  indications <- merge(exposure_keys, indications[nzchar(DRUG_SEQ) & nzchar(INDI_PT)],
    by = c("PRIMARYID", "DRUG_SEQ"), allow.cartesian = TRUE)
  indications[, TARGET := PRIMARYID %chin% core[TARGET == TRUE, PRIMARYID]]
  write_result(context, indications[, .(exposed_reports = data.table::uniqueN(PRIMARYID),
    target_reports = data.table::uniqueN(PRIMARYID[TARGET])), by = .(drug, INDI_PT)], "indication_counts")
  rm(indications)
  gc(verbose = FALSE)

  main <- signals_for_groups(core, cohort$pairs, config$drugs, "main_PS")
  main[, fisher_q_BH := stats::p.adjust(fisher_p, method = "BH")]
  write_result(context, main, "signals_main")
  sensitivity <- list()
  for (definition in c("expanded", "restricted")) {
    sensitivity[[definition]] <- signals_for_groups(core, exposure_pairs(cohort$matched, definition),
      config$drugs, paste0("names_", definition, "_PS"))
  }
  sensitivity$PS_SS <- signals_for_groups(core, exposure_pairs(cohort$matched, roles = c("PS", "SS")),
    config$drugs, "main_names_PS_SS")
  sensitivity$all_roles <- signals_for_groups(core, exposure_pairs(cohort$matched, roles = c("PS", "SS", "C", "I")),
    config$drugs, "main_names_all_known_roles")
  full_years <- context$tables$annual_report_counts[drug == "study_drug_group" & full_year_available, year]
  sensitivity$full_years <- signals_for_groups(core[SOURCE_YEAR %in% full_years], cohort$pairs,
    config$drugs, "complete_source_years_only")
  multi_ids <- cohort$pairs[, .N, by = PRIMARYID][N > 1L, PRIMARYID]
  sensitivity$single_study_drug <- signals_for_groups(core[!PRIMARYID %chin% multi_ids], cohort$pairs,
    config$drugs, "exclude_multiple_study_drugs_from_entire_universe")
  sensitivity <- data.table::rbindlist(sensitivity)
  sensitivity[, fisher_q_BH := stats::p.adjust(fisher_p, method = "BH"), by = analysis]
  write_result(context, sensitivity, "signals_sensitivity")
  subgroup <- list()
  for (variable in c("SEX_GROUP", "AGE_GROUP")) {
    for (level in sort(unique(core[[variable]]))) {
      subgroup[[length(subgroup) + 1L]] <- signals_for_groups(core[get(variable) == level], cohort$pairs,
        config$drugs, paste0("subgroup_", variable), level)
    }
  }
  subgroup <- data.table::rbindlist(subgroup)
  subgroup[, fisher_q_BH := stats::p.adjust(fisher_p, method = "BH"), by = analysis]
  write_result(context, subgroup, "signals_subgroups")
  invisible(main)
}

run_comedication <- function(context, cohort) {
  groups <- group_ids(cohort$pairs, context$config$drugs)
  target_ids <- cohort$core[TARGET == TRUE, PRIMARYID]
  frequencies <- list()
  associations <- list()
  coverage <- list()
  for (group in names(groups)) {
    ids <- groups[[group]]
    target_total <- sum(ids %chin% target_ids)
    non_target_total <- length(ids) - target_total
    records <- cohort$comed[PRIMARYID %chin% ids]
    records[, TARGET := PRIMARYID %chin% target_ids]
    summary <- records[, .(co_target_reports = data.table::uniqueN(PRIMARYID[TARGET]),
      co_other_reports = data.table::uniqueN(PRIMARYID[!TARGET])), by = .(comed_name, name_source)]
    summary[, `:=`(drug = group, target_denominator = target_total, other_denominator = non_target_total,
      target_percent = if (target_total) 100 * co_target_reports / target_total else NA_real_,
      other_percent = if (non_target_total) 100 * co_other_reports / non_target_total else NA_real_)]
    frequencies[[group]] <- summary
    coverage[[group]] <- data.table::data.table(drug = group, exposed_reports = length(ids),
      target_reports = target_total, target_with_recorded_C_drug = data.table::uniqueN(records[TARGET == TRUE, PRIMARYID]),
      all_with_recorded_C_drug = data.table::uniqueN(records$PRIMARYID))
    tested <- summary[co_target_reports + co_other_reports >= context$config$minimum_comed_reports]
    for (index in seq_len(nrow(tested))) {
      row <- tested[index]
      result <- signal_metrics(row$co_target_reports, row$co_other_reports,
        target_total - row$co_target_reports, non_target_total - row$co_other_reports)
      result[, `:=`(drug = group, comed_name = row$comed_name, name_source = row$name_source)]
      associations[[length(associations) + 1L]] <- result
    }
  }
  write_result(context, data.table::rbindlist(frequencies), "comedication_frequencies")
  write_result(context, data.table::rbindlist(coverage), "comedication_coverage")
  if (length(associations)) {
    associations <- data.table::rbindlist(associations)
    associations[, fisher_q_BH := stats::p.adjust(fisher_p, method = "BH"), by = drug]
  } else {
    associations <- signal_metrics(0, 0, 0, 0)[0]
    associations[, `:=`(drug = character(), comed_name = character(), name_source = character(), fisher_q_BH = numeric())]
  }
  write_result(context, associations, "comedication_associations")
  invisible(associations)
}
