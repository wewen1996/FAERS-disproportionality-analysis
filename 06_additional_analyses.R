make_additional_catalog <- function(context, definitions) {
  raw <- load_faers_table(context, "DRUG", c("PRIMARYID", "DRUGNAME", "ROLE_COD", "DRUG_SEQ"), "PROD_AI")
  raw[, `:=`(NAME_ID = .I, DRUGNAME = clean_text(DRUGNAME), PROD_AI = clean_text(PROD_AI),
    ROLE_COD = toupper(clean_text(ROLE_COD)), DRUG_SEQ = clean_text(DRUG_SEQ))]
  catalog <- raw[, .(NAME_ID, DRUGNAME, PROD_AI)]
  matched <- match_additional_drug_names(catalog, definitions)
  hits <- merge(matched$hits, raw[, .(NAME_ID, PRIMARYID, DRUGNAME, PROD_AI, DRUG_SEQ, ROLE_COD)],
    by = "NAME_ID", all.x = TRUE)
  hits <- unique(hits[ROLE_COD %chin% c("PS", "SS", "C", "I")])
  write_result(context, hits, "additional_drug_name_matches", detail = TRUE)
  write_result(context, matched$catalog[, .(NAME_ID, DRUGNAME, PROD_AI, combination_flag)],
    "additional_drug_name_catalog", detail = TRUE)
  rm(raw, catalog, matched)
  gc(verbose = FALSE)
  hits
}

make_additional_indication_flags <- function(context, core, definitions) {
  indi <- load_faers_table(context, "INDI", c("PRIMARYID", "INDI_PT"), "DRUG_SEQ")
  indi <- indi[PRIMARYID %chin% core$PRIMARYID]
  indi[, INDI_PT := toupper(clean_text(INDI_PT))]
  flag_names <- names(definitions$conditions)
  flags <- data.table::data.table(PRIMARYID = core$PRIMARYID)
  for (flag_name in flag_names) {
    terms <- definitions$conditions[[flag_name]]
    values <- indi[, .(value = any(INDI_PT %chin% terms)), by = PRIMARYID]
    data.table::setnames(values, "value", flag_name)
    flags <- merge(flags, values, by = "PRIMARYID", all.x = TRUE, sort = FALSE)
  }
  for (flag_name in flag_names) flags[is.na(get(flag_name)), (flag_name) := FALSE]
  write_result(context, indi, "additional_indication_records", detail = TRUE)
  write_result(context, flags, "additional_indication_flags", detail = TRUE)
  flags
}

active_comparator_one <- function(core, indication_flags, target_ids, comparator_ids,
                                  target_drug, comparator, comparator_ingredients, stratum) {
  eligible_ids <- core$PRIMARYID
  if (stratum != "all_reports") eligible_ids <- indication_flags[get(stratum) == TRUE, PRIMARYID]
  eligible <- core[PRIMARYID %chin% eligible_ids]
  target_ids <- intersect(target_ids, eligible$PRIMARYID)
  comparator_ids <- intersect(comparator_ids, eligible$PRIMARYID)
  coexposed <- intersect(target_ids, comparator_ids)
  target_only <- setdiff(target_ids, coexposed)
  comparator_only <- setdiff(comparator_ids, coexposed)
  event_ids <- eligible[TARGET == TRUE, PRIMARYID]
  result <- signal_metrics(sum(target_only %chin% event_ids),
    length(target_only) - sum(target_only %chin% event_ids),
    sum(comparator_only %chin% event_ids),
    length(comparator_only) - sum(comparator_only %chin% event_ids))
  result[, `:=`(analysis = "active_comparator", target_drug = target_drug,
    comparator = comparator, comparator_ingredients = paste(comparator_ingredients, collapse = ";"),
    stratum = stratum, target_reports = length(target_only), comparator_reports = length(comparator_only),
    coexposed_reports = length(coexposed), eligible_background_reports = nrow(eligible))]
  result
}

run_active_comparator <- function(context, cohort, definitions, additional_hits, indication_flags) {
  core <- cohort$core
  targets <- intersect(context$config$active_comparator_targets, context$config$drugs)
  strata <- c("all_reports", "diabetes_recorded", "weight_management_recorded", "metabolic_or_weight_recorded")
  output <- list()
  audit <- list()
  for (target_drug in targets) {
    target_ids <- unique(cohort$matched[main & ROLE_COD == "PS" & drug == target_drug, PRIMARYID])
    for (comparator in names(definitions$comparators)) {
      ingredients <- definitions$comparators[[comparator]]
      comparator_ids <- unique(additional_hits[ingredient %chin% ingredients & ROLE_COD == "PS" &
        !combination_flag, PRIMARYID])
      for (stratum in strata) {
        output[[length(output) + 1L]] <- active_comparator_one(core, indication_flags,
          target_ids, comparator_ids, target_drug, comparator, ingredients, stratum)
      }
      audit[[length(audit) + 1L]] <- data.table::data.table(target_drug = target_drug,
        comparator = comparator, comparator_ingredients = paste(ingredients, collapse = ";"),
        target_raw_reports = data.table::uniqueN(target_ids), comparator_raw_reports = data.table::uniqueN(comparator_ids),
        coexposed_raw_reports = data.table::uniqueN(intersect(target_ids, comparator_ids)))
    }
  }
  result <- data.table::rbindlist(output, fill = TRUE)
  result[, fisher_q_BH := stats::p.adjust(fisher_p, method = "BH"), by = .(target_drug, stratum)]
  write_result(context, result, "active_comparator_signals")
  write_result(context, data.table::rbindlist(audit), "active_comparator_exposure_audit")
  result
}

temporal_signal_one <- function(core, pairs, group, analysis, cutoff, cutoff_label, period) {
  selected <- if (period == "before") core[!is.na(FDA_DATE) & FDA_DATE < cutoff] else core[!is.na(FDA_DATE) & FDA_DATE >= cutoff]
  # The pooled group is the union of study-drug reports, not a drug name in pairs.
  # Count each report once, including reports with multiple study drugs.
  exposed_ids <- if (identical(group, "study_drug_group")) unique(pairs$PRIMARYID) else unique(pairs[drug == group, PRIMARYID])
  ids <- intersect(exposed_ids, selected$PRIMARYID)
  event_ids <- selected[TARGET == TRUE, PRIMARYID]
  a <- sum(ids %chin% event_ids)
  b <- length(ids) - a
  c <- length(setdiff(event_ids, ids))
  d <- nrow(selected) - a - b - c
  result <- signal_metrics(a, b, c, d)
  result[, `:=`(analysis = analysis, drug = group, cutoff = cutoff, cutoff_label = cutoff_label,
    period = period, background_reports = nrow(selected), missing_fda_date_reports = sum(is.na(core$FDA_DATE)))]
  result
}

run_temporal_sensitivity <- function(context, cohort, definitions, indication_flags, additional_hits) {
  core <- cohort$core
  groups <- group_ids(cohort$pairs, context$config$drugs)
  nodes <- definitions$nodes[date %in% context$config$notoriety_cutoffs]
  if (!nrow(nodes)) {
    nodes <- data.table::data.table(date = context$config$notoriety_cutoffs,
      node_id = paste0("custom_", context$config$notoriety_cutoffs), source_id = "USER_CONFIG",
      node_type = "user_configured", drug = NA_character_)
  }
  output <- list()
  counts <- list()
  for (cutoff in context$config$notoriety_cutoffs) {
    node <- nodes[date == cutoff][1L]
    cutoff_label <- if (nrow(node)) node$node_id else as.character(cutoff)
    for (period in c("before", "after")) {
      for (group in names(groups)) output[[length(output) + 1L]] <- temporal_signal_one(
        core, cohort$pairs, group, "temporal_notoriety_vs_all_background", cutoff, cutoff_label, period)
      selected <- if (period == "before") core[!is.na(FDA_DATE) & FDA_DATE < cutoff] else core[!is.na(FDA_DATE) & FDA_DATE >= cutoff]
      for (group in names(groups)) {
        ids <- intersect(groups[[group]], selected$PRIMARYID)
        counts[[length(counts) + 1L]] <- data.table::data.table(cutoff = cutoff, cutoff_label = cutoff_label,
          period = period, drug = group, reports = nrow(selected), exposed_reports = length(ids),
          target_reports = sum(selected[TARGET == TRUE, PRIMARYID] %chin% ids),
          missing_fda_date_reports = sum(is.na(core$FDA_DATE)))
      }
    }
  }
  signals <- data.table::rbindlist(output)
  signals[, fisher_q_BH := stats::p.adjust(fisher_p, method = "BH"), by = .(cutoff, period)]
  write_result(context, signals, "temporal_notoriety_signals")
  write_result(context, data.table::rbindlist(counts), "temporal_notoriety_report_counts")
  active_output <- list()
  for (cutoff in context$config$notoriety_cutoffs) {
    node <- nodes[date == cutoff][1L]
    cutoff_label <- if (nrow(node)) node$node_id else as.character(cutoff)
    for (period in c("before", "after")) {
      selected <- if (period == "before") core[!is.na(FDA_DATE) & FDA_DATE < cutoff] else core[!is.na(FDA_DATE) & FDA_DATE >= cutoff]
      for (target_drug in intersect(context$config$active_comparator_targets, context$config$drugs)) {
        target_ids <- unique(cohort$matched[main & ROLE_COD == "PS" & drug == target_drug, PRIMARYID])
        for (comparator in names(definitions$comparators)) {
          ingredients <- definitions$comparators[[comparator]]
          comparator_ids <- unique(additional_hits[ingredient %chin% ingredients & ROLE_COD == "PS" & !combination_flag, PRIMARYID])
          row <- active_comparator_one(selected, indication_flags, target_ids, comparator_ids,
            target_drug, comparator, ingredients, "all_reports")
          row[, `:=`(cutoff = cutoff, cutoff_label = cutoff_label, period = period)]
          active_output[[length(active_output) + 1L]] <- row
        }
      }
    }
  }
  temporal_active <- data.table::rbindlist(active_output, fill = TRUE)
  temporal_active[, fisher_q_BH := stats::p.adjust(fisher_p, method = "BH"), by = .(cutoff, period, target_drug)]
  write_result(context, temporal_active, "temporal_active_comparator_signals")
  write_result(context, definitions$nodes, "temporal_nodes")
  invisible(signals)
}

run_target_comorbidity <- function(context, cohort, definitions, additional_hits, indication_flags) {
  core <- cohort$core
  groups <- group_ids(cohort$pairs, context$config$drugs)
  target_ids <- core[TARGET == TRUE, PRIMARYID]
  output <- list()
  denominators <- list()
  condition_names <- names(definitions$conditions)
  for (group in names(groups)) {
    ids <- intersect(groups[[group]], target_ids)
    denominators[[length(denominators) + 1L]] <- data.table::data.table(drug = group, target_case_reports = length(ids))
    for (condition_name in condition_names) {
      n <- data.table::uniqueN(indication_flags[PRIMARYID %chin% ids & get(condition_name) == TRUE, PRIMARYID])
      interval <- proportion_recorded(n, length(ids))
      output[[length(output) + 1L]] <- data.table::data.table(drug = group, source = "INDI",
        role_definition = "INDI_PT exact term set", category = condition_name, subcategory = NA_character_,
        n_recorded = n, denominator = length(ids), percent = interval[["percent"]],
        percent_low = interval[["low"]], percent_high = interval[["high"]])
    }
    for (role_definition in c("C_only", "non_primary_known_roles")) {
      roles <- if (role_definition == "C_only") "C" else c("SS", "C", "I")
      for (drug_category in c("diuretic", "anticoagulant", "glucose_lowering")) {
        matched_ids <- unique(additional_hits[PRIMARYID %chin% ids & ROLE_COD %chin% roles &
          category == drug_category & !combination_flag, PRIMARYID])
        interval <- proportion_recorded(length(matched_ids), length(ids))
        output[[length(output) + 1L]] <- data.table::data.table(drug = group,
          source = "DRUG", role_definition = role_definition, category = drug_category, subcategory = NA_character_,
          n_recorded = length(matched_ids), denominator = length(ids), percent = interval[["percent"]],
          percent_low = interval[["low"]], percent_high = interval[["high"]])
      }
      for (ingredient_name in unique(additional_hits$ingredient)) {
        matched_ids <- unique(additional_hits[PRIMARYID %chin% ids & ROLE_COD %chin% roles &
          ingredient == ingredient_name & !combination_flag, PRIMARYID])
        if (!length(matched_ids)) next
        interval <- proportion_recorded(length(matched_ids), length(ids))
        output[[length(output) + 1L]] <- data.table::data.table(drug = group,
          source = "DRUG", role_definition = role_definition, category = "specific_ingredient",
          subcategory = ingredient_name, n_recorded = length(matched_ids), denominator = length(ids),
          percent = interval[["percent"]], percent_low = interval[["low"]], percent_high = interval[["high"]])
      }
    }
  }
  result <- data.table::rbindlist(output, fill = TRUE)
  write_result(context, result, "target_case_comorbidity_and_comedication")
  write_result(context, data.table::rbindlist(denominators), "target_case_denominators")
  result
}

hill_strength <- function(rows) {
  rows <- rows[status == "ok" & is.finite(ROR_low)]
  if (!nrow(rows)) return(list(assessment = "Not estimable", score = NA_real_, evidence = "No estimable active-comparator ROR"))
  positive <- rows[ROR_low > 1 & a >= 3]
  if (nrow(positive) == nrow(rows)) return(list(assessment = "Supports association", score = 2, evidence = paste(nrow(positive), "of", nrow(rows), "comparisons positive")))
  if (nrow(positive)) return(list(assessment = "Mixed", score = 1, evidence = paste(nrow(positive), "of", nrow(rows), "comparisons positive")))
  list(assessment = "Does not support", score = 0, evidence = paste("0 of", nrow(rows), "comparisons positive"))
}

build_bradford_hill <- function(context, cohort, definitions) {
  main <- context$tables$signals_main
  active <- context$tables$active_comparator_signals
  temporal <- context$tables$temporal_notoriety_signals
  tto <- context$tables$tto_summary_weibull
  targets <- intersect(context$config$active_comparator_targets, context$config$drugs)
  rows <- list()
  for (target in targets) {
    active_all <- active[target_drug == target & stratum == "all_reports" &
      comparator %chin% names(definitions$comparators)]
    strength <- hill_strength(active_all)
    other <- main[drug %chin% setdiff(context$config$drugs, target)]
    target_main <- main[drug == target]
    target_positive <- nrow(target_main) && isTRUE(target_main$ROR_signal[1L])
    other_positive <- sum(other$ROR_signal, na.rm = TRUE)
    consistency <- if (!nrow(active_all)) c("Not estimable", NA_real_, "No active-comparator rows") else if (all(active_all$status == "ok") && length(unique(active_all$ROR_signal)) == 1L) c(if (isTRUE(active_all$ROR_signal[1L])) "Consistent positive" else "Consistent null", if (isTRUE(active_all$ROR_signal[1L])) 2 else 0, "All configured active-comparator rows agree") else c("Inconsistent", 1, "Configured active-comparator rows do not agree")
    specificity <- if (target_positive && !other_positive) c("Possible drug-specific pattern", 1, "Target signal present while other selected GLP-1 drugs are not flagged") else if (target_positive && other_positive) c("No specificity", 0, "Other selected GLP-1 drugs also show a signal") else c("Not supported", 0, "Target drug is not flagged in the primary comparison")
    temporal_target <- temporal[drug == target]
    temporal_assessment <- if (nrow(temporal_target) && any(temporal_target$period == "after" & temporal_target$ROR_signal, na.rm = TRUE)) c("Partially met", 1, "Post-cutoff temporal signal is present; this can reflect reporting attention") else c("Not demonstrated", 0, "No post-cutoff temporal signal in the available windows")
    tto_target <- tto[drug == target]
    tto_evidence <- if (nrow(tto_target) && tto_target$positive_TTO > 0) paste0("Positive TTO available (n = ", tto_target$positive_TTO, ")") else "No positive TTO evidence available"
    criteria <- list(
      list(name = "Strength of the association", assessment = strength$assessment, score = strength$score,
        description = paste0("Active-comparator ROR assessment: ", strength$evidence, "."), source = "active_comparator_signals.csv"),
      list(name = "Consistency", assessment = consistency[[1L]], score = as.numeric(consistency[[2L]]),
        description = consistency[[3L]], source = "active_comparator_signals.csv"),
      list(name = "Analogy", assessment = if (other_positive) "Supportive analogy" else "Limited analogy",
        score = if (other_positive) 2 else 0, description = "Compared with other selected GLP-1 drugs in the same report-level background.", source = "signals_main.csv"),
      list(name = "Biological plausibility/empirical evidence", assessment = "Manual evidence review required", score = NA_real_,
        description = "Use the target-case INDI and co-reported drug proportions to structure review of constipation, dehydration, vascular disease and heart failure; FAERS co-reporting is not mechanistic evidence.", source = "target_case_comorbidity_and_comedication.csv plus external literature"),
      list(name = "Exclusion of biases/confounders", assessment = "Not met", score = 0,
        description = "Indication bias, notoriety bias, missing exposure denominators and residual co-reporting bias cannot be excluded by spontaneous-report data.", source = "Study design limitation; active comparator and temporal sensitivity analyses"),
      list(name = "Specificity", assessment = specificity[[1L]], score = as.numeric(specificity[[2L]]),
        description = specificity[[3L]], source = "signals_main.csv"),
      list(name = "Temporal relationship", assessment = temporal_assessment[[1L]], score = as.numeric(temporal_assessment[[2L]]),
        description = paste(temporal_assessment[[3L]], tto_evidence, sep = " "), source = "temporal_notoriety_signals.csv and tto_summary_weibull.csv"),
      list(name = "Reversibility", assessment = "Not applicable", score = NA_real_,
        description = "FAERS does not reliably encode treatment discontinuation and dechallenge for this analysis.", source = "Not available in the supplied tables"),
      list(name = "Coherence", assessment = "Manual evidence review required", score = NA_real_,
        description = "Compare the observed signal with clinical, mechanistic and case-report evidence for ischemic colitis; do not infer coherence from FAERS alone.", source = "External clinical and mechanistic literature"))
    for (criterion in criteria) rows[[length(rows) + 1L]] <- data.table::data.table(drug = target,
      criteria = criterion$name, assessment = criterion$assessment, score_0_to_2 = criterion$score,
      description = criterion$description, source_method = criterion$source,
      data_evidence = if (criterion$name == "Strength of the association") strength$evidence else if (criterion$name == "Temporal relationship") temporal_assessment[[3L]] else NA_character_)
  }
  result <- data.table::rbindlist(rows)
  summary <- result[, .(observed_score = sum(score_0_to_2, na.rm = TRUE),
    criteria_scored = sum(!is.na(score_0_to_2)), maximum_observed_score = 2 * sum(!is.na(score_0_to_2))), by = drug]
  write_result(context, result, "bradford_hill_framework")
  write_result(context, summary, "bradford_hill_score_summary")
  write_result(context, data.table::rbindlist(list(definitions$sources, data.table::data.table(
    source_id = "HILL_SCORE_RULE", title = "Exploratory 0-2 score rule used by this script",
    identifier = "0 = no support; 1 = partial/mixed; 2 = support; NA = not data-derived",
    url = NA_character_, checked_on = format(Sys.Date())))), "bradford_hill_sources")
  result
}

run_additional_analyses <- function(context, cohort) {
  definitions <- additional_definitions()
  write_result(context, definitions$sources, "additional_analysis_sources")
  write_result(context, definitions$nodes, "additional_temporal_nodes")
  additional_hits <- make_additional_catalog(context, definitions)
  indication_flags <- make_additional_indication_flags(context, cohort$core, definitions)
  active <- run_active_comparator(context, cohort, definitions, additional_hits, indication_flags)
  temporal <- run_temporal_sensitivity(context, cohort, definitions, indication_flags, additional_hits)
  comorbidity <- run_target_comorbidity(context, cohort, definitions, additional_hits, indication_flags)
  hill <- build_bradford_hill(context, cohort, definitions)
  write_result(context, data.table::data.table(
    module = c("active_comparator", "temporal_notoriety", "bradford_hill", "target_comorbidity_comedication"),
    primary_output = c("active_comparator_signals.csv", "temporal_notoriety_signals.csv",
      "bradford_hill_framework.csv", "target_case_comorbidity_and_comedication.csv"),
    interpretation = c("Report-level target versus comparator ROR; not incidence or causal effect.",
      "Before/after report-date ROR; temporal change may reflect reporting attention.",
      "Exploratory causal-evidence framework; scored rows are data-derived only.",
      "Proportion of target reports with INDI terms or matched co-reported drugs.")), "additional_analysis_manifest")
  invisible(list(active = active, temporal = temporal, comorbidity = comorbidity, hill = hill))
}
