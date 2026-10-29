build_cohort <- function(context) {
  config <- context$config
  mapping <- read_drug_mapping(context)
  demo <- load_faers_table(context, "rawDEMO", c("PRIMARYID", "CASEID", "FDA_DT", "GETDATAYEAR"),
    c("EVENT_DT", "AGE", "AGE_COD", "SEX", "OCCP_COD", "REPORTER_COUNTRY", "WT", "WT_COD"))
  initial_n <- nrow(demo)
  flow_record(context, "01_supplied_rawDEMO_rows", initial_n, 0, initial_n, "DEMO rows",
    "Supplied cleaned dataset, before this script's deduplication; not FDA download row counts")
  demo[, QUARTER := quarter_index(GETDATAYEAR)]
  write_result(context, demo[, .(source_demo_rows = .N), by = .(GETDATAYEAR, QUARTER)], "source_quarters")
  keep <- !is.na(demo$QUARTER) & demo$QUARTER >= config$first_quarter & demo$QUARTER <= config$last_quarter
  demo <- demo[keep]
  flow_record(context, "02_source_quarter_window", initial_n, initial_n - nrow(demo), nrow(demo),
    "DEMO rows", "Quarter selected by GetDataYear, not by FDA_DT")
  covered <- sort(unique(demo$QUARTER))
  missing_quarters <- setdiff(seq(config$first_quarter, config$last_quarter), covered)
  write_result(context, data.table::data.table(quarter = quarter_label(seq(config$first_quarter, config$last_quarter)),
    present = seq(config$first_quarter, config$last_quarter) %in% covered), "quarter_coverage")
  if (length(missing_quarters)) warning("Missing source quarters: ", paste(quarter_label(missing_quarters), collapse = ", "))
  previous_n <- nrow(demo)
  demo <- demo[nzchar(clean_text(PRIMARYID)) & nzchar(clean_text(CASEID))]
  flow_record(context, "03_valid_identifiers", previous_n, previous_n - nrow(demo), nrow(demo),
    "DEMO rows", "Remove missing PRIMARYID or CASEID")
  previous_n <- nrow(demo)
  deduplicated <- deduplicate_demo(demo)
  core <- deduplicated$retained
  write_result(context, deduplicated$primaryid_case_conflicts, "dedup_primaryid_case_conflicts")
  flow_record(context, "04_unique_report_versions", previous_n, deduplicated$duplicate_rows,
    previous_n - deduplicated$duplicate_rows, "report versions", "One row per PRIMARYID")
  flow_record(context, "05_latest_case_version", previous_n - deduplicated$duplicate_rows,
    deduplicated$older_versions, nrow(core), "case reports",
    "Latest source quarter per PRIMARYID, then latest valid FDA_DT per CASEID; numeric ID ties use digit-aware order and other IDs use character order")
  rm(demo, deduplicated, keep)
  gc(verbose = FALSE)
  write_result(context, core[, .(reports = .N), by = .(fda_date_missing = is.na(FDA_DATE))], "dedup_date_quality")

  reactions <- load_faers_table(context, "REAC", c("PRIMARYID", "PT"))
  reactions <- reactions[PRIMARYID %chin% core$PRIMARYID]
  reactions[, PT := toupper(clean_text(PT))]
  with_reaction <- unique(reactions[nzchar(PT), PRIMARYID])
  target_ids <- unique(reactions[PT == config$target_pt, PRIMARYID])
  write_result(context, data.table::data.table(target_pt = config$target_pt,
    matched_rows = reactions[PT == config$target_pt, .N], unique_report_versions = length(target_ids)), "pt_match_audit")
  previous_n <- nrow(core)
  missing_reaction <- core[!PRIMARYID %chin% with_reaction, .(CASEID, PRIMARYID)]
  write_result(context, missing_reaction, "excluded_no_reaction", detail = TRUE)
  core <- core[PRIMARYID %chin% with_reaction]
  core[, TARGET := PRIMARYID %chin% target_ids]
  flow_record(context, "06_nonempty_reaction", previous_n, previous_n - nrow(core), nrow(core),
    "case reports", "Latest retained PRIMARYID must have at least one nonempty REAC.PT")
  rm(reactions, with_reaction, missing_reaction)
  gc(verbose = FALSE)

  drugs <- load_faers_table(context, "DRUG", c("PRIMARYID", "DRUGNAME", "ROLE_COD", "DRUG_SEQ"), "PROD_AI")
  drugs <- drugs[PRIMARYID %chin% core$PRIMARYID]
  drugs[, `:=`(DRUGNAME = as.character(DRUGNAME), ROLE_COD = toupper(clean_text(ROLE_COD)),
                DRUG_SEQ = clean_text(DRUG_SEQ))]
  with_drug <- unique(drugs[!is.na(DRUGNAME) & nzchar(trimws(DRUGNAME)), PRIMARYID])
  previous_n <- nrow(core)
  write_result(context, core[!PRIMARYID %chin% with_drug, .(CASEID, PRIMARYID)], "excluded_no_drug", detail = TRUE)
  core <- core[PRIMARYID %chin% with_drug]
  flow_record(context, "07_analysis_universe", previous_n, previous_n - nrow(core), nrow(core),
    "case reports", "Latest retained PRIMARYID must have at least one named drug; universe for all 2x2 tables")
  if (!nrow(core)) stop("No eligible report-level analysis universe.")
  drugs <- drugs[PRIMARYID %chin% core$PRIMARYID]
  matched <- drugs[DRUGNAME %chin% mapping$original_drug_name,
                   .(PRIMARYID, DRUG_SEQ, DRUGNAME, ROLE_COD)]
  mapping_index <- match(matched$DRUGNAME, mapping$original_drug_name)
  matched[, `:=`(drug = mapping$standardized_ingredient[mapping_index],
    main = mapping$include_main[mapping_index], expanded = mapping$include_sensitivity_expanded[mapping_index],
    restricted = mapping$include_sensitivity_restricted[mapping_index])]
  matched <- unique(matched)
  write_result(context, matched[, .(report_n = data.table::uniqueN(PRIMARYID), drug_record_n = .N),
    by = .(drug, DRUGNAME, ROLE_COD, main, expanded, restricted)], "drug_name_matches")
  names_in_source <- unique(drugs$DRUGNAME)
  audit_pattern <- paste(c(config$drugs, "OZEMPIC", "WEGOVY", "RYBELSUS", "VICTOZA", "SAXENDA",
    "MOUNJARO", "ZEPBOUND", "TRULICITY", "BYETTA", "BYDUREON", "ADLYXIN", "LYXUMIA", "TANZEUM", "EPERZAN"), collapse = "|")
  candidate_names <- names_in_source[!is.na(names_in_source) & grepl(audit_pattern, names_in_source, ignore.case = TRUE)]
  unmatched_names <- setdiff(candidate_names, mapping$original_drug_name)
  write_result(context, drugs[DRUGNAME %chin% unmatched_names,
    .(report_n = data.table::uniqueN(PRIMARYID)), by = .(DRUGNAME, ROLE_COD)][order(-report_n)], "unmapped_name_candidates")
  if (length(unmatched_names)) warning("Unmapped candidate names were audited, not automatically included.")
  primary_pairs <- exposure_pairs(matched)
  exposed_ids <- unique(primary_pairs$PRIMARYID)
  comed <- drugs[PRIMARYID %chin% exposed_ids & ROLE_COD == "C" & !DRUGNAME %chin% mapping$original_drug_name]
  comed[, `:=`(PROD_AI = toupper(clean_text(PROD_AI)), name_fallback = toupper(clean_text(DRUGNAME)))]
  comed[, `:=`(comed_name = ifelse(nzchar(PROD_AI), PROD_AI, name_fallback),
                name_source = ifelse(nzchar(PROD_AI), "PROD_AI", "DRUGNAME_unstandardized"))]
  comed <- unique(comed[nzchar(comed_name), .(PRIMARYID, comed_name, name_source)])
  rm(drugs, with_drug, names_in_source, mapping_index)
  gc(verbose = FALSE)

  core[, `:=`(AGE_YEARS = suppressWarnings(as.numeric(as.character(AGE))),
    SEX_GROUP = toupper(clean_text(SEX)), COUNTRY = toupper(clean_text(REPORTER_COUNTRY)),
    REPORTER = toupper(clean_text(OCCP_COD)), SOURCE_YEAR = QUARTER %/% 4L)]
  core[, AGE_STATUS := data.table::fcase(is.na(AGE_YEARS), "missing_or_nonnumeric",
    !is.finite(AGE_YEARS) | AGE_YEARS < 0 | AGE_YEARS > config$maximum_age, "outside_allowed_range", default = "valid")]
  core[AGE_STATUS != "valid", AGE_YEARS := NA_real_]
  core[, AGE_GROUP := as.character(cut(AGE_YEARS, breaks = c(-Inf, 18, 45, 65, Inf),
    right = FALSE, labels = c("<18", "18-44", "45-64", "65+")))]
  core[is.na(AGE_GROUP), AGE_GROUP := "Unknown"]
  core[!SEX_GROUP %chin% c("F", "M"), SEX_GROUP := "Unknown"]
  core[!nzchar(COUNTRY), COUNTRY := "Unknown"]
  core[!nzchar(REPORTER), REPORTER := "Unknown"]
  core[, MAIN_EXPOSED := PRIMARYID %chin% exposed_ids]
  write_result(context, core[, .(report_n = .N), by = .(AGE_STATUS, AGE_COD)], "age_quality")
  write_result(context, core[MAIN_EXPOSED == TRUE], "exposed_case_demographics", detail = TRUE)
  write_result(context, matched, "matched_drug_records", detail = TRUE)
  write_result(context, primary_pairs, "primary_exposure_pairs", detail = TRUE)
  write_result(context, comed, "comedication_records", detail = TRUE)
  saveRDS(core, file.path(context$out, "data", "analysis_core.rds"), compress = "gzip")
  groups_any <- group_ids(exposure_pairs(matched, roles = c("PS", "SS", "C", "I")), config$drugs)
  groups_ps <- group_ids(primary_pairs, config$drugs)
  counts <- data.table::rbindlist(lapply(names(groups_ps), function(group) {
    selected <- core[PRIMARYID %chin% groups_ps[[group]]]
    data.table::data.table(drug = group, main_any_role_reports = length(groups_any[[group]]),
      main_PS_reports = nrow(selected), target_PT_reports = sum(selected$TARGET),
      without_target_PT = sum(!selected$TARGET))
  }))
  write_result(context, counts, "flow_drug_counts")
  overlap <- primary_pairs[, .(number_of_study_drugs = .N), by = PRIMARYID]
  write_result(context, overlap[, .(report_n = .N, target_report_n = sum(PRIMARYID %chin% target_ids)),
    by = number_of_study_drugs], "study_drug_overlap")
  stopifnot(sum(core$MAIN_EXPOSED) == counts[drug == "study_drug_group", main_PS_reports],
    sum(core$MAIN_EXPOSED & core$TARGET) == counts[drug == "study_drug_group", target_PT_reports])
  list(core = core, matched = matched, pairs = primary_pairs, comed = comed, counts = counts)
}
