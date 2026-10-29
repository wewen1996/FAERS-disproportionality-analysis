clean_text <- function(value) {
  value <- trimws(as.character(value))
  value[is.na(value)] <- ""
  value
}

identifier_text <- function(value) {
  # R may display numeric FAERS identifiers such as 5000000 as "5e+06".
  # Format numeric IDs without scientific notation before converting to text.
  if (is.numeric(value)) value <- format(value, scientific = FALSE, trim = TRUE, digits = 22L)
  clean_text(value)
}

exact_date <- function(value) {
  value <- clean_text(value)
  distinct <- unique(value)
  parsed <- as.Date(rep(NA_character_, length(distinct)))
  usable <- grepl("^[0-9]{8}$", distinct)
  parsed[usable] <- suppressWarnings(as.Date(distinct[usable], "%Y%m%d"))
  valid <- !is.na(parsed) & format(parsed, "%Y%m%d") == distinct
  parsed[!valid] <- as.Date(NA)
  parsed[match(value, distinct)]
}

quarter_index <- function(value) {
  code <- suppressWarnings(as.integer(as.character(value)))
  year <- 2000L + code %/% 10L
  quarter <- code %% 10L
  result <- year * 4L + quarter - 1L
  result[is.na(code) | !quarter %in% 1:4 | year < 2004L | year > 2099L] <- NA_integer_
  result
}

quarter_label <- function(index) paste0(index %/% 4L, "Q", index %% 4L + 1L)

write_result <- function(context, table, name, detail = FALSE) {
  stopifnot(is.data.frame(table), !anyDuplicated(names(table)))
  directory <- if (detail) "data" else "tables"
  data.table::fwrite(table, file.path(context$out, directory, paste0(name, ".csv")), bom = TRUE, na = "")
  if (!detail) context$tables[[name]] <- table
  invisible(table)
}

flow_record <- function(context, stage, input, excluded, retained, unit, rule) {
  stopifnot(is.finite(input), is.finite(excluded), is.finite(retained),
            input >= 0, excluded >= 0, retained >= 0, input == excluded + retained)
  context$flow <- data.table::rbindlist(list(context$flow, data.table::data.table(
    stage = stage, input_n = input, excluded_n = excluded, retained_n = retained, unit = unit, rule = rule)))
  write_result(context, context$flow, "flow_attrition")
}

load_faers_table <- function(context, table_name, required, optional = character()) {
  path <- file.path(context$config$data_dir, paste0(table_name, ".RData"))
  message("Reading ", table_name, " at ", format(Sys.time(), "%H:%M:%S"))
  store <- new.env(parent = baseenv())
  loaded <- load(path, envir = store)
  candidates <- loaded[vapply(loaded, function(name) is.data.frame(store[[name]]), logical(1))]
  selected <- if (table_name %in% candidates) table_name else candidates
  if (length(selected) != 1L) stop("Expected one table in: ", path)
  table <- data.table::as.data.table(store[[selected]])
  original_names <- names(table)
  upper_names <- toupper(original_names)
  if (anyDuplicated(upper_names)) stop("Ambiguous column names in ", table_name)
  schema <- data.table::data.table(table = table_name, column = original_names,
    storage_class = vapply(table, function(column) paste(class(column), collapse = "/"), character(1)),
    source_rows = nrow(table))
  context$schema <- data.table::rbindlist(list(context$schema, schema))
  write_result(context, context$schema, "input_schema")
  data.table::setnames(table, upper_names)
  missing <- setdiff(required, names(table))
  if (length(missing)) stop(table_name, " lacks: ", paste(missing, collapse = ", "))
  for (column in setdiff(optional, names(table))) data.table::set(table, j = column, value = NA_character_)
  table[, setdiff(names(table), c(required, optional)) := NULL]
  if ("PRIMARYID" %in% names(table)) table[, PRIMARYID := identifier_text(PRIMARYID)]
  table
}

deduplicate_demo <- function(demo) {
  demo[, `:=`(PRIMARYID = identifier_text(PRIMARYID), CASEID = identifier_text(CASEID))]
  # PRIMARYID is an identifier, not a measurement.  Some cleaned exports contain
  # alphanumeric IDs, so retain the full string and use a deterministic tie-breaker
  # instead of coercing it to numeric (which could also lose precision for long IDs).
  demo[, pid_numeric := grepl("^[0-9]+$", PRIMARYID)]
  if (any(!demo$pid_numeric)) {
    warning("PRIMARYID contains non-numeric identifiers; retaining them as character values and using character order for date ties.")
  }
  duplicate_ids <- unique(demo$PRIMARYID[duplicated(demo$PRIMARYID)])
  conflicting <- demo[PRIMARYID %chin% duplicate_ids,
                      .(case_count = data.table::uniqueN(CASEID),
                        cases = paste(sort(unique(CASEID)), collapse = " | "),
                        source_quarters = paste(sort(unique(QUARTER)), collapse = " | "),
                        rows = .N), by = PRIMARYID][case_count > 1L]
  if (nrow(conflicting)) {
    warning("Found ", nrow(conflicting),
      " PRIMARYID values linked to multiple CASEIDs; retaining the row from the latest source quarter. See dedup_primaryid_case_conflicts.csv.")
  }
  demo[, `:=`(FDA_DATE = exact_date(FDA_DT), pid_digits = fifelse(pid_numeric, nchar(PRIMARYID), NA_integer_))]
  demo[, fda_order := as.numeric(FDA_DATE)]
  demo[is.na(fda_order), fda_order := -Inf]
  demo[, source_row := .I]
  # First resolve repeated PRIMARYID values. Source quarter is the primary
  # tie-breaker because the cleaned files can reuse an ID with another CASEID
  # in a later quarterly extract.
  data.table::setorderv(demo,
    c("PRIMARYID", "QUARTER", "fda_order", "pid_numeric", "pid_digits", "source_row"))
  unique_reports <- demo[!duplicated(PRIMARYID, fromLast = TRUE)]
  # Then resolve old versions within each CASEID.
  data.table::setorderv(unique_reports,
    c("CASEID", "fda_order", "pid_numeric", "pid_digits", "PRIMARYID", "QUARTER", "source_row"))
  retained <- unique_reports[!duplicated(CASEID, fromLast = TRUE)]
  duplicate_rows <- nrow(demo) - nrow(unique_reports)
  older_versions <- nrow(unique_reports) - nrow(retained)
  retained[, c("pid_numeric", "pid_digits", "fda_order", "source_row") := NULL]
  stopifnot(!anyDuplicated(retained$CASEID), !anyDuplicated(retained$PRIMARYID))
  list(retained = retained, duplicate_rows = duplicate_rows, older_versions = older_versions,
    primaryid_case_conflicts = conflicting)
}

read_drug_mapping <- function(context) {
  config <- context$config
  mapping <- data.table::as.data.table(readRDS(file.path(config$mapping_dir, "drug_name_mapping.rds")))
  definitions <- c(main = "include_main", expanded = "include_sensitivity_expanded",
                   restricted = "include_sensitivity_restricted")
  required <- c("original_drug_name", "standardized_ingredient", unname(definitions))
  if (length(setdiff(required, names(mapping)))) stop("Mapping schema is incomplete.")
  if (anyDuplicated(mapping$original_drug_name)) stop("Mapping contains duplicate original names.")
  for (column in unname(definitions)) {
    if (!is.logical(mapping[[column]]) || anyNA(mapping[[column]])) stop("Invalid mapping decision: ", column)
  }
  lists <- readRDS(file.path(config$mapping_dir, "drug_lists_for_analysis.rds"))
  list_names <- c(main = "main", expanded = "sensitivity_expanded", restricted = "sensitivity_restricted")
  for (definition in names(definitions)) {
    for (ingredient in config$drugs) {
      if (!ingredient %in% names(lists[[list_names[[definition]]]])) stop("Drug missing from lists: ", ingredient)
      expected <- mapping[standardized_ingredient == ingredient & get(definitions[[definition]]), original_drug_name]
      if (!setequal(expected, lists[[list_names[[definition]]]][[ingredient]])) {
        stop("Mapping and drug list disagree: ", ingredient, "/", definition, ". Re-export both mapping files.")
      }
    }
  }
  mapping <- mapping[standardized_ingredient %chin% config$drugs]
  write_result(context, mapping, "mapping_snapshot", detail = TRUE)
  mapping
}

signal_metrics <- function(a, b, c, d) {
  cells <- as.numeric(c(a, b, c, d))
  stopifnot(length(cells) == 4L, all(is.finite(cells)), all(cells >= 0), all(cells == floor(cells)))
  total <- sum(cells)
  degenerate <- any(c(a + b, c + d, a + c, b + d) == 0)
  answer <- data.table::data.table(a = a, b = b, c = c, d = d, total_reports = total,
    zero_cell_correction = any(cells == 0) && !degenerate,
    ROR = NA_real_, ROR_low = NA_real_, ROR_high = NA_real_,
    PRR = NA_real_, PRR_low = NA_real_, PRR_high = NA_real_,
    pearson_chisq = NA_real_, fisher_p = NA_real_, ROR_signal = FALSE, PRR_signal = FALSE,
    status = if (degenerate) "not_estimable_zero_margin" else "ok")
  if (degenerate) return(answer)
  adjusted <- if (any(cells == 0)) cells + 0.5 else cells
  names(adjusted) <- c("a", "b", "c", "d")
  ror <- adjusted[["a"]] * adjusted[["d"]] / (adjusted[["b"]] * adjusted[["c"]])
  ror_se <- sqrt(sum(1 / adjusted))
  prr <- (adjusted[["a"]] / sum(adjusted[c("a", "b")])) /
    (adjusted[["c"]] / sum(adjusted[c("c", "d")]))
  prr_se <- sqrt(1 / adjusted[["a"]] - 1 / sum(adjusted[c("a", "b")]) +
                 1 / adjusted[["c"]] - 1 / sum(adjusted[c("c", "d")]))
  chisq <- total * (as.numeric(a) * d - as.numeric(b) * c)^2 /
    (as.numeric(a + b) * (c + d) * (a + c) * (b + d))
  fisher <- stats::fisher.test(matrix(cells, nrow = 2L, byrow = TRUE))$p.value
  z_value <- stats::qnorm(0.975)
  answer[, `:=`(ROR = ror, ROR_low = exp(log(ror) - z_value * ror_se),
    ROR_high = exp(log(ror) + z_value * ror_se), PRR = prr,
    PRR_low = exp(log(prr) - z_value * prr_se), PRR_high = exp(log(prr) + z_value * prr_se),
    pearson_chisq = chisq, fisher_p = fisher)]
  answer[, `:=`(ROR_signal = a >= 3 & ROR_low > 1,
    PRR_signal = a >= 3 & PRR >= 2 & pearson_chisq >= 4)]
  answer
}

exposure_pairs <- function(drug_rows, definition = "main", roles = "PS") {
  unique(drug_rows[get(definition) & ROLE_COD %chin% roles, .(PRIMARYID, drug)])
}

group_ids <- function(pairs, drugs) {
  groups <- lapply(drugs, function(ingredient) unique(pairs[drug == ingredient, PRIMARYID]))
  names(groups) <- drugs
  c(list(study_drug_group = unique(pairs$PRIMARYID)), groups)
}

signals_for_groups <- function(core, pairs, drugs, analysis, stratum = "all") {
  targets <- core[TARGET == TRUE, PRIMARYID]
  groups <- group_ids(pairs, drugs)
  data.table::rbindlist(lapply(names(groups), function(group) {
    exposed <- intersect(groups[[group]], core$PRIMARYID)
    event_count <- sum(exposed %chin% targets)
    result <- signal_metrics(event_count, length(exposed) - event_count,
      length(targets) - event_count, nrow(core) - length(exposed) - length(targets) + event_count)
    result[, `:=`(analysis = analysis, stratum = stratum, drug = group)]
    result
  }))
}

summarize_numeric <- function(value) {
  value <- value[is.finite(value)]
  data.table::data.table(n_observed = length(value),
    mean = if (length(value)) mean(value) else NA_real_,
    sd = if (length(value) > 1L) stats::sd(value) else NA_real_,
    median = if (length(value)) stats::median(value) else NA_real_,
    q1 = if (length(value)) unname(stats::quantile(value, 0.25)) else NA_real_,
    q3 = if (length(value)) unname(stats::quantile(value, 0.75)) else NA_real_,
    min = if (length(value)) min(value) else NA_real_,
    max = if (length(value)) max(value) else NA_real_)
}

weibull_fit <- function(days, minimum_n = 10L) {
  days <- days[is.finite(days) & days > 0]
  result <- data.table::data.table(n_fit = 0L, scale_days = NA_real_,
    scale_low = NA_real_, scale_high = NA_real_, shape = NA_real_, shape_low = NA_real_,
    shape_high = NA_real_, distribution_pattern = NA_character_, fit_status = "insufficient_data")
  if (length(days) < minimum_n || length(unique(days)) < 2L) return(list(summary = result, model = NULL))
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(
    survival::survreg(survival::Surv(days, rep(1L, length(days))) ~ 1, dist = "weibull"),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }), error = function(condition) condition)
  if (inherits(fit, "error") || length(warnings)) {
    result$fit_status <- if (inherits(fit, "error")) conditionMessage(fit) else paste(warnings, collapse = "; ")
    return(list(summary = result, model = NULL))
  }
  covariance <- stats::vcov(fit)
  shape_index <- match("Log(scale)", rownames(covariance))
  if (is.na(shape_index) || any(!is.finite(covariance)) || any(diag(covariance) <= 0)) {
    result$fit_status <- "invalid_covariance"
    return(list(summary = result, model = NULL))
  }
  log_scale <- unname(stats::coef(fit)[1L])
  log_shape <- -log(fit$scale)
  scale_error <- sqrt(covariance[1L, 1L])
  shape_error <- sqrt(covariance[shape_index, shape_index])
  z_value <- stats::qnorm(0.975)
  result[, `:=`(n_fit = length(days), scale_days = exp(log_scale), scale_low = exp(log_scale - z_value * scale_error),
    scale_high = exp(log_scale + z_value * scale_error), shape = exp(log_shape),
    shape_low = exp(log_shape - z_value * shape_error), shape_high = exp(log_shape + z_value * shape_error),
    fit_status = "ok")]
  result[, distribution_pattern := ifelse(shape_high < 1, "decreasing_pattern",
    ifelse(shape_low > 1, "increasing_pattern", "CI_includes_one"))]
  list(summary = result, model = fit)
}

draw_files <- function(context, name, width, height, draw) {
  for (extension in c("pdf", "png", "tiff")) {
    path <- file.path(context$out, "figures", paste0(name, ".", extension))
    if (extension == "pdf") grDevices::pdf(path, width = width, height = height, family = "Times")
    if (extension == "png") grDevices::png(path, width = width, height = height, units = "in", res = 600)
    if (extension == "tiff") grDevices::tiff(path, width = width, height = height, units = "in", res = 600, compression = "lzw")
    graphics::par(family = "serif")
    tryCatch(draw(), finally = grDevices::dev.off())
  }
}
