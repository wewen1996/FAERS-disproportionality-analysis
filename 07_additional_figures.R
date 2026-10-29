# Figures from existing summary tables; no raw FAERS data or new estimation.
# Design reference: Exploring potential associations ... (2025), Table 2,
# Figure 3 and the July 2023 temporal analysis (doi:10.1016/j.eclinm.2025.103385).
# Base R graphics keeps this compatible with the existing R 4.1 workflow.

additional_label <- function(x) {
  dictionary <- c(study_drug_group = "Study-drug group", albiglutide = "Albiglutide",
    dulaglutide = "Dulaglutide", exenatide = "Exenatide", liraglutide = "Liraglutide",
    lixisenatide = "Lixisenatide", semaglutide = "Semaglutide", tirzepatide = "Tirzepatide",
    metformin = "Metformin", empagliflozin = "Empagliflozin",
    sglt2_inhibitors = "SGLT2 inhibitors", orlistat = "Orlistat",
    all_reports = "All eligible reports", diabetes_recorded = "Recorded diabetes indication",
    weight_management_recorded = "Recorded weight-management indication",
    metabolic_or_weight_recorded = "Recorded diabetes or weight-management indication")
  value <- unname(dictionary[as.character(x)])
  value[is.na(value)] <- gsub("_", " ", as.character(x)[is.na(value)], fixed = TRUE)
  value
}

additional_date <- function(x) {
  # for (cutoff in Date_vector) in R loses the Date class; existing CSVs
  # therefore contain days since 1970. Also accept corrected ISO date CSVs.
  text <- as.character(x)
  numeric_date <- grepl("^[0-9]+([.]0+)?$", text)
  result <- as.Date(rep(NA_character_, length(text)))
  result[numeric_date] <- as.Date(as.numeric(text[numeric_date]), origin = "1970-01-01")
  result[!numeric_date] <- as.Date(text[!numeric_date], format = "%Y-%m-%d")
  if (anyNA(result)) stop("Invalid temporal cutoff date in the result tables.")
  result
}

additional_estimable <- function(x) {
  x$status == "ok" & is.finite(x$ROR) & is.finite(x$ROR_low) &
    is.finite(x$ROR_high) & x$ROR_low > 0 & x$ROR > 0 &
    x$ROR_low <= x$ROR & x$ROR <= x$ROR_high
}

additional_limits <- function(x) {
  good <- additional_estimable(x)
  values <- c(1, x$ROR_low[good], x$ROR_high[good])
  limits <- range(log10(values))
  if (diff(limits) < 0.3) limits <- mean(limits) + c(-0.3, 0.3)
  10^(limits + c(-1, 1) * max(0.08, 0.04 * diff(limits)))
}

additional_forest <- function(x, title, subtitle, series_note, limits) {
  n <- nrow(x)
  y <- rev(seq_len(n))
  good <- additional_estimable(x)
  correction <- !is.na(x$zero_cell_correction) & x$zero_cell_correction
  cols <- c("#0072B2", "#A64B00")[x$plot_series]
  pch <- c(21L, 24L)[x$plot_series]
  graphics::par(oma = c(4.3, 0.4, 4.6, 0.4), family = "sans", cex = 0.96)
  graphics::layout(matrix(1:3, nrow = 1), widths = c(4.4, 3.4, 5.6))
  panel <- function(log = FALSE) {
    graphics::par(mar = c(3.2, 0.1, 1.7, 0.1))
    graphics::plot.new()
    graphics::plot.window(xlim = if (log) limits else c(0, 1),
      ylim = c(0.5, n + 1), log = if (log) "x" else "", xaxs = "i", yaxs = "i")
    fill <- rep(c("#F1F5F9", "white"), length.out = length(unique(x$plot_group)))
    bands <- fill[match(x$plot_group, unique(x$plot_group))]
    graphics::rect(if (log) limits[1L] else 0, y - 0.5,
      if (log) limits[2L] else 1, y + 0.5, border = NA, col = bands)
  }
  panel()
  graphics::text(0.02, n + 0.75, "Study drug / reference", adj = 0, font = 2, cex = 0.86)
  graphics::text(0.02, y, x$plot_label, adj = 0, cex = 0.83)

  panel(TRUE)
  graphics::abline(v = 1, lty = 2, col = "#666666", lwd = 1.1)
  graphics::segments(x$ROR_low[good], y[good], x$ROR_high[good], y[good], col = cols[good], lwd = 1.8)
  graphics::segments(x$ROR_low[good], y[good] - 0.08, x$ROR_low[good], y[good] + 0.08, col = cols[good])
  graphics::segments(x$ROR_high[good], y[good] - 0.08, x$ROR_high[good], y[good] + 0.08, col = cols[good])
  graphics::points(x$ROR[good], y[good], pch = pch[good], bg = cols[good], col = cols[good], cex = 1.1)
  if (any(!good)) graphics::text(sqrt(prod(limits)), y[!good], "NE", cex = 0.8, col = "#666666")
  candidates <- sort(unique(c(10^seq(floor(log10(limits[1L])), ceiling(log10(limits[2L]))),
    as.vector(outer(c(2, 5), 10^seq(floor(log10(limits[1L])), ceiling(log10(limits[2L]))))))))
  ticks <- candidates[candidates >= limits[1L] & candidates <= limits[2L]]
  # Keep labels readable for very wide sparse-data intervals.
  if (length(ticks) > 8L) ticks <- ticks[seq(1L, length(ticks), length.out = 7L) |> round() |> unique()]
  graphics::axis(1, at = ticks, labels = format(ticks, scientific = FALSE, trim = TRUE), cex.axis = 0.78)
  graphics::mtext("Reporting odds ratio (log scale)", side = 1, line = 2, cex = 0.84)

  panel()
  graphics::text(c(0.20, 0.65, 0.95), n + 0.78,
    c("Events / reports", "ROR (95% CI)", "BH q"), font = 2, cex = 0.80)
  graphics::text(0.20, n + 0.46, "Study drug | Reference", cex = 0.67, col = "#555555")
  counts <- sprintf("%.0f/%.0f | %.0f/%.0f", x$a, x$a + x$b, x$c, x$c + x$d)
  graphics::text(0.20, y, counts, cex = 0.74)
  estimates <- rep("Not estimable", n)
  estimates[good] <- sprintf("%.2f (%.2f-%.2f)%s", x$ROR[good], x$ROR_low[good],
    x$ROR_high[good], ifelse(correction[good], "*", ""))
  graphics::text(0.65, y, estimates, cex = 0.78)
  q <- ifelse(is.na(x$fisher_q_BH), "--", ifelse(x$fisher_q_BH < 0.001, "<0.001",
    sprintf("%.3f", x$fisher_q_BH)))
  graphics::text(0.95, y, q, cex = 0.77)
  graphics::mtext(title, side = 3, outer = TRUE, line = 3.0, cex = 1.16, font = 2)
  graphics::mtext(subtitle, side = 3, outer = TRUE, line = 1.6, cex = 0.87)
  graphics::mtext(series_note, side = 3, outer = TRUE, line = 0.25, cex = 0.79)
  graphics::mtext("95% log-Wald CIs from saved results. * Zero-cell correction (+0.5); displayed counts are uncorrected.",
    side = 1, outer = TRUE, line = 0.5, cex = 0.76)
  graphics::mtext("NE: not estimable (e.g. zero margin). BH q: saved adjusted Fisher P; grouping follows the source analysis.",
    side = 1, outer = TRUE, line = 1.7, cex = 0.75)
  graphics::mtext("Spontaneous-report comparisons do not estimate incidence or causality; before/after significance is not an interaction test.",
    side = 1, outer = TRUE, line = 2.9, cex = 0.74)
}

plot_additional_results <- function(run_dir, output_dir = NULL,
                                    formats = c("pdf", "png", "tiff"), dpi = 600L) {
  formats <- unique(match.arg(formats, c("pdf", "png", "tiff"), several.ok = TRUE))
  stopifnot(length(dpi) == 1L, is.finite(dpi), dpi >= 72)
  files <- c(active = "active_comparator_signals.csv", temporal = "temporal_notoriety_signals.csv",
    temporal_active = "temporal_active_comparator_signals.csv")
  paths <- file.path(run_dir, "tables", files)
  if (any(!file.exists(paths))) stop("Missing plotting inputs: ", paste(paths[!file.exists(paths)], collapse = "; "))
  read_table <- function(path) {
    x <- utils::read.csv(path, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE, check.names = FALSE)
    required <- c("a", "b", "c", "d", "ROR", "ROR_low", "ROR_high", "status", "zero_cell_correction", "fisher_q_BH")
    if (!all(required %in% names(x))) stop("Incomplete result schema in ", path)
    for (column in c("a", "b", "c", "d", "ROR", "ROR_low", "ROR_high", "fisher_q_BH")) {
      x[[column]] <- as.numeric(x[[column]])
    }
    if (any(!is.finite(as.matrix(x[c("a", "b", "c", "d")])) | as.matrix(x[c("a", "b", "c", "d")]) < 0)) {
      stop("Invalid four-cell counts in ", path)
    }
    if (any(x$status == "ok" & !additional_estimable(x))) stop("Invalid saved ROR/CI in ", path)
    x
  }
  tables <- lapply(paths, read_table)
  names(tables) <- names(files)
  for (name in c("temporal", "temporal_active")) {
    tables[[name]]$cutoff <- additional_date(tables[[name]]$cutoff)
    if (any(!tables[[name]]$period %in% c("before", "after"))) stop("Unexpected temporal period label.")
  }
  if (is.null(output_dir)) output_dir <- file.path(run_dir, "figures",
    paste0("additional_", format(Sys.time(), "%Y%m%d_%H%M%S", tz = "Asia/Shanghai"), "_", Sys.getpid()))
  if (dir.exists(output_dir)) stop("Output folder exists; select a new folder to preserve prior figures: ", output_dir)
  if (!dir.create(output_dir, recursive = TRUE)) stop("Cannot create plotting output directory.")
  dir.create(file.path(output_dir, "source_data"))
  manifest <- list()
  save_plot <- function(x, name, title, subtitle, series_note, limits, source) {
    height <- max(5.4, 2.8 + 0.32 * nrow(x))
    utils::write.csv(x, file.path(output_dir, "source_data", paste0(name, ".csv")), row.names = FALSE, na = "")
    for (extension in formats) {
      path <- file.path(output_dir, paste0(name, ".", extension))
      if (extension == "pdf") grDevices::pdf(path, width = 13.4, height = height, family = "Helvetica", useDingbats = FALSE)
      if (extension == "png") grDevices::png(path, width = 13.4, height = height, units = "in", res = dpi, bg = "white")
      if (extension == "tiff") grDevices::tiff(path, width = 13.4, height = height, units = "in", res = dpi, compression = "lzw", bg = "white")
      tryCatch(additional_forest(x, title, subtitle, series_note, limits), finally = grDevices::dev.off())
      manifest[[length(manifest) + 1L]] <<- data.frame(figure = name, format = extension,
        width_inches = 13.4, height_inches = height, dpi = if (extension == "pdf") NA else dpi,
        source = source, displayed_rows = nrow(x), not_estimable_rows = sum(!additional_estimable(x)))
    }
    message("Plotted: ", name)
  }
  active <- tables$active
  required_active <- c("stratum", "target_drug", "comparator")
  if (!all(required_active %in% names(active))) stop("Active-comparator labels are missing.")
  targets <- unique(active$target_drug)
  comparators <- unique(active$comparator)
  strata <- unique(active$stratum)
  if (length(targets) > 2L) stop("This layout supports up to two active-comparator study drugs.")
  prepare_active <- function(x, temporal = FALSE) {
    order_id <- order(match(x$target_drug, targets), match(x$comparator, comparators),
      if (temporal) match(x$period, c("before", "after")) else seq_len(nrow(x)))
    x <- x[order_id, , drop = FALSE]
    x$plot_group <- paste(x$target_drug, x$comparator)
    x$plot_label <- paste(additional_label(x$target_drug), "vs", additional_label(x$comparator))
    x$plot_series <- if (temporal) match(x$period, c("before", "after")) else match(x$target_drug, targets)
    if (temporal) x$plot_label <- paste0(x$plot_label, " / ", ifelse(x$period == "before", "Before", "On/after"))
    x
  }
  active_limits <- additional_limits(active)
  for (stratum in strata) {
    x <- prepare_active(active[active$stratum == stratum, , drop = FALSE])
    save_plot(x, paste0("active_comparator_", stratum), "Active-comparator analysis: ischaemic colitis",
      additional_label(stratum), paste("Study-drug PS vs reference-drug PS; dual exposure excluded. Circle:",
        additional_label(targets[1L]), " | Triangle:", additional_label(targets[min(2L, length(targets))])),
      active_limits, files[["active"]])
  }
  # The saved pooled temporal rows in older runs were computed with drug ==
  # 'study_drug_group', which is absent from the exposure pairs. Do not plot that
  # erroneous zero as a measured class effect; retain all individual-drug rows.
  pooled <- tables$temporal[tables$temporal$drug == "study_drug_group", , drop = FALSE]
  temporal <- tables$temporal[tables$temporal$drug != "study_drug_group", , drop = FALSE]
  drugs <- unique(temporal$drug)
  temporal_limits <- additional_limits(temporal)
  temporal_active_limits <- additional_limits(tables$temporal_active)
  for (cutoff_text in unique(as.character(temporal$cutoff))) {
    cutoff <- as.Date(cutoff_text)
    x <- temporal[temporal$cutoff == cutoff, , drop = FALSE]
    x <- x[order(match(x$drug, drugs), match(x$period, c("before", "after"))), , drop = FALSE]
    x$plot_group <- x$drug
    x$plot_series <- match(x$period, c("before", "after"))
    x$plot_label <- paste0(additional_label(x$drug), " / ", ifelse(x$period == "before", "Before", "On/after"))
    node <- unique(x$cutoff_label)
    node_label <- switch(node[1L], reference_temporal_proxy = "Reference-paper attention proxy (exploratory for colitis)",
      sema_weight_approval = "Semaglutide weight-management approval proxy",
      tirz_weight_approval = "Tirzepatide weight-management approval proxy", gsub("_", " ", node[1L]))
    series_note <- "Circle: FDA receipt date < cutoff | Triangle: receipt date >= cutoff; separate, overlapping cutoff analyses."
    save_plot(x, paste0("temporal_notoriety_", cutoff_text),
      paste0("Before/after analysis: ", cutoff_text), node_label, series_note,
      temporal_limits, files[["temporal"]])
    x <- tables$temporal_active[tables$temporal_active$cutoff == cutoff, , drop = FALSE]
    if (nrow(x)) save_plot(prepare_active(x, TRUE), paste0("temporal_active_comparator_", cutoff_text),
      paste0("Active comparators before/after ", cutoff_text), paste(node_label, "- all eligible reports"),
      series_note, temporal_active_limits, files[["temporal_active"]])
  }
  utils::write.csv(do.call(rbind, manifest), file.path(output_dir, "figure_manifest.csv"), row.names = FALSE)
  utils::write.csv(data.frame(path = paths, md5 = unname(tools::md5sum(paths))),
    file.path(output_dir, "input_checksums.csv"), row.names = FALSE)
  if (nrow(pooled)) utils::write.csv(pooled, file.path(output_dir, "source_data", "pooled_temporal_rows_not_plotted.csv"), row.names = FALSE)
  notes <- c(
    "Supplemental forest plots: existing summary estimates, no new statistical fitting.",
    "Reference: doi:10.1016/j.eclinm.2025.103385, Table 2, Figure 3 and July 2023 temporal analysis.",
    "The paper's active comparisons are a table; these figures adapt its ROR forest-plot design.",
    "Outcome: COLITIS ISCHAEMIC. This is not the depressive-disorder endpoint in the reference paper.",
    "July 2023 is a reference-paper proxy concerning neuropsychiatric attention, not a validated colitis-specific event.",
    "Dates in old CSVs are converted from days since 1970-01-01 to ISO dates.",
    "Temporal panels are FDA receipt-date splits (< cutoff; >= cutoff), not event-onset date splits.",
    "Each cutoff partitions the same study data independently; these periods must not be added across figures.",
    "All individual-drug rows are shown, including NE. ROR=1 is the null reference line.",
    "A zero margin is NE; a zero cell with nonzero margins uses the saved +0.5 correction, marked *.",
    "Counts are original a/(a+b) and c/(c+d); ROR confidence intervals and BH q values are taken unchanged from CSVs.",
    "Common full-range log axes within each plot family; no CI clipping, no estimates invented for missing rows.",
    "The paper's IC was not computed by this project and is not invented for these figures.",
    "No formal before/after interaction test is available; visual significance differences do not establish an interaction.",
    "All_reports refers to the comparison eligibility stratum, not a population exposure denominator.",
    "Pooled temporal rows are excluded: older code yields zero exposure for the synthetic group label; see preserved rows.",
    "Only the source tables of the completed run are read. Original tables, workbook and previous figures are untouched.",
    "Export defaults: PDF vector, PNG/TIFF 600 dpi, TIFF LZW; final journal sizing remains to be selected.",
    "Software-method acknowledgement: Kassis T, Agarwal V, He Y, Patel D, Brueckner AM (2026).",
    "Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents. doi:10.48550/arXiv.2609.00065.")
  writeLines(notes, file.path(output_dir, "FIGURE_NOTES.txt"), useBytes = TRUE)
  file.copy(file.path("F:/GLP1-RA/new", "07_additional_figures.R"), output_dir)
  capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  message("Additional figures saved in: ", output_dir)
  invisible(output_dir)
}

write_additional_figures <- function(context) {
  if (isTRUE(context$config$make_figures)) plot_additional_results(context$out)
  invisible(NULL)
}
