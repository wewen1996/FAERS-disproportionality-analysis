write_flow_and_figures <- function(context, cohort, tto) {
  flow <- context$flow
  counts <- merge(cohort$counts, tto$summary[, .(drug, valid_nonnegative_TTO, same_day_reports,
    positive_TTO, excluded_date_reports, fit_status)], by = "drug", all.x = TRUE, sort = FALSE)
  write_result(context, counts, "flow_counts_with_TTO")
  pooled <- counts[drug == "study_drug_group"]
  stage_labels <- c("Supplied cleaned FAERS DEMO rows", "Reports within source-quarter window",
    "Rows with valid identifiers", "Unique report versions", "Latest report per CASEID",
    "Latest reports with a recorded reaction", "Eligible report-level background")
  nodes <- data.table::data.table(node_id = flow$stage, label = stage_labels,
    n = flow$retained_n, unit = flow$unit)
  nodes <- data.table::rbindlist(list(nodes,
    data.table::data.table(node_id = c("08_study_drug_PS", "09_target_PT", "10_TTO_valid", "11_TTO_positive"),
      label = c("Study-drug reports, primary suspect", context$config$target_pt,
        "Reports with valid nonnegative TTO", "Reports with positive TTO"),
      n = c(pooled$main_PS_reports, pooled$target_PT_reports, pooled$valid_nonnegative_TTO, pooled$positive_TTO),
      unit = "case reports")))
  edges <- data.table::data.table(from = head(nodes$node_id, -1L), to = tail(nodes$node_id, -1L),
    retained_n = tail(nodes$n, -1L), excluded_n = head(nodes$n, -1L) - tail(nodes$n, -1L))
  stopifnot(all(edges$excluded_n >= 0))
  write_result(context, nodes, "flow_nodes")
  write_result(context, edges, "flow_edges")
  labels <- paste0(nodes$label, ": n = ", format(nodes$n, big.mark = ",", scientific = FALSE, trim = TRUE))
  writeLines(enc2utf8(labels), file.path(context$out, "tables", "flow_labels.txt"), useBytes = TRUE)
  if (!context$config$make_figures) return(invisible(nodes))
  draw_files(context, "study_flowchart", 11, 11, function() {
    graphics::par(mar = rep(0.5, 4))
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
    box <- function(xcenter, ycenter, width, height, label, color = "#D8E1F5", size = 0.86) {
      graphics::rect(xcenter - width / 2, ycenter - height / 2,
        xcenter + width / 2, ycenter + height / 2, col = color, border = "#777777")
      graphics::text(xcenter, ycenter, label, cex = size)
    }
    connect <- function(y1, y2) graphics::arrows(0.38, y1, 0.38, y2, length = 0.08, col = "#777777")
    period <- paste(quarter_label(context$config$first_quarter), quarter_label(context$config$last_quarter), sep = " - ")
    box(0.38, 0.92, 0.54, 0.11, paste("Supplied cleaned FAERS dataset", period,
      paste0("DEMO rows: ", flow$retained_n[1L]), sep = "\n"))
    connect(0.865, 0.825)
    box(0.38, 0.775, 0.54, 0.10, paste0("Latest report per CASEID\nn = ", flow$retained_n[5L]))
    box(0.84, 0.78, 0.28, 0.16, paste0("Excluded DEMO rows\nOutside window: ", flow$excluded_n[2L],
      "\nMissing IDs: ", flow$excluded_n[3L], "\nRepeated PRIMARYID: ", flow$excluded_n[4L],
      "\nOlder case versions: ", flow$excluded_n[5L]), "#E4E4E4", 0.76)
    connect(0.725, 0.69)
    box(0.38, 0.64, 0.54, 0.10, paste0("Eligible report-level background\nn = ", nrow(cohort$core)))
    box(0.84, 0.64, 0.28, 0.10, paste0("Missing REAC: ", flow$excluded_n[6L],
      "\nMissing named DRUG: ", flow$excluded_n[7L]), "#E4E4E4", 0.78)
    connect(0.59, 0.555)
    drug_lines <- counts[drug != "study_drug_group", paste0(drug, ": ", main_PS_reports)]
    box(0.38, 0.465, 0.54, 0.18, paste(c(paste0("Primary-suspect study drugs; unique reports = ", pooled$main_PS_reports), drug_lines), collapse = "\n"), size = 0.78)
    box(0.84, 0.465, 0.28, 0.12, "Exact mapped drug names\nROLE_COD = PS\nDrug groups may overlap", "#FFD7B5", 0.8)
    connect(0.375, 0.335)
    target_lines <- counts[drug != "study_drug_group", paste0(drug, ": ", target_PT_reports)]
    box(0.38, 0.25, 0.54, 0.17, paste(c(paste0(context$config$target_pt, "; unique reports = ", pooled$target_PT_reports), target_lines), collapse = "\n"), size = 0.78)
    box(0.84, 0.25, 0.28, 0.12, paste0("TTO: valid n = ", pooled$valid_nonnegative_TTO,
      "\nPositive n = ", pooled$positive_TTO, "\nSame day n = ", pooled$same_day_reports), "#E4E4E4", 0.8)
    for (position in c(0.78, 0.64, 0.25)) graphics::arrows(0.65, position, 0.70, position, length = 0.06, col = "#777777")
    graphics::arrows(0.70, 0.465, 0.65, 0.465, length = 0.06, col = "#777777")
    centers <- c(0.13, 0.38, 0.63, 0.88)
    titles <- c("Descriptive\nanalysis", "Disproportionality\nanalysis", "Time to onset\nanalysis", "Co-medication\nanalysis")
    for (index in seq_along(centers)) {
      graphics::arrows(0.38, 0.165, centers[index], 0.11, length = 0.06, col = "#777777")
      box(centers[index], 0.065, 0.22, 0.09, titles[index], "#D2E5BC", 0.86)
    }
  })
  signal <- context$tables$signals_main[is.finite(ROR_low) & is.finite(ROR_high) & ROR_low > 0]
  if (nrow(signal)) draw_files(context, "ROR_forest", 9, 6, function() {
    graphics::par(mar = c(5, 10, 3, 2))
    limits <- range(c(1, signal$ROR_low, signal$ROR_high))
    graphics::plot(signal$ROR, seq_len(nrow(signal)), log = "x", xlim = limits,
      ylim = c(0.5, nrow(signal) + 0.5), yaxt = "n", xlab = "Reporting odds ratio (95% CI)", ylab = "", pch = 19)
    graphics::axis(2, at = seq_len(nrow(signal)), labels = signal$drug, las = 1, cex.axis = 0.85)
    graphics::segments(signal$ROR_low, seq_len(nrow(signal)), signal$ROR_high, seq_len(nrow(signal)))
    graphics::abline(v = 1, lty = 2, col = "grey50")
  })
  annual <- context$tables$annual_report_counts[drug == "study_drug_group"]
  draw_files(context, "annual_target_reports", 9, 5, function() {
    graphics::par(mar = c(6, 5, 3, 1))
    graphics::barplot(annual$target_reports, names.arg = paste0(annual$year, ifelse(annual$full_year_available, "", "*")),
      las = 2, col = ifelse(annual$full_year_available, "#5E81AC", "#D08770"),
      ylab = "Target-PT reports", main = "* Incomplete source year")
  })
  comed <- context$tables$comedication_frequencies
  if (!is.null(comed)) {
    comed <- comed[drug == "study_drug_group" & co_target_reports > 0][order(-co_target_reports)]
    comed <- head(comed, 20L)
    if (nrow(comed)) draw_files(context, "comedication_target_top20", 10, 7, function() {
      graphics::par(mar = c(5, 18, 3, 1))
      graphics::barplot(rev(comed$co_target_reports), names.arg = rev(comed$comed_name),
        horiz = TRUE, las = 1, cex.names = 0.65, col = "#88A97A",
        xlab = "Target-PT reports with a recorded concomitant drug (C)")
    })
  }
  for (group in c("study_drug_group", context$config$drugs)) {
    days <- tto$records[drug == group & status == "valid", TTO_DAYS]
    if (!length(days)) next
    fit <- tto$summary[drug == group]
    draw_files(context, paste0("TTO_", group), 8, 6, function() {
      graphics::par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
      graphics::boxplot(days, horizontal = TRUE, main = group, xlab = "Days from drug start to reported event")
      breaks <- seq(0, max(days) + 20, by = 20)
      graphics::hist(days, breaks = breaks, include.lowest = TRUE, col = "#88C0D0",
        xlab = "Days (same-day reports retained)", main = "Observed reported times; 20-day bins")
      if (fit$fit_status == "ok") {
        grid <- seq(0.1, max(days) + 20, length.out = 500)
        graphics::lines(grid, stats::dweibull(grid, shape = fit$shape, scale = fit$scale_days) * fit$n_fit * 20,
          col = "#BF616A", lwd = 2)
      }
    })
  }
  invisible(nodes)
}
