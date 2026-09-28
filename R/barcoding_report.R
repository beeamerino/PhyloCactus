# ------------------------------------------------------------------------------
# Molecular diagnostic branch (11_barcoding/), the identification report. PhyloCactus 0.5.0.
#
# For every query, one self-contained HTML written by the package (no pandoc) and a set of CSV
# files: the answer per locus, the alternatives under the threshold, a descriptive reading across
# loci (never a combined call, E8), the sampling of the named genera in the library, and the
# provenance of the run (input, library, models, software). The figures follow Zeng et al. (2026,
# New Phytologist): the route of the query, the answer per locus and the sampling of the genus.
#
# Governing documents: 01_minutes_and_phases/04_barcoding_0.5.0/ of the audit repository, Phase 9
# proposal of 2026-09-28 (decisions IR1 to IR8 of BMM) and the review of maturity after Phase 8.
# ------------------------------------------------------------------------------

#' Base64 of a raw vector (RFC 4648), for the images embedded in the report
#' @noRd
.bc_base64 <- function(x) {
  n <- length(x)
  if (n == 0L) return("")
  chars <- c(LETTERS, letters, 0:9, "+", "/")
  pad <- (3L - n %% 3L) %% 3L
  v <- c(as.integer(x), rep(0L, pad))
  m <- matrix(v, nrow = 3L)
  w <- m[1L, ] * 65536 + m[2L, ] * 256 + m[3L, ]
  idx <- rbind(w %/% 262144, (w %/% 4096) %% 64, (w %/% 64) %% 64, w %% 64)
  s <- chars[as.vector(idx) + 1L]
  if (pad > 0L) s[(length(s) - pad + 1L):length(s)] <- "="
  paste(s, collapse = "")
}

#' Text safe inside HTML
#' @noRd
.bc_html_escape <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}

#' A data frame as an HTML table
#' @noRd
.bc_html_table <- function(d, class = "tab") {
  if (is.null(d) || nrow(d) == 0L) return("<p>None.</p>")
  fmt <- function(v) if (is.numeric(v)) ifelse(is.na(v), "", format(round(v, 3), trim = TRUE)) else .bc_html_escape(v)
  head <- paste0("<tr>", paste0("<th>", .bc_html_escape(names(d)), "</th>", collapse = ""), "</tr>")
  cols <- lapply(d, fmt)
  body <- vapply(seq_len(nrow(d)), function(i) {
    paste0("<tr>", paste0("<td>", vapply(cols, `[`, "", i), "</td>", collapse = ""), "</tr>")
  }, "")
  paste0("<div style=\"overflow-x:auto\"><table class=\"", class, "\">", head, paste(body, collapse = ""), "</table></div>")
}

#' A ggplot as an embedded PNG
#' @noRd
.bc_html_plot <- function(p, width = 8, height = 4) {
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f), add = TRUE)
  type <- if (isTRUE(capabilities("cairo"))) "cairo" else getOption("bitmapType")
  grDevices::png(f, width = width, height = height, units = "in", res = 110, type = type)
  print(p)
  grDevices::dev.off()
  paste0("<img alt=\"figure\" src=\"data:image/png;base64,", .bc_base64(readBin(f, "raw", file.size(f))), "\"/>")
}

#' Colours of the three states, used in every panel
#' @noRd
.bc_state_colours <- function() c("1" = "#1B7837", "2" = "#2166AC", "3" = "#9E9E9E")

#' Nuclear loci of the library; every other locus is plastid
#' @noRd
.bc_nuclear_loci <- function() c("ITS", "phyC", "pepC_like")

#' Versions of the software of a run
#' @noRd
.bc_software_record <- function(tools = character(0)) {
  v <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  sha <- tryCatch(utils::packageDescription("PhyloCactus")$RemoteSha, error = function(e) NULL)
  d <- data.frame(kind = "software", name = c("R", "PhyloCactus", "DECIPHER", "Biostrings"),
                  value = c(R.version.string, v("PhyloCactus"), v("DECIPHER"), v("Biostrings")), stringsAsFactors = FALSE)
  if (!is.null(sha) && !is.na(sha)) d <- rbind(d, data.frame(kind = "software", name = "PhyloCactus_commit", value = sha))
  if (length(tools)) d <- rbind(d, data.frame(kind = "software", name = names(tools), value = unname(tools)))
  d
}

#' First line of the version of an external tool
#' @noRd
.bc_tool_version <- function(exe, arg = "--version") {
  out <- tryCatch(suppressWarnings(system2(exe, arg, stdout = TRUE, stderr = TRUE)), error = function(e) character(0))
  out <- trimws(out[nzchar(trimws(out))])
  if (length(out)) out[1] else NA_character_
}

#' The run record of one identification (IR1), written as RUN_<run_name>.csv
#' @noRd
.bc_write_run_record <- function(dir, run_name, input_path, input_md5, input_type, route, library_dir, loci,
                                 models_dir, threshold, seed, min_overlap, tools = character(0), extra = NULL,
                                 loci_declared = NULL) {
  lib_files <- file.path(library_dir, paste0("LIB_", loci, ".fasta"))
  mod_files <- file.path(models_dir, paste0("MODEL_idtaxa_", loci, ".rds"))
  mod_files <- mod_files[file.exists(mod_files)]
  d <- rbind(
    data.frame(kind = "input", name = c("path", "md5", "input_type"), value = c(input_path, input_md5, input_type)),
    data.frame(kind = "route", name = "steps", value = paste(route, collapse = "; ")),
    data.frame(kind = "setting", name = c("run_name", "threshold", "seed", "min_overlap", "loci_declared", "date"),
               value = c(run_name, format(threshold), format(seed), format(min_overlap),
                         if (length(loci_declared)) paste(loci_declared, collapse = ", ") else "none",
                         format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))),
    data.frame(kind = "machine", name = c("machine", "operating_system", "platform"),
               value = c(unname(Sys.info()[["nodename"]]), paste(Sys.info()[["sysname"]], Sys.info()[["release"]]),
                         R.version$platform)),
    .bc_software_record(tools),
    data.frame(kind = "library", name = basename(lib_files), value = unname(tools::md5sum(lib_files))),
    if (length(mod_files)) data.frame(kind = "model", name = basename(mod_files), value = unname(tools::md5sum(mod_files))),
    extra)
  rownames(d) <- NULL
  utils::write.csv(d, file.path(dir, paste0("RUN_", run_name, ".csv")), row.names = FALSE)
  invisible(d)
}

#' Per-genus discrimination of each locus in the 6B folds
#'
#' From the IdTaxa predictions of the validation folds (step 7), cut at the operating threshold of
#' step 9 (60 by default), counts per locus, scheme and true genus the queries whose species was
#' named correctly, the wrong species, the right and wrong genera, and the unassigned, with their
#' rates. It says how well each locus tells the species of one genus apart in the library, which the
#' identification report shows beside the answer (decision IR4 of BMM, 28-09).
#'
#' @param classifier_dir Character. Directory of step 7, with
#'   `TABLE_barcoding_predictions_{species,genus}_idtaxa.csv`.
#' @param threshold_dir Character. Directory of step 9, with
#'   `TABLE_barcoding_threshold_operating_idtaxa.csv`.
#' @param output_dir Character. Where `TABLE_barcoding_genus_discrimination_idtaxa.csv` is written.
#' @param idtaxa_rule Character. Row of the operating table used, `"default_60"` by default.
#' @return Invisibly, the table.
#' @seealso [summarise_barcoding_metrics()], [report_barcoding_identification()].
#' @examples
#' \dontrun{
#' summarise_barcoding_genus_discrimination()
#' }
#' @export
summarise_barcoding_genus_discrimination <- function(classifier_dir = file.path("11_barcoding", "7_classifier"),
                                                     threshold_dir = file.path("11_barcoding", "9_threshold"),
                                                     output_dir = file.path("11_barcoding", "11_metrics"),
                                                     idtaxa_rule = "default_60") {
  .bc_assert_output_dir(output_dir)
  f_pred <- file.path(classifier_dir, sprintf("TABLE_barcoding_predictions_%s_idtaxa.csv", c("species", "genus")))
  if (!all(file.exists(f_pred))) {
    stop("No ", basename(f_pred[!file.exists(f_pred)][1]), " in ", classifier_dir, ". Run classify_barcoding_folds() first.",
         call. = FALSE)
  }
  f_op <- file.path(threshold_dir, "TABLE_barcoding_threshold_operating_idtaxa.csv")
  if (!file.exists(f_op)) stop("No ", basename(f_op), " in ", threshold_dir, ". Run sweep_barcoding_threshold() first.", call. = FALSE)
  pred <- do.call(rbind, lapply(f_pred, utils::read.csv, stringsAsFactors = FALSE))
  p <- .bc_metrics_operating(pred, utils::read.csv(f_op, stringsAsFactors = FALSE), "idtaxa", idtaxa_rule)
  p$category <- NA_character_
  for (sc in unique(p$scheme)) p$category[p$scheme == sc] <- .bc_outcome_category(p[p$scheme == sc, , drop = FALSE], sc)
  key <- paste(p$locus, p$scheme, p$true_genus, sep = "\r")
  rows <- lapply(split(seq_len(nrow(p)), key), function(i) {
    x <- p[i, , drop = FALSE]
    n <- nrow(x)
    cnt <- vapply(.bc_outcome_levels(), function(k) sum(x$category == k), integer(1))
    data.frame(locus = x$locus[1], scheme = x$scheme[1], genus = x$true_genus[1], threshold = x$operating_threshold[1],
               queries = n, correct_species = cnt[["correct_species"]], wrong_species = cnt[["wrong_species"]],
               correct_genus = cnt[["correct_genus"]], wrong_genus = cnt[["wrong_genus"]], unassigned = cnt[["unassigned"]],
               correct_species_rate = cnt[["correct_species"]] / n, wrong_species_rate = cnt[["wrong_species"]] / n,
               correct_genus_rate = cnt[["correct_genus"]] / n, unassigned_rate = cnt[["unassigned"]] / n,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$locus, out$scheme, out$genus, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(output_dir, "TABLE_barcoding_genus_discrimination_idtaxa.csv")
  utils::write.csv(out, f, row.names = FALSE)
  message("Per-genus discrimination: ", nrow(out), " rows (locus, scheme, genus) written to ", f, ".")
  invisible(out)
}

#' Library counts per species and locus, from the headers of LIB_<locus>.fasta
#' @noRd
.bc_report_library_counts <- function(library_dir) {
  files <- list.files(library_dir, pattern = "^LIB_.*\\.fasta$", full.names = TRUE)
  do.call(rbind, lapply(files, function(f) {
    h <- names(Biostrings::readDNAStringSet(f))
    p <- .bc_parse_header(h)
    data.frame(locus = sub("^LIB_(.*)\\.fasta$", "\\1", basename(f)), species = p$species, sid = p$sid,
               source = ifelse(grepl("__", p$sid, fixed = TRUE), "genbank_plastome", "phylotaR"), stringsAsFactors = FALSE)
  }))
}

#' Candidates of state 2, one row per locus, for the HTML
#' @noRd
.bc_report_candidates_summary <- function(al) {
  c2 <- al[al$kind == "candidate_of_named_genus", , drop = FALSE]
  if (!nrow(c2)) return(NULL)
  key <- paste(if ("query" %in% names(c2)) c2$query else "", c2$locus, sep = "\r")
  do.call(rbind, lapply(split(seq_len(nrow(c2)), factor(key, levels = unique(key))), function(i) {
    x <- c2[i, , drop = FALSE]
    data.frame(locus = x$locus[1], genus = x$genus[1], species_in_library = nrow(x),
               species = paste(sort(x$species, method = "radix"), collapse = ", "), stringsAsFactors = FALSE)
  }))
}

#' Route diagram of one query (Fig. 1 of Zeng et al. 2026, as run for this query)
#' @noRd
.bc_report_route_plot <- function(input_type, t) {
  first <- switch(input_type,
                  assembly = c("Assembly", "minimap2", "samtools faidx"),
                  reads = c("Paired reads", "fastp", "GetOrganelle"),
                  genbank = c("GenBank record", "Strand rule"),
                  c("Sequences", "Strand rule"))
  has_region <- any(!t$path %in% "none")
  steps <- c(first, "Crop", "Cut to locus core (J1)", "IdTaxa", "Cut to genus core (J3b)", "State per locus")
  run <- c(rep(TRUE, length(first)), any(t$path %in% "crop"), any(t$core_trimmed > 0 & !is.na(t$core_trimmed)),
           has_region, any(t$genus_core_trimmed > 0 & !is.na(t$genus_core_trimmed)), TRUE)
  d <- data.frame(x = seq_along(steps), step = steps, run = ifelse(run, "run", "not run"), stringsAsFactors = FALSE)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = 1)) +
    ggplot2::geom_segment(data = d[-nrow(d), ], ggplot2::aes(x = .data$x + 0.42, xend = .data$x + 0.58, y = 1, yend = 1),
                          arrow = grid::arrow(length = grid::unit(0.08, "in"))) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$run), width = 0.8, height = 0.6, colour = "grey30") +
    ggplot2::geom_text(ggplot2::aes(label = gsub(" ", "\n", .data$step)), size = 2.6) +
    ggplot2::scale_fill_manual(values = c(run = "#D9F0D3", `not run` = "#F0F0F0"), name = NULL) +
    ggplot2::theme_void() + ggplot2::theme(legend.position = "bottom")
}

#' Answer per locus: species as a solid bar, genus as a striped bar, colour by state (Fig. 3c style)
#'
#' The two ranks are told apart by the fill (solid, striped), not by a shade, so that no rank can be
#' read as a state (BMM, 28-09). Drawn with rectangles and segments: no pattern package is needed.
#' @noRd
.bc_report_answer_plot <- function(t, threshold) {
  t <- t[order(t$compartment != "plastid", t$locus, method = "radix"), , drop = FALSE]
  # One position per locus across both panels, so that no break of one panel carries a label of the other
  t$i <- seq_len(nrow(t))
  t$state <- factor(as.character(t$state), levels = c("1", "2", "3"))
  w <- 0.36
  bars <- rbind(data.frame(compartment = t$compartment, state = t$state, rank = "genus", xmin = t$i - w - 0.02, xmax = t$i - 0.02,
                           ymax = ifelse(is.na(t$genus_confidence), 0, t$genus_confidence)),
                data.frame(compartment = t$compartment, state = t$state, rank = "species", xmin = t$i + 0.02, xmax = t$i + w + 0.02,
                           ymax = ifelse(is.na(t$species_confidence), 0, t$species_confidence)))
  gen <- bars[bars$rank == "genus" & bars$ymax > 0, , drop = FALSE]
  stripes <- do.call(rbind, c(list(data.frame(compartment = character(0), state = factor(character(0), levels = c("1", "2", "3")),
                                              x = numeric(0), xend = numeric(0), y = numeric(0))),
                              lapply(seq_len(nrow(gen)), function(k) {
    if (gen$ymax[k] < 4) return(NULL)
    y <- seq(4, gen$ymax[k], by = 4)
    data.frame(compartment = gen$compartment[k], state = gen$state[k], x = gen$xmin[k], xend = gen$xmax[k], y = y)
  })))
  lab <- unique(t[, c("compartment", "i", "locus")])
  # Zero-height rectangles of the three states, so that the legend always shows all three
  keys <- data.frame(compartment = t$compartment[1], state = factor(c("1", "2", "3"), levels = c("1", "2", "3")),
                     xmin = t$i[1], xmax = t$i[1], ymax = 0)
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = keys, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0, ymax = .data$ymax,
                                                 fill = .data$state)) +
    ggplot2::geom_rect(data = bars[bars$rank == "species", ], ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0,
                                                                           ymax = .data$ymax, fill = .data$state), colour = "grey20") +
    ggplot2::geom_rect(data = gen, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0, ymax = .data$ymax,
                                                colour = .data$state), fill = "white", linewidth = 0.6) +
    ggplot2::geom_segment(data = stripes, ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y,
                                                       colour = .data$state), linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed") +
    ggplot2::scale_fill_manual(values = .bc_state_colours(), name = "State", drop = FALSE) +
    ggplot2::scale_colour_manual(values = .bc_state_colours(), guide = "none", drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = lab$i, labels = lab$locus, expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::coord_cartesian(ylim = c(0, 100)) +
    ggplot2::facet_grid(. ~ compartment, scales = "free_x", space = "free_x") +
    ggplot2::labs(x = NULL, y = "IdTaxa confidence", caption = "Striped bar: genus. Solid bar: species. Dashed line: threshold.") +
    ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Sampling of the named genera in the library: species by locus (Fig. 4c style)
#' @noRd
.bc_report_sampling_plot <- function(s, named_species) {
  s$label <- ifelse(s$sequences > 0, s$sequences, "")
  s$named <- s$species %in% named_species
  ggplot2::ggplot(s, ggplot2::aes(x = .data$locus, y = .data$species)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$sequences > 0), colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$label, fontface = ifelse(.data$named, "bold", "plain")), size = 2.6) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#A6DBA0", `FALSE` = "#F7F7F7"), labels = c(`TRUE` = "sequences",
                               `FALSE` = "none"), name = NULL) +
    ggplot2::facet_grid(genus ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                                         strip.text.y = ggplot2::element_text(angle = 0, face = "italic"))
}

#' The kind of data of a query, named at the start of the report
#' @noRd
.bc_report_data_kind <- function(input_type, query_length, query = NULL) {
  if (identical(input_type, "assembly")) return("Nuclear genome assembly (FASTA); loci located with minimap2 and cut with samtools")
  if (identical(input_type, "reads")) {
    return("Genome skimming or whole-genome sequencing reads (Illumina, paired-end); plastome and nrDNA assembled with GetOrganelle")
  }
  if (!is.null(query) && length(unique(query)) > 1L) {
    len <- tapply(query_length, query, max)
    return(sprintf("%d sequences, %s%s", length(len),
                   if (min(len) == max(len)) sprintf("%d bases each", min(len)) else sprintf("%d to %d bases", min(len), max(len)),
                   if (identical(input_type, "genbank")) ", read from a GenBank flat file" else ""))
  }
  n <- suppressWarnings(max(query_length, na.rm = TRUE))
  kind <- if (is.finite(n) && n >= 100000) sprintf("Complete plastome (%s bases)", format(n, big.mark = " "))
          else if (is.finite(n) && n >= 5000) sprintf("Contig or partial plastome (%s bases)", format(n, big.mark = " "))
          else sprintf("Sanger-type sequence (%s bases; one or a few loci)", if (is.finite(n)) n else "unknown")
  if (identical(input_type, "genbank")) paste0(kind, ", read from a GenBank flat file") else kind
}

#' The sources of the library, in words
#' @noRd
.bc_report_library_words <- function(libc) {
  ex <- libc[libc$source == "genbank_plastome", , drop = FALSE]
  w <- "Sanger records of GenBank (phylotaR)"
  if (nrow(ex)) w <- paste0(w, " and ", nrow(ex), " loci cut from ", length(unique(sub("__.*$", "", ex$sid))), " GenBank plastomes")
  w
}

#' Banner of the report: the logo of the package and its name
#' @noRd
.bc_report_banner <- function(subtitle) {
  f <- system.file("report", "logo.png", package = "PhyloCactus")
  logo <- if (nzchar(f)) paste0("<img alt=\"PhyloCactus logo\" src=\"data:image/png;base64,",
                                .bc_base64(readBin(f, "raw", file.size(f))), "\"/>") else ""
  paste0("<div class=\"banner\">", logo, "<div><h1>PhyloCactus \U0001F335</h1>",
         "<p><strong>Molecular identification report</strong></p><p>", .bc_html_escape(subtitle), "</p></div></div>")
}

#' The loci found in the query, and whether they were declared or detected
#' @noRd
.bc_report_loci_found <- function(t, declared) {
  found <- sort(unique(t$locus[!t$reason %in% c("no_overlap", "assembly_failed", "single_species_library")]), method = "radix")
  how <- if (is.na(declared) || identical(declared, "none")) {
    "detected by overlap with the library; no locus was declared, so every locus of the library was tested"
  } else paste0("declared by the user: ", declared)
  paste0(if (length(found)) paste(found, collapse = ", ") else "none", " (", how, ")")
}

#' The answers of a sample, counted over the loci found, and the loci not found
#' @noRd
.bc_report_answers_words <- function(t) {
  f <- t[!t$reason %in% c("no_overlap", "assembly_failed", "single_species_library"), , drop = FALSE]
  nf <- setdiff(sort(unique(t$locus), method = "radix"), f$locus)
  paste0(sprintf("%d loci found: %d named to species (state 1), %d to genus (state 2), %d not assignable (state 3)",
                 nrow(f), sum(f$state == 1L), sum(f$state == 2L), sum(f$state == 3L)),
         "; not found: ", if (length(nf)) paste(nf, collapse = ", ") else "none")
}

#' Footer of the report: version of PhyloCactus and how to cite it
#' @noRd
.bc_report_footer <- function(run_date) {
  v <- tryCatch(as.character(utils::packageVersion("PhyloCactus")), error = function(e) "unknown")
  year <- if (!is.na(run_date) && grepl("^[0-9]{4}", run_date)) substr(run_date, 1, 4) else format(Sys.Date(), "%Y")
  cit <- tryCatch(suppressWarnings(format(utils::citation("PhyloCactus"), style = "text")), error = function(e) character(0))
  cit <- if (length(cit)) paste(cit, collapse = " ") else
    paste0("Meri\u00f1o, B. M. (", year, "). PhyloCactus. R package version ", v, ". https://github.com/beeamerino/PhyloCactus")
  cit <- gsub("(????)", paste0("(", year, ")"), cit, fixed = TRUE)
  cit <- gsub("\\s+", " ", gsub("[_<>]", "", cit))
  paste0("<footer><p>Report written by PhyloCactus ", .bc_html_escape(v), " \U0001F335</p>",
         "<p><strong>How to cite:</strong> ", .bc_html_escape(cit),
         " Until a publication describes PhyloCactus, cite the repository as above.</p></footer>")
}

#' Build the identification report of one run
#'
#' Reads `TABLE_barcoding_identify_<run_name>.csv` and `RUN_<run_name>.csv` from `results_dir` and
#' writes, for each query, `REPORT_<run_name>_<query>.html`, a self-contained page written by the
#' package (no `pandoc`); an assembly, a set of reads, or several sequences declared as one sample
#' (argument `sample` of [identify_barcoding_query()]) get one `REPORT_<run_name>.html` that reads
#' them together. For the run, a set of CSV files: `REPORT_<run_name>_answer.csv`,
#' `_alternatives.csv`, `_agreement.csv`, `_sampling.csv` and `_provenance.csv`, and
#' `REPORT_<run_name>_checksums.csv` with the md5 of the table, the run record and every file written.
#'
#' The report has nine sections: (1) query and route; (2) answer per locus; (3) reading across loci,
#' descriptive only, since the loci are never combined into one call (rule E8 of the validation
#' plan); (4) alternatives under the threshold, marked not assigned: IdTaxa gives one path per query,
#' so the alternatives are that path when it falls under the threshold and, in state 2, the species
#' of the named genus in the library; (5) how far to trust each locus: the 6D wrong-species rates and
#' the per-genus discrimination of [summarise_barcoding_genus_discrimination()]; (6) sampling of the
#' named genera in the library; (7) reads; (8) library, models and software; (9) limits. The figures
#' follow Zeng et al. (2026, New Phytologist). The three identification functions call it at the end;
#' a report can be rebuilt from saved tables.
#'
#' @param run_name Character. The run.
#' @param results_dir Character. Where the table and the run record are.
#' @param library_dir Character. Directory of step 4, with the `LIB_<locus>.fasta` files.
#' @param metrics_dir Character or `NULL`. Directory of step 11, where
#'   `TABLE_barcoding_genus_discrimination_idtaxa.csv` is looked for.
#' @param output_dir Character. Where the report is written; `results_dir` by default.
#' @return Invisibly, the paths of the HTML reports.
#' @seealso [identify_barcoding_query()], [identify_barcoding_assembly()], [identify_barcoding_reads()].
#' @examples
#' \dontrun{
#' report_barcoding_identification("query")
#' }
#' @export
report_barcoding_identification <- function(run_name,
                                            results_dir = file.path("11_barcoding", "10_identify"),
                                            library_dir = file.path("11_barcoding", "4_library"),
                                            metrics_dir = file.path("11_barcoding", "11_metrics"),
                                            output_dir = results_dir) {
  f_tab <- file.path(results_dir, paste0("TABLE_barcoding_identify_", run_name, ".csv"))
  f_rec <- file.path(results_dir, paste0("RUN_", run_name, ".csv"))
  for (f in c(f_tab, f_rec)) if (!file.exists(f)) stop("Not found: ", f, ".", call. = FALSE)
  tab <- utils::read.csv(f_tab, stringsAsFactors = FALSE)
  rec <- utils::read.csv(f_rec, stringsAsFactors = FALSE, colClasses = "character")
  get <- function(n) { v <- rec$value[rec$name == n]; if (length(v)) v[1] else NA_character_ }
  input_type <- get("input_type")
  threshold <- suppressWarnings(as.numeric(get("threshold")))
  if (is.na(threshold)) threshold <- tab$threshold[1]
  tab$compartment <- ifelse(tab$locus %in% .bc_nuclear_loci(), "nuclear", "plastid")
  # One report per sample: a query of sequences is a sample; an assembly or a set of reads is one
  # sample whose regions (scaffolds, contigs) are rows of the same table
  sample_name <- get("sample")
  per_run <- input_type %in% c("assembly", "reads") || !is.na(sample_name)
  tab$unit <- if (!is.na(sample_name)) sample_name else if (per_run) run_name else tab$query
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Answer
  if (!"identical_to_library" %in% names(tab)) tab$identical_to_library <- NA_character_
  extra <- intersect(c("identical_to_library", "scaffold", "scaffold_length", "mapq", "scaffolds_hit", "low_mapq", "contig", "organelle"),
                     names(tab))
  answer <- tab[, c("query", "locus", "compartment", "state", "predicted_species", "predicted_genus", "genus_confidence",
                    "species_confidence", "reason", "region_length", "validation_ws_rate_species_present",
                    "validation_ws_rate_species_absent", extra)]

  # Alternatives, never assigned
  alt <- list()
  for (i in seq_len(nrow(tab))) {
    r <- tab[i, ]
    if (r$state %in% 2:3 && !is.na(r$species_idtaxa) && (is.na(r$species_confidence) || r$species_confidence < threshold)) {
      alt[[length(alt) + 1L]] <- data.frame(query = r$query, locus = r$locus, kind = "best_path_under_threshold",
                                            genus = r$genus_idtaxa, species = r$species_idtaxa,
                                            genus_confidence = r$genus_confidence, species_confidence = r$species_confidence,
                                            assigned = FALSE, stringsAsFactors = FALSE)
    }
    if (r$state %in% 2L && !is.na(r$candidates) && nzchar(r$candidates)) {
      cs <- strsplit(r$candidates, "|", fixed = TRUE)[[1]]
      alt[[length(alt) + 1L]] <- data.frame(query = r$query, locus = r$locus, kind = "candidate_of_named_genus",
                                            genus = r$predicted_genus, species = cs, genus_confidence = NA_real_,
                                            species_confidence = NA_real_, assigned = FALSE, stringsAsFactors = FALSE)
    }
  }
  alternatives <- if (length(alt)) do.call(rbind, alt) else
    data.frame(query = character(0), locus = character(0), kind = character(0), genus = character(0), species = character(0),
               genus_confidence = numeric(0), species_confidence = numeric(0), assigned = logical(0))

  # Reading across loci, descriptive
  agree <- list()
  for (q in unique(tab$unit)) {
    x <- tab[tab$unit == q, , drop = FALSE]
    for (rk in c("genus", "species")) {
      col <- paste0(rk, "_idtaxa"); cf <- paste0(rk, "_confidence")
      for (tx in sort(unique(x[[col]][!is.na(x[[col]])]), method = "radix")) {
        k <- which(x[[col]] %in% tx)
        assigned <- if (rk == "genus") x$state[k] %in% 1:2 & x$predicted_genus[k] %in% tx else x$state[k] %in% 1L & x$predicted_species[k] %in% tx
        agree[[length(agree) + 1L]] <- data.frame(query = q, rank = rk, taxon = tx, n_loci = length(k), n_assigned = sum(assigned),
                                                  loci = paste(sort(x$locus[k], method = "radix"), collapse = "; "),
                                                  max_confidence = max(x[[cf]][k], na.rm = TRUE), stringsAsFactors = FALSE)
      }
    }
    # Loci with no answer in any sequence of the sample, each counted once
    none <- setdiff(sort(unique(x$locus), method = "radix"), x$locus[!is.na(x$genus_idtaxa)])
    if (length(none)) {
      agree[[length(agree) + 1L]] <- data.frame(query = q, rank = "none", taxon = NA_character_, n_loci = length(none), n_assigned = 0L,
                                                loci = paste(none, collapse = "; "), max_confidence = NA_real_, stringsAsFactors = FALSE)
    }
  }
  agreement <- do.call(rbind, agree)

  # Sampling of the named genera in the library
  libc <- .bc_report_library_counts(library_dir)
  genera <- sort(unique(tab$genus_idtaxa[!is.na(tab$genus_idtaxa)]), method = "radix")
  loci_lib <- sort(unique(libc$locus), method = "radix")
  sampling <- do.call(rbind, c(list(data.frame(genus = character(0), species = character(0), locus = character(0),
                                               sequences = integer(0))),
                               lapply(genera, function(g) {
    sp <- sort(unique(libc$species[sub("_.*$", "", libc$species) == g]), method = "radix")
    if (!length(sp)) return(NULL)
    grid <- expand.grid(locus = loci_lib, species = sp, stringsAsFactors = FALSE)
    grid$sequences <- vapply(seq_len(nrow(grid)), function(i) sum(libc$species == grid$species[i] & libc$locus == grid$locus[i]), integer(1))
    data.frame(genus = g, species = grid$species, locus = grid$locus, sequences = grid$sequences, stringsAsFactors = FALSE)
  })))

  # Per-genus discrimination, if written
  f_gd <- if (is.null(metrics_dir)) NA_character_ else file.path(metrics_dir, "TABLE_barcoding_genus_discrimination_idtaxa.csv")
  gd <- if (!is.na(f_gd) && file.exists(f_gd)) utils::read.csv(f_gd, stringsAsFactors = FALSE) else NULL

  provenance <- rec
  base <- file.path(output_dir, paste0("REPORT_", run_name, "_"))
  utils::write.csv(answer, paste0(base, "answer.csv"), row.names = FALSE)
  utils::write.csv(alternatives, paste0(base, "alternatives.csv"), row.names = FALSE)
  utils::write.csv(agreement, paste0(base, "agreement.csv"), row.names = FALSE)
  utils::write.csv(sampling, paste0(base, "sampling.csv"), row.names = FALSE)
  utils::write.csv(provenance, paste0(base, "provenance.csv"), row.names = FALSE)

  css <- paste0("<style>body{font-family:Helvetica,Arial,sans-serif;max-width:60em;margin:2em auto;padding:0 1em;color:#222}",
                "h1{font-size:1.5em}h2{font-size:1.15em;border-bottom:1px solid #ccc;margin-top:1.6em}",
                "table{border-collapse:collapse;font-size:0.85em;margin:0.5em 0}th,td{border:1px solid #ccc;padding:2px 6px;text-align:left}",
                "th{background:#f3f3f3}img{max-width:100%}footer{margin-top:2.5em;padding-top:0.8em;border-top:2px solid #231640;font-size:0.85em;color:#444}.note{color:#555;font-size:0.9em}",
                ".banner{display:flex;align-items:center;gap:1.2em;padding:0.4em 0 0.9em 0;color:#231640;",
                "border-bottom:2px solid #231640;margin-bottom:1.2em}",
                ".banner img{height:120px;width:auto}.banner h1{margin:0;font-size:1.8em;letter-spacing:0.02em}",
                ".banner p{margin:0.25em 0 0 0;color:#444}table.glance td:first-child{font-weight:bold;width:7em}</style>")
  lib_sum <- if (nrow(libc)) {
    s <- do.call(rbind, lapply(split(libc, libc$locus), function(x) data.frame(locus = x$locus[1], species = length(unique(x$species)),
                   sequences = nrow(x), from_phylotaR = sum(x$source == "phylotaR"), from_genbank_plastomes = sum(x$source == "genbank_plastome"))))
    s[order(s$locus, method = "radix"), , drop = FALSE]
  } else NULL
  paths <- character(0)
  for (q in unique(tab$unit)) {
    t <- tab[tab$unit == q, , drop = FALSE]
    t <- t[order(t$compartment != "plastid", t$locus, method = "radix"), , drop = FALSE]
    named_genera <- sort(unique(t$genus_idtaxa[!is.na(t$genus_idtaxa)]), method = "radix")
    named_species <- unique(c(t$species_idtaxa[!is.na(t$species_idtaxa)], alternatives$species[alternatives$query %in% t$query]))
    specimen <- if (any(rec$kind == "sample")) {
      dec <- get("declared_species"); vou <- get("voucher")
      paste0("declared species: ", if (is.na(dec)) "not declared" else dec,
             "; voucher: ", if (is.na(vou)) "none (no voucher: the leakage rule per specimen cannot be applied to this sample)" else vou,
             "; answer against the declared species: ", .bc_declared_comparison(t, dec))
    } else NULL
    idl <- t[!is.na(t$identical_to_library) & !t$reason %in% c("no_overlap", "assembly_failed", "single_species_library"), , drop = FALSE]
    identical <- if (nrow(idl)) paste0(paste(sprintf("%s identical to library accession %s", idl$locus, idl$identical_to_library),
                                             collapse = "; "),
                                       ": the answer of ", if (nrow(idl) > 1L) "these loci" else "this locus",
                                       " is the query finding itself in the library, not an independent identification") else NULL
    ran <- paste0(if (is.na(get("date"))) "date not recorded" else get("date"), " on ",
                  if (is.na(get("machine"))) "machine not recorded" else
                    paste0(get("machine"), " (", get("operating_system"), "; ", get("platform"), ")"))
    glance <- .bc_html_table(data.frame(
      item = c("Identification run", "Data", if (!is.null(specimen)) "Specimen", "Loci found", if (!is.null(identical)) "Identical to the library",
               "Library", "Answers", "Reading"),
      value = c(ran, .bc_report_data_kind(input_type, t$query_length, t$query), specimen, .bc_report_loci_found(t, get("loci_declared")),
                identical,
                .bc_report_library_words(libc),
                .bc_report_answers_words(t),
                "Each locus is answered on its own; no call is made across loci (section 3)")), class = "glance")
    sec <- function(i, title, body) sprintf("<h2 id=\"section-%d\">%d. %s</h2>\n%s", i, i, title, body)
    s1 <- paste0(.bc_html_table(data.frame(item = c("Query", "Input", "Input type", "md5 of the input", "Run"),
                                           value = c(q, get("path"), .bc_report_data_kind(input_type, t$query_length, t$query), get("md5"), run_name))),
                 "<p>Steps of the route: ", if (is.na(get("steps"))) "not recorded" else .bc_html_escape(get("steps")), ". Steps in grey were not needed for this query.</p>",
                 .bc_html_plot(.bc_report_route_plot(input_type, t), width = 9, height = 1.8))
    # The loci found only; the loci absent from the query are named once (every row stays in the CSV)
    fnd <- !t$reason %in% c("no_overlap", "assembly_failed", "single_species_library")
    tf <- t[fnd, , drop = FALSE]
    nf <- setdiff(sort(unique(t$locus), method = "radix"), tf$locus)
    af <- answer[paste(answer$query, answer$locus) %in% paste(tf$query, tf$locus), if (per_run) names(answer) else setdiff(names(answer), "query"),
                 drop = FALSE]
    s2 <- paste0("<p>One answer per locus, each from its own model; state 1 names a species, state 2 a genus, state 3 ",
                 "is not assignable. The dashed line is the threshold of ", threshold, ".</p>",
                 if (nrow(tf)) .bc_html_plot(.bc_report_answer_plot(tf, threshold)) else "<p>No locus of the library was found.</p>",
                 .bc_html_table(af),
                 "<p>Loci not found in the ", if (per_run) "sample" else "query", ": ",
                 if (length(nf)) .bc_html_escape(paste(nf, collapse = ", ")) else "none", ".</p>")
    ag <- agreement[agreement$query == q, setdiff(names(agreement), "query"), drop = FALSE]
    s3 <- paste0("<p>The loci are read side by side; no call is made across them (rule E8 of the validation plan). ",
                 "For each genus and species named by any locus, at any confidence: the loci that name it and how many of them ",
                 "assign it at the threshold. Loci with no answer are listed under rank <em>none</em>.</p>", .bc_html_table(ag))
    al <- alternatives[alternatives$query %in% t$query, if (per_run) names(alternatives) else setdiff(names(alternatives), "query"),
                       drop = FALSE]
    s4 <- paste0("<p>Alternatives under the threshold of ", threshold, ", not assigned. IdTaxa gives one path per query and ",
                 "locus: the alternatives are that path when it falls under the threshold, and in state 2 the species of the ",
                 "named genus in the library. They are not identifications; they say where the data point when they do not ",
                 "decide.</p>",
                 "<p>Best path under the threshold, per locus:</p>",
                 .bc_html_table(al[al$kind == "best_path_under_threshold", setdiff(names(al), c("kind", "assigned")), drop = FALSE]),
                 "<p>Species of the named genus in the library, for the loci in state 2 (all of them in ",
                 "<code>REPORT_", .bc_html_escape(run_name), "_alternatives.csv</code>):</p>",
                 .bc_html_table(.bc_report_candidates_summary(al)))
    rates <- unique(t[, c("locus", "validation_ws_rate_species_present", "validation_ws_rate_species_absent")])
    gd_q <- if (is.null(gd)) NULL else gd[gd$genus %in% named_genera, , drop = FALSE]
    s5 <- paste0("<p>Wrong-species rate of each locus in the validation (6D): with the species of the query in the library ",
                 "(scheme E) and without it (scheme G).</p>", .bc_html_table(rates),
                 if (is.null(gd)) "<p>Per-genus discrimination not available: run summarise_barcoding_genus_discrimination().</p>"
                 else paste0("<p>Per-genus discrimination of each locus in the validation folds, for the genera named above.</p>",
                             .bc_html_table(gd_q)))
    sq <- sampling[sampling$genus %in% named_genera, , drop = FALSE]
    s6 <- if (nrow(sq)) paste0("<p>Sequences per species and locus in the library for the genera named by the query; empty cells ",
                               "are species with no sequence of that locus; species named by the query in bold.</p>",
                               .bc_html_plot(.bc_report_sampling_plot(sq, named_species), width = 8,
                                             height = 1 + 0.25 * length(unique(sq$species)) + 0.3 * length(unique(sq$genus))))
          else "<p>No genus was named by any locus.</p>"
    s7 <- if (identical(input_type, "reads")) .bc_html_table(rec[rec$kind %in% c("output", "tool", "command"), , drop = FALSE])
          else "<p>Not applicable: the query is not a set of reads.</p>"
    s8 <- paste0("<p>Names follow the Cactaceae checklist of Korotkova et al. (2021, Willdenowia 51(2): 251-270). ",
                 "Sources of the library: phylotaR clusters of GenBank records and, where present, loci cut from GenBank plastomes.</p>",
                 .bc_html_table(lib_sum), .bc_html_table(rec[rec$kind %in% c("library", "model", "software", "setting"), , drop = FALSE]))
    s9 <- paste0("<p>Standard plant barcodes resolve about half of the species, complete plastomes a little more, and a part of ",
                 "species is not resolvable by any of them, because of hybridisation, incomplete lineage sorting and species that ",
                 "are evolutionary grades rather than clades (Zeng et al., 2026, New Phytologist). In Cactaceae the plastid loci ",
                 "vary little within genera and many species have one sequence or none in the library (section 6). A state 3 is ",
                 "not evidence that the species is absent from the library; a state 1 relies on the species being in it. If the ",
                 "query comes from a specimen already in the library, its answer is not an independent test.</p>")
    html <- paste0("<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\"/><title>PhyloCactus identification report \U0001F335: ",
                   .bc_html_escape(q), "</title>", css, "</head><body>\n",
                   .bc_report_banner(paste0(if (per_run) "Sample " else "Query ", q, " \u00b7 run ", run_name)), "\n", glance, "\n",
                   sec(1, "Query and route", s1), "\n", sec(2, "Answer per locus", s2), "\n",
                   sec(3, "Reading across loci", s3), "\n", sec(4, "Alternatives (not assigned)", s4), "\n",
                   sec(5, "How far to trust each locus", s5), "\n", sec(6, "Sampling of the named genera in the library", s6), "\n",
                   sec(7, "Reads", s7), "\n", sec(8, "Library, models and software", s8), "\n", sec(9, "Limits", s9),
                   "\n", .bc_report_footer(get("date")), "\n</body></html>\n")
    p <- file.path(output_dir, if (per_run) paste0("REPORT_", run_name, ".html")
                               else paste0("REPORT_", run_name, "_", gsub("[^A-Za-z0-9._-]", "_", q), ".html"))
    writeLines(html, p, useBytes = TRUE)
    paths <- c(paths, p)
  }
  # Checksums of what the report rests on and of what it wrote (the HTML cannot hold its own md5)
  ck_files <- c(f_tab, f_rec, paths, paste0(base, c("answer", "alternatives", "agreement", "sampling", "provenance"), ".csv"))
  utils::write.csv(data.frame(file = basename(ck_files), md5 = unname(tools::md5sum(ck_files)), stringsAsFactors = FALSE),
                   paste0(base, "checksums.csv"), row.names = FALSE)
  message("Identification report: ", paste(paths, collapse = ", "))
  invisible(paths)
}
