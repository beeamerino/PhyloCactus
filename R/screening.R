#' Screen Locus Alignments for Substitution Saturation and Informativeness
#'
#' Evaluates phylogenetic informativeness, sequence coverage, alignment length, and substitution saturation across individual locus alignments.
#' Filtering out loci exhibiting high substitution saturation or severe site length anomalies prevents systematic noise
#' and long-branch attraction (LBA) artifacts from distorting maximum-likelihood supermatrix inference.
#'
#' @param fasta_folder Character. Directory path containing aligned locus FASTA files.
#' @param out_base Character. Base destination directory for diagnostic plots and screened FASTA outputs (`filtered_markers/`).
#' @param min_cols_to_evaluate Integer. Minimum number of alignment columns required to compute saturation metrics. Defaults to `50L`.
#' @param min_aln_len_to_retain Integer. Minimum alignment length in base pairs required to retain a locus. Defaults to `200L`.
#' @param min_nseq_to_retain Integer. Minimum number of sequences required per locus alignment. Defaults to `100L`.
#' @param max_marker_missing Numeric. Maximum allowable missing data fraction per locus. Defaults to `0.7`.
#' @param saturation_flag_cutoff Numeric. Uncorrected p-distance vs. raw distance slope threshold to flag substitution saturation. Defaults to `0.3`.
#' @param saturation_keep_cutoff Numeric. Saturation slope cutoff threshold below which saturated loci are excluded. Defaults to `0.5`.
#' @param iqr_multiplier Numeric. Interquartile range (IQR) multiplier for identifying site-length outlier bounds. Defaults to `1.5`.
#' @return A data frame summarizing saturation statistics, alignment dimensions, and retention decisions across screened loci.
#' @examples
#' \dontrun{
#' run_marker_screening(
#'   fasta_folder = "2_MAFFT_Cactaceae/alignments",
#'   out_base = "3_Screening_Ingroup",
#'   min_aln_len_to_retain = 200L,
#'   min_nseq_to_retain = 50L
#' )
#' }
#' @export
run_marker_screening <- function(
  fasta_folder,
  out_base,
  min_cols_to_evaluate = 50L,
  min_aln_len_to_retain = 200L,
  min_nseq_to_retain = 100L,
  max_marker_missing = 0.7,
  saturation_flag_cutoff = 0.3,
  saturation_keep_cutoff = 0.5,
  iqr_multiplier = 1.5
) {
  dir_markers  <- file.path(out_base, "filtered_markers")
  dir.create(out_base, showWarnings = FALSE, recursive = TRUE)
  dir.create(dir_markers, showWarnings = FALSE, recursive = TRUE)
  
  out_pdf             <- file.path(out_base, "SUPP_FIG_marker_screening_diagnostics.pdf")
  out_summary_csv     <- file.path(out_base, "SUPP_TABLE_marker_screening_summary.csv")
  out_log             <- file.path(out_base, "LOG_marker_screening.txt")
  out_sessioninfo_txt <- file.path(out_base, "LOG_marker_screening_sessionInfo.txt")
  
  log_message <- function(...) {
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
    cat(msg, "\n")
    cat(msg, "\n", file = out_log, append = TRUE)
  }
  
  if (file.exists(out_log)) file.remove(out_log)
  
  log_message("Starting ingroup marker screening.")
  log_message("Input directory: ", fasta_folder)
  
  fasta_files <- list.files(
    fasta_folder,
    pattern = "^ALN_masked_final_.*\\.fasta$",
    full.names = TRUE
  )
  fasta_files <- sort(fasta_files)
  
  if (length(fasta_files) == 0) {
    stop("No input FASTA files found in: ", fasta_folder, call. = FALSE)
  }
  
  safe_as_matrix <- function(dna_bin) {
    mat <- as.character(dna_bin)
    if (is.null(dim(mat))) mat <- matrix(mat, nrow = 1)
    rownames(mat) <- rownames(dna_bin)
    mat
  }
  
  plot_pairwise_distances <- function(d, marker, title_suffix = "") {
    df <- data.frame(dist = as.numeric(d), index = seq_along(as.numeric(d)))
    ggplot2::ggplot(df, ggplot2::aes(index, dist)) +
      ggplot2::geom_point(size = 1.2) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::ggtitle(paste(marker, title_suffix)) +
      ggplot2::xlab("Pair index") + ggplot2::ylab("ML distance") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"), axis.title = ggplot2::element_text(size = 9))
  }
  
  plot_saturation_proxy <- function(dna_bin, marker, title_suffix = "") {
    dist_k80 <- ape::dist.dna(dna_bin, model = "K80", pairwise.deletion = TRUE)
    dist_k80_gamma <- ape::dist.dna(dna_bin, model = "K80", pairwise.deletion = TRUE, gamma = TRUE)
    
    df <- data.frame(k80 = as.numeric(dist_k80), k80_gamma = as.numeric(dist_k80_gamma))
    
    ggplot2::ggplot(df, ggplot2::aes(k80, k80_gamma)) +
      ggplot2::geom_point(size = 1.2) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::ggtitle(paste(marker, title_suffix)) +
      ggplot2::xlab("K80 pairwise distance") + ggplot2::ylab("K80 pairwise distance (gamma)") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"), axis.title = ggplot2::element_text(size = 9))
  }
  
  detect_outlier_sequences <- function(d, seq_names, iqr_multiplier = 1.5) {
    dm <- as.matrix(d)
    if (nrow(dm) < 2 || ncol(dm) < 2) return(list(outliers = character(0), upper_cut = NA_real_, iqr = NA_real_, mean_dist = rep(NA_real_, length(seq_names))))
    mean_dist <- rowMeans(dm, na.rm = TRUE)
    if (all(is.na(mean_dist))) return(list(outliers = character(0), upper_cut = NA_real_, iqr = NA_real_, mean_dist = mean_dist))
    q1 <- as.numeric(stats::quantile(mean_dist, 0.25, na.rm = TRUE))
    q3 <- as.numeric(stats::quantile(mean_dist, 0.75, na.rm = TRUE))
    iqr_value <- q3 - q1
    upper_cut <- q3 + iqr_multiplier * iqr_value
    list(outliers = seq_names[mean_dist > upper_cut], upper_cut = upper_cut, iqr = iqr_value, mean_dist = mean_dist)
  }
  
  extract_marker_name <- function(path) sub("^ALN_masked_final_", "", tools::file_path_sans_ext(basename(path)))
  
  test_saturation_proxy <- function(dna_bin, saturation_flag_cutoff = 0.3) {
    dist_1 <- tryCatch(as.numeric(ape::dist.dna(dna_bin, model = "K80", pairwise.deletion = TRUE)), error = function(e) NA_real_)
    dist_2 <- tryCatch(as.numeric(ape::dist.dna(dna_bin, model = "K80", pairwise.deletion = TRUE, gamma = TRUE)), error = function(e) NA_real_)
    ok <- !(is.na(dist_1) | is.na(dist_2))
    dist_1 <- dist_1[ok]; dist_2 <- dist_2[ok]
    if (length(dist_1) < 5) return(list(slope = NA_real_, saturated = NA, reason = "insufficient_data"))
    df <- data.frame(x = dist_1, y = dist_2)
    if (is.na(stats::var(df$x)) || is.na(stats::var(df$y))) return(list(slope = NA_real_, saturated = NA, reason = "variance_na"))
    fit <- tryCatch(stats::lm(y ~ x, data = df), error = function(e) NULL)
    if (is.null(fit)) return(list(slope = NA_real_, saturated = NA, reason = "regression_failed"))
    slope <- unname(stats::coef(fit)[2])
    list(slope = slope, saturated = isTRUE(slope < saturation_flag_cutoff), reason = "ok")
  }
  
  alignment_stats <- function(dna_bin) {
    aln_mat <- toupper(safe_as_matrix(dna_bin))
    gaps_mask_orig <- aln_mat == "-"
    nseq <- nrow(aln_mat); total_cols <- ncol(aln_mat)
    aln_base_only <- aln_mat
    aln_base_only[!aln_base_only %in% c("A", "C", "G", "T")] <- NA
    informative_col <- apply(aln_base_only, 2, function(col) any(!is.na(col)))
    aln_len_effective <- sum(informative_col)
    missing_prop <- sum(is.na(aln_base_only)) / (nseq * total_cols)
    gap_prop <- sum(gaps_mask_orig, na.rm = TRUE) / (nseq * total_cols)
    n_bases <- sum(aln_base_only %in% c("A", "C", "G", "T"), na.rm = TRUE)
    gc_prop <- if (n_bases == 0) NA_real_ else sum(aln_base_only == "G" | aln_base_only == "C", na.rm = TRUE) / n_bases
    var_sites <- sum(apply(aln_base_only, 2, function(col) length(unique(stats::na.omit(col))) > 1))
    pi_sites <- sum(apply(aln_base_only, 2, function(col) { tt <- table(stats::na.omit(col)); any(tt >= 2) && length(tt) >= 2 }))
    ent_vals <- apply(aln_base_only, 2, function(col) {
      tt <- table(stats::na.omit(col))
      if (length(tt) == 0) return(NA_real_)
      p <- tt / sum(tt); -sum(p * log(p))
    })
    list(aln_len = aln_len_effective, total_cols = total_cols, missing = missing_prop, gaps = gap_prop, gc = gc_prop, var_sites = var_sites, pi_sites = pi_sites, entropy = mean(ent_vals, na.rm = TRUE))
  }
  
  make_summary_row <- function(marker, status, decision, decision_reason, n_original = NA_integer_, n_filtered = NA_integer_, n_removed = NA_integer_, aln_len_before = NA_real_, aln_len_after = NA_real_, missing_before = NA_real_, missing_after = NA_real_, gc_before = NA_real_, gc_after = NA_real_, var_sites_before = NA_real_, var_sites_after = NA_real_, pi_sites_before = NA_real_, pi_sites_after = NA_real_, entropy_before = NA_real_, entropy_after = NA_real_, slope_before = NA_real_, slope_after = NA_real_, saturated_before = NA, saturated_after = NA, reason_before = NA_character_, reason_after = NA_character_, removed = "-") {
    data.frame(marker = marker, status = status, decision = decision, decision_reason = decision_reason, n_original = n_original, n_filtered = n_filtered, n_removed = n_removed, aln_len_before = aln_len_before, aln_len_after = aln_len_after, missing_before = missing_before, missing_after = missing_after, gc_before = gc_before, gc_after = gc_after, var_sites_before = var_sites_before, var_sites_after = var_sites_after, pi_sites_before = pi_sites_before, pi_sites_after = pi_sites_after, entropy_before = entropy_before, entropy_after = entropy_after, slope_before = slope_before, slope_after = slope_after, saturated_before = saturated_before, saturated_after = saturated_after, reason_before = reason_before, reason_after = reason_after, removed = removed, stringsAsFactors = FALSE)
  }
  
  results <- list(); filtered_alignments <- list(); pw_before_list <- list(); sat_before_list <- list(); pw_after_list <- list(); sat_after_list <- list()
  
  for (f in fasta_files) {
    marker <- extract_marker_name(f)
    log_message("Processing marker: ", marker)
    res <- tryCatch({
      aln <- ape::read.dna(f, format = "fasta")
      aln_char <- toupper(safe_as_matrix(aln))
      aln_char[!(aln_char %in% c("A", "C", "G", "T"))] <- "-"
      aln <- ape::as.DNAbin(aln_char)
      seq_names <- rownames(aln)
      if (is.null(seq_names) || length(seq_names) != nrow(aln)) seq_names <- paste0("seq_", seq_len(nrow(aln))); rownames(aln) <- seq_names
      if (ncol(aln) < min_cols_to_evaluate) {
        log_message("Marker skipped for short alignment (<", min_cols_to_evaluate, " columns): ", marker)
        return(list(summary_row = make_summary_row(marker = marker, status = "skipped_short_alignment", decision = "NO", decision_reason = paste0("alignment_shorter_than_", min_cols_to_evaluate, "_columns"), n_original = nrow(aln), n_filtered = nrow(aln), n_removed = 0L, aln_len_before = ncol(aln), aln_len_after = ncol(aln), removed = "-"), pw_before = NULL, sat_before = NULL, pw_after = NULL, sat_after = NULL, filtered_alignment = NULL))
      }
      stats_before <- alignment_stats(aln)
      d_ml <- tryCatch(phangorn::dist.ml(aln), error = function(e) stop("dist.ml failed before filtering: ", e$message, call. = FALSE))
      sat_before <- test_saturation_proxy(aln, saturation_flag_cutoff = saturation_flag_cutoff)
      pw_before <- tryCatch(plot_pairwise_distances(d_ml, marker), error = function(e) NULL)
      sat_before_plot <- tryCatch(plot_saturation_proxy(aln, marker), error = function(e) NULL)
      outinfo <- detect_outlier_sequences(d_ml, seq_names, iqr_multiplier = iqr_multiplier)
      outliers <- outinfo$outliers
      aln_filtered <- aln
      if (length(outliers) > 0) aln_filtered <- aln[!(seq_names %in% outliers), , drop = FALSE]
      stats_after <- alignment_stats(aln_filtered)
      d_ml_filt <- tryCatch(phangorn::dist.ml(aln_filtered), error = function(e) stop("dist.ml failed after filtering: ", e$message, call. = FALSE))
      sat_after <- test_saturation_proxy(aln_filtered, saturation_flag_cutoff = saturation_flag_cutoff)
      pw_after <- tryCatch(plot_pairwise_distances(d_ml_filt, marker), error = function(e) NULL)
      sat_after_plot <- tryCatch(plot_saturation_proxy(aln_filtered, marker), error = function(e) NULL)
      keep_marker <- (!isTRUE(sat_after$saturated) && !is.na(sat_after$slope) && sat_after$slope > saturation_keep_cutoff && stats_after$missing < max_marker_missing && stats_after$var_sites > 0 && stats_after$aln_len >= min_aln_len_to_retain && nrow(aln_filtered) >= min_nseq_to_retain)
      decision_reason <- if (keep_marker) "passed_all_thresholds" else paste(c(if (isTRUE(sat_after$saturated)) "flagged_as_saturated" else NULL, if (is.na(sat_after$slope)) "slope_na" else NULL, if (!is.na(sat_after$slope) && sat_after$slope <= saturation_keep_cutoff) paste0("slope_le_", saturation_keep_cutoff) else NULL, if (stats_after$missing >= max_marker_missing) paste0("missing_ge_", max_marker_missing) else NULL, if (stats_after$var_sites <= 0) "no_variable_sites" else NULL, if (stats_after$aln_len < min_aln_len_to_retain) paste0("aln_len_lt_", min_aln_len_to_retain) else NULL, if (nrow(aln_filtered) < min_nseq_to_retain) paste0("nseq_lt_", min_nseq_to_retain) else NULL), collapse = ";")
      summary_row <- make_summary_row(marker = marker, status = if (keep_marker) "processed_retained" else "processed_rejected", decision = if (keep_marker) "YES" else "NO", decision_reason = decision_reason, n_original = nrow(aln), n_filtered = nrow(aln_filtered), n_removed = length(outliers), aln_len_before = stats_before$aln_len, aln_len_after = stats_after$aln_len, missing_before = stats_before$missing, missing_after = stats_after$missing, gc_before = stats_before$gc, gc_after = stats_after$gc, var_sites_before = stats_before$var_sites, var_sites_after = stats_after$var_sites, pi_sites_before = stats_before$pi_sites, pi_sites_after = stats_after$pi_sites, entropy_before = stats_before$entropy, entropy_after = stats_after$entropy, slope_before = sat_before$slope, slope_after = sat_after$slope, saturated_before = sat_before$saturated, saturated_after = sat_after$saturated, reason_before = sat_before$reason, reason_after = sat_after$reason, removed = if (length(outliers) == 0) "-" else paste(outliers, collapse = ";"))
      list(summary_row = summary_row, pw_before = pw_before, sat_before = sat_before_plot, pw_after = pw_after, sat_after = sat_after_plot, filtered_alignment = if (keep_marker) aln_filtered else NULL)
    }, error = function(e) {
      log_message("Processing failed for marker ", marker, ": ", e$message)
      list(summary_row = make_summary_row(marker = marker, status = "failed_processing", decision = "NO", decision_reason = paste0("processing_error:", e$message)), pw_before = NULL, sat_before = NULL, pw_after = NULL, sat_after = NULL, filtered_alignment = NULL)
    })
    results[[marker]] <- res$summary_row
    pw_before_list[[marker]] <- res$pw_before; sat_before_list[[marker]] <- res$sat_before
    pw_after_list[[marker]] <- res$pw_after; sat_after_list[[marker]] <- res$sat_after
    if (!is.null(res$filtered_alignment)) filtered_alignments[[marker]] <- res$filtered_alignment
  }
  
  summary_table <- dplyr::bind_rows(results)
  if (nrow(summary_table) == 0) stop("Summary table is empty; no markers were processed.", call. = FALSE)
  summary_table <- dplyr::arrange(summary_table, marker)
  write.csv(summary_table, out_summary_csv, row.names = FALSE)
  
  grDevices::cairo_pdf(filename = out_pdf, width = 11, height = 8.5, pointsize = 12)
  arrange_plots_on_page <- function(plot_list, ncol = 2, top_title = "") {
    plots <- unname(plot_list[!vapply(plot_list, is.null, logical(1))])
    if (length(plots) == 0) {
      grid::grid.newpage(); grid::grid.text(paste("No plots available:", top_title), gp = grid::gpar(fontsize = 16, fontface = "bold"))
      return(invisible(NULL))
    }
    print(patchwork::wrap_plots(plots, ncol = ncol) + patchwork::plot_annotation(title = top_title, theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5))))
  }
  arrange_plots_on_page(pw_before_list, ncol = 3, top_title = "Pairwise ML distances (before outlier filtering)")
  arrange_plots_on_page(pw_after_list, ncol = 3, top_title = "Pairwise ML distances (after outlier filtering)")
  arrange_plots_on_page(sat_before_list, ncol = 3, top_title = "Saturation-style proxy plots (before outlier filtering)")
  arrange_plots_on_page(sat_after_list, ncol = 3, top_title = "Saturation-style proxy plots (after outlier filtering)")
  grDevices::dev.off()
  
  selected_markers <- dplyr::pull(dplyr::filter(summary_table, decision == "YES"), marker)
  for (marker in names(filtered_alignments)) {
    out_fasta <- file.path(dir_markers, paste0(marker, ".fasta"))
    ape::write.dna(filtered_alignments[[marker]], file = out_fasta, format = "fasta", nbcol = -1, colw = 9999)
  }
  
  saveRDS(list(pw_before = pw_before_list, pw_after = pw_after_list, sat_before = sat_before_list, sat_after = sat_after_list), file = file.path(out_base, "RDS_marker_screening_plots.rds"))
  writeLines(capture.output(utils::sessionInfo()), con = out_sessioninfo_txt)
  
  message("\nMarker screening and quality check completed. ")
  return(summary_table)
}

#' Final Marker Integration and Taxonomic Cleaning
#'
#' @param ingroup_dir Character.
#' @param outgroup_dir Character.
#' @param output_dir Character.
#' @param accepted_list_file Character.
#' @param metadata_in_file Character.
#' @param metadata_out_file Character.
#' @return NULL
#' @export
integrate_and_clean_markers <- function(
  ingroup_dir,
  outgroup_dir,
  output_dir,
  accepted_list_file,
  metadata_in_file,
  metadata_out_file
) {
  
  out_marker_dir <- file.path(output_dir, "cleaned_markers")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_marker_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Starting marker integration and cleaning pipeline... ")
  message(sprintf("Input ingroup directory: %s", ingroup_dir))
  message(sprintf("Input outgroup directory: %s", outgroup_dir))
  message(sprintf("Output directory: %s", output_dir))
  
  assert_exists <- function(path, type = c("file", "dir")) {
    type <- match.arg(type)
    ok <- if (type == "file") file.exists(path) else dir.exists(path)
    if (!ok) stop(sprintf("Missing %s: %s", type, path), call. = FALSE)
  }

  assert_cols <- function(df, required_cols, object_name) {
    missing_cols <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0) {
      stop(sprintf("%s is missing required columns: %s", object_name, paste(missing_cols, collapse = ", ")), call. = FALSE)
    }
  }

  strip_marker_prefix <- function(x) {
    x |> basename() |> stringr::str_remove("\\.fasta$") |> stringr::str_remove("^ALN_clean_input_FASTA_MARKER_") |> stringr::str_remove("^ALN_masked_final_") |> stringr::str_remove("^MARKER_") |> stringr::str_remove("^Masked_")
  }

  normalize_marker <- function(x) {
    x |> strip_marker_prefix() |> stringr::str_trim() |> stringr::str_replace_all("[[:space:]/\\\\\\.-]+", "_") |> stringr::str_replace_all("^_+|_+$", "")
  }

  normalize_species <- function(x) {
    x |> stringr::str_trim() |> stringr::str_replace_all("[ -]+", "_")
  }

  extract_species <- function(header) {
    if (is.na(header) || !nzchar(header)) return(NA_character_)
    first_field <- stringr::str_split(header, "\\|", simplify = TRUE)[1]
    first_field <- stringr::str_trim(first_field)
    if (!nzchar(first_field)) return(NA_character_)
    normalize_species(first_field)
  }

  read_fasta_safe <- function(paths) {
    if (length(paths) == 0) return(Biostrings::DNAStringSet())
    fasta_list <- lapply(paths, function(p) {
      tryCatch(Biostrings::readDNAStringSet(p), error = function(e) { warning(sprintf("Could not read FASTA: %s -> %s", p, e$message), call. = FALSE); Biostrings::DNAStringSet() })
    })
    fasta_list <- fasta_list[lengths(fasta_list) > 0]
    if (length(fasta_list) == 0) return(Biostrings::DNAStringSet())
    do.call(c, fasta_list)
  }

  write_csv_dual <- function(df, legacy_name, standard_name) {
    legacy_path <- file.path(output_dir, legacy_name)
    standard_path <- file.path(output_dir, standard_name)
    readr::write_csv(df, legacy_path, na = "")
    if (!identical(legacy_name, standard_name)) readr::write_csv(df, standard_path, na = "")
  }

  collapse_unique <- function(x, sep = ";") {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    x <- sort(unique(x))
    if (length(x) == 0) return(NA_character_)
    paste(x, collapse = sep)
  }

  assert_exists(ingroup_dir, "dir")
  assert_exists(outgroup_dir, "dir")
  assert_exists(accepted_list_file, "file")
  assert_exists(metadata_in_file, "file")
  assert_exists(metadata_out_file, "file")

  message("Indexing FASTA files...")
  ingroup_files <- sort(list.files(ingroup_dir, full.names = TRUE, pattern = "\\.fasta$", ignore.case = TRUE))
  outgroup_files <- sort(list.files(outgroup_dir, full.names = TRUE, pattern = "^ALN_masked_final_.*\\.fasta$", ignore.case = TRUE))
  
  if (length(ingroup_files) == 0) stop(sprintf("No FASTA files found in ingroup_dir: %s", ingroup_dir), call. = FALSE)
  
  ingroup_index <- dplyr::tibble(path = ingroup_files, basename = basename(ingroup_files), marker_key = normalize_marker(basename))
  outgroup_index <- dplyr::tibble(path = outgroup_files, basename = basename(outgroup_files), marker_key = normalize_marker(basename))
  markers_to_process <- sort(unique(ingroup_index$marker_key))

  message("Reading and validating metadata...")
  meta_in <- read.csv(metadata_in_file, stringsAsFactors = FALSE)
  meta_out <- read.csv(metadata_out_file, stringsAsFactors = FALSE)

  assert_cols(meta_in, c("species", "Marker_std", "sid"), "meta_in")
  assert_cols(meta_out, c("species", "Marker_std", "sid"), "meta_out")

  metadata_all <- dplyr::bind_rows(
    dplyr::mutate(meta_in, is_outgroup = FALSE),
    dplyr::mutate(meta_out, is_outgroup = TRUE)
  ) |>
    dplyr::mutate(species = normalize_species(species), marker_key = normalize_marker(as.character(Marker_std)), sid = as.character(sid)) |>
    dplyr::distinct(species, marker_key, sid, .keep_all = TRUE)

  message("Building raw integration table...")
  tabla_raw_list <- vector("list", length(markers_to_process))

  for (i in seq_along(markers_to_process)) {
    mk <- markers_to_process[i]
    message(sprintf("Reading marker [%d/%d]: %s", i, length(markers_to_process), mk))
    
    ingroup_paths <- dplyr::pull(dplyr::filter(ingroup_index, marker_key == mk), path)
    outgroup_paths <- dplyr::pull(dplyr::filter(outgroup_index, marker_key == mk), path)
    
    ingroup_seq <- read_fasta_safe(ingroup_paths)
    outgroup_seq <- read_fasta_safe(outgroup_paths)
    
    df_in <- dplyr::tibble(marker_source = mk, marker_key = mk, source_branch = "ingroup", species = vapply(names(ingroup_seq), extract_species, character(1)), fasta_name = names(ingroup_seq))
    df_out <- dplyr::tibble(marker_source = mk, marker_key = mk, source_branch = "outgroup", species = vapply(names(outgroup_seq), extract_species, character(1)), fasta_name = names(outgroup_seq))
    
    combined_df <- dplyr::bind_rows(df_in, df_out)
    if (nrow(combined_df) == 0) {
      warning(sprintf("Marker %s contains no sequences. Skipping.", mk), call. = FALSE)
      next
    }
    tabla_raw_list[[i]] <- combined_df
  }

  tabla_raw <- dplyr::bind_rows(tabla_raw_list)
  if (nrow(tabla_raw) == 0) stop("No sequences were recovered from ingroup/outgroup FASTA files.", call. = FALSE)

  tabla_raw <- dplyr::left_join(tabla_raw, metadata_all, by = c("species", "marker_key"))
  unmatched_metadata <- dplyr::distinct(dplyr::filter(tabla_raw, is.na(sid)), marker_key, source_branch, species, fasta_name)

  message("Exporting TABLE_raw_species_metadata.csv ...")
  write_csv_dual(tabla_raw, legacy_name = "Tabla_raw_species_metadata.csv", standard_name = "TABLE_raw_species_metadata.csv")
  if (nrow(unmatched_metadata) > 0) write_csv_dual(unmatched_metadata, legacy_name = "unmatched_metadata_rows.csv", standard_name = "TABLE_unmatched_metadata_rows.csv")

  message("Reading checklist and defining accepted taxa...")
  full_checklist <- if (grepl("\\.xlsx?$", accepted_list_file, ignore.case = TRUE)) readxl::read_excel(accepted_list_file) else read.csv(accepted_list_file, stringsAsFactors = FALSE)

  if (!"pureName" %in% names(full_checklist)) {
    nm <- tolower(names(full_checklist))
    candidates <- intersect(nm, c("scientificname", "species", "name", "taxon", "purename"))
    if (length(candidates) == 0) stop("Checklist missing a scientific name column.", call. = FALSE)
    names(full_checklist)[nm == candidates[1]] <- "pureName"
  }

  accepted_species <- full_checklist |>
    dplyr::mutate(pureName = stringr::str_squish(as.character(pureName))) |>
    dplyr::filter(!stringr::str_detect(pureName, stringr::regex("\\b(subsp|ssp|var|forma|f\\.|subg|sect|ser|cf\\.|aff\\.|sp\\.|spp\\.|nr\\.)\\b|x|\\bx\\b", ignore_case = TRUE))) |>
    dplyr::filter(stringr::str_detect(pureName, "^[A-Z][a-z]+\\s[a-z]+$")) |>
    dplyr::mutate(species = normalize_species(pureName)) |>
    dplyr::pull(species) |> unique()

  outgroup_species <- unique(dplyr::pull(dplyr::filter(metadata_all, is_outgroup), species))

  tabla_raw <- tabla_raw |>
    dplyr::mutate(
      matched_to_metadata = !is.na(sid),
      accepted_in_checklist = species %in% accepted_species,
      is_outgroup_species = species %in% outgroup_species,
      accepted_or_outgroup = accepted_in_checklist | is_outgroup_species,
      matched_to_checklist = dplyr::if_else(source_branch == "ingroup", accepted_in_checklist, NA),
      record_status = dplyr::case_when(
        !matched_to_metadata                       ~ "unmatched_metadata",
        source_branch == "outgroup"                ~ "accepted_outgroup",
        accepted_in_checklist                      ~ "accepted_ingroup",
        source_branch == "ingroup"                 ~ "rejected_ingroup",
        TRUE                                       ~ "other"
      )
    )

  message("Performing duplicate audit and master registry...")
  duplicate_groups <- tabla_raw |>
    dplyr::filter(!is.na(species), !is.na(marker_key)) |>
    dplyr::group_by(marker_key, species) |>
    dplyr::summarise(
      n_records_total    = dplyr::n(),
      n_records_ingroup  = sum(source_branch == "ingroup", na.rm = TRUE),
      n_records_outgroup = sum(source_branch == "outgroup", na.rm = TRUE),
      n_records_accepted = sum(accepted_or_outgroup %in% TRUE, na.rm = TRUE),
      n_records_rejected = sum(source_branch == "ingroup" & accepted_or_outgroup %in% FALSE, na.rm = TRUE),
      n_classes          = dplyr::n_distinct(accepted_or_outgroup),
      classification_conflict = n_classes > 1,
      duplicate_class = dplyr::case_when(
        dplyr::n() <= 1                            ~ "unique",
        n_records_ingroup > 0 & n_records_outgroup > 0 ~ "duplicate_mixed",
        n_records_ingroup > 1                      ~ "duplicate_ingroup",
        n_records_outgroup > 1                     ~ "duplicate_outgroup",
        TRUE                                       ~ "duplicate_other"
      ),
      is_duplicate_species_marker = dplyr::n() > 1,
      .groups = "drop"
    )

  tabla_raw <- tabla_raw |>
    dplyr::left_join(
      duplicate_groups |> dplyr::select(marker_key, species, n_records_total, n_records_ingroup, n_records_outgroup, classification_conflict, duplicate_class, is_duplicate_species_marker),
      by = c("marker_key", "species")
    )

  message("Exporting TABLE_raw_with_acceptance.csv ...")
  write_csv_dual(tabla_raw, legacy_name = "Tabla_raw_with_acceptance.csv", standard_name = "TABLE_raw_with_acceptance.csv")

  sequence_registry <- tabla_raw |>
    dplyr::select(marker_key, source_branch, species, fasta_name, sid, is_outgroup, matched_to_metadata, accepted_in_checklist, is_outgroup_species, accepted_or_outgroup, matched_to_checklist, record_status, is_duplicate_species_marker, n_records_total, n_records_ingroup, n_records_outgroup, duplicate_class, classification_conflict) |>
    dplyr::arrange(marker_key, species, source_branch, fasta_name)

  message("Exporting sequence_registry_with_acceptance.csv ...")
  write_csv_dual(sequence_registry, legacy_name = "TABLE_sequence_registry_with_acceptance.csv", standard_name = "TABLE_sequence_registry_with_acceptance.csv")

  duplicates <- dplyr::arrange(dplyr::filter(duplicate_groups, is_duplicate_species_marker), dplyr::desc(n_records_total), marker_key, species)
  message(sprintf("Duplicate species \U00d7 marker combinations: %d", nrow(duplicates)))

  message("Selecting one record per species X marker...")
  tabla_selected <- tabla_raw |>
    dplyr::filter(!is.na(species), !is.na(marker_key), !is.na(sid)) |>
    dplyr::mutate(sid = as.character(sid), fasta_name = as.character(fasta_name), accepted_rank = dplyr::if_else(accepted_or_outgroup, 1L, 0L), outgroup_rank = dplyr::if_else(is_outgroup_species, 1L, 0L)) |>
    dplyr::arrange(marker_key, species, dplyr::desc(accepted_rank), dplyr::desc(outgroup_rank), sid, fasta_name) |>
    dplyr::distinct(species, marker_key, .keep_all = TRUE) |>
    dplyr::select(species, marker_key, sid, fasta_name, source_branch, accepted_in_checklist, is_outgroup_species, accepted_or_outgroup)

  tabla_selected_final <- dplyr::filter(tabla_selected, accepted_or_outgroup)

  tabla_wide <- tabla_selected_final |>
    dplyr::mutate(species_class = dplyr::if_else(is_outgroup_species, "outgroup", "ingroup")) |>
    dplyr::select(species, species_class, accepted_or_outgroup, marker_key, sid) |>
    dplyr::distinct(species, marker_key, .keep_all = TRUE) |>
    tidyr::pivot_wider(names_from = marker_key, values_from = sid) |>
    dplyr::arrange(species_class, species)

  message("Exporting species_marker_sid_matrix.csv ...")
  write_csv_dual(tabla_wide, legacy_name = "Species_marker_SID_matrix.csv", standard_name = "TABLE_species_marker_sid_matrix.csv")

  duplicate_resolution <- tabla_raw |>
    dplyr::filter(!is.na(species), !is.na(marker_key)) |>
    dplyr::group_by(marker_key, species) |>
    dplyr::summarise(
      n_records_total    = dplyr::n(),
      n_records_ingroup  = sum(source_branch == "ingroup", na.rm = TRUE),
      n_records_outgroup = sum(source_branch == "outgroup", na.rm = TRUE),
      n_records_accepted = sum(accepted_or_outgroup %in% TRUE, na.rm = TRUE),
      n_records_rejected = sum(source_branch == "ingroup" & accepted_or_outgroup %in% FALSE, na.rm = TRUE),
      sid_candidates              = collapse_unique(sid),
      fasta_candidates            = collapse_unique(fasta_name),
      source_candidates           = collapse_unique(source_branch),
      accepted_status_candidates  = collapse_unique(ifelse(accepted_or_outgroup, "accepted_or_outgroup", "rejected")),
      classification_conflict     = dplyr::n_distinct(accepted_or_outgroup) > 1,
      duplicate_class = dplyr::case_when(
        dplyr::n() <= 1                            ~ "unique",
        n_records_ingroup > 0 & n_records_outgroup > 0 ~ "mixed_ingroup_outgroup",
        n_records_ingroup > 1                      ~ "ingroup_only",
        n_records_outgroup > 1                     ~ "outgroup_only",
        TRUE                                       ~ "other"
      ),
      .groups = "drop"
    ) |>
    dplyr::filter(n_records_total > 1) |>
    dplyr::left_join(
      dplyr::transmute(tabla_selected, marker_key, species, sid_selected = sid, fasta_selected = fasta_name, source_selected = source_branch, selected_accepted_or_outgroup = accepted_or_outgroup),
      by = c("marker_key", "species")
    ) |>
    dplyr::arrange(dplyr::desc(n_records_total), marker_key, species)

  message("Exporting duplicate_resolution_species_marker.csv ...")
  write_csv_dual(duplicate_resolution, legacy_name = "Duplicated_species_marker.csv", standard_name = "TABLE_duplicate_resolution_species_marker.csv")

  message("Mapping selected sequences back to FASTA headers...")
  tabla_fasta_map <- tabla_raw |>
    dplyr::select(marker_key, sid, fasta_name, species, source_branch) |>
    dplyr::filter(!is.na(sid), !is.na(fasta_name), !is.na(species)) |>
    dplyr::mutate(sid = as.character(sid)) |> dplyr::distinct()

  selected_fasta <- tabla_selected_final |>
    dplyr::select(species, marker_key, sid, source_branch) |>
    dplyr::left_join(tabla_fasta_map, by = c("species", "marker_key", "sid", "source_branch")) |>
    dplyr::group_by(species, marker_key, sid, source_branch) |>
    dplyr::summarise(fasta_name = dplyr::first(stats::na.omit(fasta_name)), .groups = "drop")

  missing_map <- dplyr::filter(selected_fasta, is.na(fasta_name))
  if (nrow(missing_map) > 0) {
    warning(sprintf("Some selected SIDs have no matching FASTA header (%d rows).", nrow(missing_map)), call. = FALSE)
    write_csv_dual(missing_map, legacy_name = "missing_map_selected_sids.csv", standard_name = "TABLE_missing_sid_header_map.csv")
  }

  message("Exporting cleaned FASTAs...")
  exported_markers <- character(0)
  markers_without_sequences <- character(0)
  exported_registry_list <- vector("list", length(markers_to_process))

  for (i in seq_along(markers_to_process)) {
    mk <- markers_to_process[i]
    ingroup_paths <- dplyr::pull(dplyr::filter(ingroup_index, marker_key == mk), path)
    outgroup_paths <- dplyr::pull(dplyr::filter(outgroup_index, marker_key == mk), path)
    
    ingroup_seq <- read_fasta_safe(ingroup_paths)
    outgroup_seq <- read_fasta_safe(outgroup_paths)
    combined <- c(ingroup_seq, outgroup_seq)
    
    if (length(combined) == 0) {
      markers_without_sequences <- c(markers_without_sequences, mk)
      next
    }
    
    keep_fasta_names <- unique(dplyr::pull(dplyr::filter(selected_fasta, marker_key == mk, !is.na(fasta_name)), fasta_name))
    if (length(keep_fasta_names) == 0) { message(sprintf("No accepted sequences for %s", mk)); next }
    
    keep_idx <- which(names(combined) %in% keep_fasta_names)
    if (length(keep_idx) == 0) { warning(sprintf("No matching FASTA headers found for marker %s", mk), call. = FALSE); next }
    
    aln_out <- combined[keep_idx]
    sp_order <- vapply(names(aln_out), extract_species, character(1))
    aln_out <- aln_out[order(sp_order, names(aln_out))]
    
    out_file <- file.path(out_marker_dir, paste0(mk, ".fasta"))
    Biostrings::writeXStringSet(aln_out, out_file)
    exported_markers <- c(exported_markers, mk)
    message(sprintf("Exported: %s (%d sequences)", basename(out_file), length(aln_out)))
    
    exported_registry_list[[i]] <- dplyr::tibble(marker_key = mk, fasta_name = names(aln_out), species = vapply(names(aln_out), extract_species, character(1)))
  }

  exported_registry <- dplyr::bind_rows(exported_registry_list) |>
    dplyr::left_join(dplyr::select(selected_fasta, marker_key, species, fasta_name, source_branch), by = c("marker_key", "species", "fasta_name"))

  message("Exporting marker taxon composition table...")
  raw_counts_group <- tabla_raw |> dplyr::group_by(marker_key, source_branch) |> dplyr::summarise(n_species = dplyr::n_distinct(species[!is.na(species)]), n_records = dplyr::n_distinct(fasta_name[!is.na(fasta_name)]), .groups = "drop") |> dplyr::mutate(stage = "raw", group = source_branch) |> dplyr::select(marker_key, stage, group, n_species, n_records)
  raw_counts_total <- tabla_raw |> dplyr::group_by(marker_key) |> dplyr::summarise(n_species = dplyr::n_distinct(species[!is.na(species)]), n_records = dplyr::n_distinct(fasta_name[!is.na(fasta_name)]), .groups = "drop") |> dplyr::mutate(stage = "raw", group = "total") |> dplyr::select(marker_key, stage, group, n_species, n_records)
  selected_counts_group <- tabla_selected_final |> dplyr::mutate(group = dplyr::if_else(is_outgroup_species, "outgroup", "ingroup")) |> dplyr::group_by(marker_key, group) |> dplyr::summarise(n_species = dplyr::n_distinct(species[!is.na(species)]), n_records = dplyr::n_distinct(sid[!is.na(sid)]), .groups = "drop") |> dplyr::mutate(stage = "selected") |> dplyr::select(marker_key, stage, group, n_species, n_records)
  selected_counts_total <- tabla_selected_final |> dplyr::group_by(marker_key) |> dplyr::summarise(n_species = dplyr::n_distinct(species[!is.na(species)]), n_records = dplyr::n_distinct(sid[!is.na(sid)]), .groups = "drop") |> dplyr::mutate(stage = "selected", group = "total") |> dplyr::select(marker_key, stage, group, n_species, n_records)
  exported_counts_group <- exported_registry |> dplyr::mutate(group = dplyr::if_else(source_branch == "outgroup", "outgroup", "ingroup")) |> dplyr::group_by(marker_key, group) |> dplyr::summarise(n_species = dplyr::n_distinct(species[!is.na(species)]), n_records = dplyr::n_distinct(fasta_name[!is.na(fasta_name)]), .groups = "drop") |> dplyr::mutate(stage = "exported") |> dplyr::select(marker_key, stage, group, n_species, n_records)
  exported_counts_total <- exported_registry |> dplyr::group_by(marker_key) |> dplyr::summarise(n_species = dplyr::n_distinct(species[!is.na(species)]), n_records = dplyr::n_distinct(fasta_name[!is.na(fasta_name)]), .groups = "drop") |> dplyr::mutate(stage = "exported", group = "total") |> dplyr::select(marker_key, stage, group, n_species, n_records)

  marker_taxon_composition <- dplyr::bind_rows(raw_counts_group, raw_counts_total, selected_counts_group, selected_counts_total, exported_counts_group, exported_counts_total) |> dplyr::arrange(marker_key, factor(stage, levels = c("raw", "selected", "exported")), group)
  write_csv_dual(marker_taxon_composition, legacy_name = "TABLE_marker_taxon_composition.csv", standard_name = "TABLE_marker_taxon_composition.csv")

  message("Exporting marker summary table...")
  raw_summary <- tabla_raw |> dplyr::group_by(marker_key) |> dplyr::summarise(n_records_raw_total = dplyr::n_distinct(fasta_name[!is.na(fasta_name)]), n_records_raw_ingroup = dplyr::n_distinct(fasta_name[source_branch == "ingroup" & !is.na(fasta_name)]), n_records_raw_outgroup = dplyr::n_distinct(fasta_name[source_branch == "outgroup" & !is.na(fasta_name)]), n_species_raw_total = dplyr::n_distinct(species[!is.na(species)]), n_species_raw_ingroup = dplyr::n_distinct(species[source_branch == "ingroup" & !is.na(species)]), n_species_raw_outgroup = dplyr::n_distinct(species[source_branch == "outgroup" & !is.na(species)]), n_species_raw_accepted_ingroup = dplyr::n_distinct(species[source_branch == "ingroup" & accepted_in_checklist & !is.na(species)]), n_species_raw_rejected_ingroup = dplyr::n_distinct(species[source_branch == "ingroup" & !accepted_in_checklist & !is.na(species)]), .groups = "drop")
  duplicate_summary <- duplicate_groups |> dplyr::group_by(marker_key) |> dplyr::summarise(n_species_marker_duplicated = sum(n_records_total > 1, na.rm = TRUE), n_duplicate_records_excess_total = sum(pmax(n_records_total - 1L, 0L), na.rm = TRUE), n_duplicate_records_excess_ingroup = sum(pmax(n_records_ingroup - 1L, 0L), na.rm = TRUE), n_duplicate_records_excess_outgroup = sum(pmax(n_records_outgroup - 1L, 0L), na.rm = TRUE), n_duplicate_records_excess_mixed = sum(dplyr::if_else(n_records_ingroup > 0 & n_records_outgroup > 0, n_records_total - 1L, 0L), na.rm = TRUE), n_classification_conflicts = sum(classification_conflict, na.rm = TRUE), .groups = "drop")
  selected_summary <- tabla_selected_final |> dplyr::group_by(marker_key) |> dplyr::summarise(n_species_selected_total = dplyr::n_distinct(species[!is.na(species)]), n_species_selected_ingroup = dplyr::n_distinct(species[!is.na(species) & !is_outgroup_species]), n_species_selected_outgroup = dplyr::n_distinct(species[!is.na(species) & is_outgroup_species]), n_records_selected_total = dplyr::n_distinct(sid[!is.na(sid)]), .groups = "drop")
  exported_summary <- exported_registry |> dplyr::left_join(dplyr::distinct(dplyr::select(tabla_selected_final, marker_key, species, source_branch, is_outgroup_species)), by = c("marker_key", "species", "source_branch")) |> dplyr::group_by(marker_key) |> dplyr::summarise(n_sequences_exported = dplyr::n_distinct(fasta_name[!is.na(fasta_name)]), n_species_exported_total = dplyr::n_distinct(species[!is.na(species)]), n_species_exported_ingroup = dplyr::n_distinct(species[!is.na(species) & !is_outgroup_species]), n_species_exported_outgroup = dplyr::n_distinct(species[!is.na(species) & is_outgroup_species]), .groups = "drop")

  tabla_summary <- raw_summary |> dplyr::left_join(duplicate_summary, by = "marker_key") |> dplyr::left_join(selected_summary, by = "marker_key") |> dplyr::left_join(exported_summary, by = "marker_key") |> dplyr::mutate(n_species_marker_duplicated = dplyr::coalesce(n_species_marker_duplicated, 0L), n_duplicate_records_excess_total = dplyr::coalesce(n_duplicate_records_excess_total, 0L), n_duplicate_records_excess_ingroup = dplyr::coalesce(n_duplicate_records_excess_ingroup, 0L), n_duplicate_records_excess_outgroup = dplyr::coalesce(n_duplicate_records_excess_outgroup, 0L), n_duplicate_records_excess_mixed = dplyr::coalesce(n_duplicate_records_excess_mixed, 0L), n_classification_conflicts = dplyr::coalesce(n_classification_conflicts, 0L), n_species_selected_total = dplyr::coalesce(n_species_selected_total, 0L), n_species_selected_ingroup = dplyr::coalesce(n_species_selected_ingroup, 0L), n_species_selected_outgroup = dplyr::coalesce(n_species_selected_outgroup, 0L), n_records_selected_total = dplyr::coalesce(n_records_selected_total, 0L), n_sequences_exported = dplyr::coalesce(n_sequences_exported, 0L), n_species_exported_total = dplyr::coalesce(n_species_exported_total, 0L), n_species_exported_ingroup = dplyr::coalesce(n_species_exported_ingroup, 0L), n_species_exported_outgroup = dplyr::coalesce(n_species_exported_outgroup, 0L), pct_species_raw_accepted_ingroup = dplyr::if_else(n_species_raw_ingroup > 0, round(100 * n_species_raw_accepted_ingroup / n_species_raw_ingroup, 2), NA_real_), pct_species_retained_from_raw = dplyr::if_else(n_species_raw_total > 0, round(100 * n_species_exported_total / n_species_raw_total, 2), NA_real_), pct_records_retained_from_raw = dplyr::if_else(n_records_raw_total > 0, round(100 * n_sequences_exported / n_records_raw_total, 2), NA_real_)) |> dplyr::arrange(dplyr::desc(n_species_raw_total), marker_key)
  write_csv_dual(tabla_summary, legacy_name = "Marker_summary_table.csv", standard_name = "TABLE_marker_summary.csv")

  log_lines <- c("FINAL MARKER INTEGRATION LOG", sprintf("Input ingroup directory: %s", ingroup_dir), sprintf("Input outgroup directory: %s", outgroup_dir), sprintf("Output directory: %s", output_dir), sprintf("Markers in ingroup index: %d", nrow(ingroup_index)), sprintf("Markers in outgroup index: %d", nrow(outgroup_index)), sprintf("Unique ingroup-driven markers processed: %d", length(markers_to_process)), sprintf("Rows in raw integration table: %d", nrow(tabla_raw)), sprintf("Unmatched metadata rows: %d", nrow(unmatched_metadata)), sprintf("Duplicate species \U00d7 marker combinations: %d", nrow(duplicates)), sprintf("Selected unique species \U00d7 marker rows (accepted final): %d", nrow(tabla_selected_final)), sprintf("Missing SID-to-header mappings: %d", nrow(missing_map)), sprintf("Exported FASTA markers: %d", length(exported_markers)), sprintf("Markers with zero combined sequences: %d", length(markers_without_sequences)))
  writeLines(log_lines, con = file.path(output_dir, "LOG_clean_integration_summary.txt"))
  
  message("\nMarker integration and cleaning pipeline completed. ")
}
