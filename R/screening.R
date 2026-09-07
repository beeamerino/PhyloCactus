#' Test a Locus Alignment for Substitution Saturation via Corrected-vs-Uncorrected Distance Regression
#'
#' Internal helper used by `run_marker_screening()`. Regresses uncorrected pairwise sequence distance
#' against a model-corrected distance; a regression slope below `saturation_flag_cutoff` indicates that
#' additional substitutions are being masked by multiple hits at the same site (saturation), reducing
#' phylogenetic informativeness at deeper divergences. Extracted as a standalone, independently testable
#' package-internal function (previously a closure nested inside `run_marker_screening()`).
#'
#' @param dna_bin An object of class `DNAbin` (`ape`) representing the locus alignment to test.
#' @param saturation_method Character. `"corrected"` regresses uncorrected p-distance (`raw`) against a
#'   Gamma-corrected K80 distance; `"legacy"` regresses Gamma-corrected K80 against equal-rates K80
#'   (retained for backward comparability; sensitive to a known mathematical artifact inflating the slope
#'   toward >= 1.0). Defaults to `"corrected"`.
#' @param saturation_flag_cutoff Numeric. Slope threshold below which the locus is flagged as saturated. Defaults to `0.3`.
#' @return A list with `slope` (numeric regression slope or `NA`), `saturated` (logical or `NA`), and
#'   `reason` (character diagnostic code: `"ok"`, `"insufficient_data"`, `"variance_na"`, or `"regression_failed"`).
#' @noRd
.test_saturation_proxy <- function(dna_bin, saturation_method = c("corrected", "legacy"), saturation_flag_cutoff = 0.3) {
  saturation_method <- match.arg(saturation_method)

  if (saturation_method == "legacy") {
    dist_1 <- tryCatch(as.numeric(ape::dist.dna(dna_bin, model = "K80", pairwise.deletion = TRUE)), error = function(e) NA_real_)
    dist_2 <- tryCatch(as.numeric(ape::dist.dna(dna_bin, model = "K80", pairwise.deletion = TRUE, gamma = TRUE)), error = function(e) NA_real_)
  } else {
    dist_1 <- tryCatch(as.numeric(ape::dist.dna(dna_bin, model = "K80", pairwise.deletion = TRUE, gamma = TRUE)), error = function(e) NA_real_)
    dist_2 <- tryCatch(as.numeric(ape::dist.dna(dna_bin, model = "raw", pairwise.deletion = TRUE)), error = function(e) NA_real_)
  }

  ok <- !(is.na(dist_1) | is.na(dist_2))
  dist_corrected <- dist_1[ok]
  dist_uncorrected <- dist_2[ok]

  if (length(dist_corrected) < 5) return(list(slope = NA_real_, saturated = NA, reason = "insufficient_data"))
  df <- data.frame(x = dist_corrected, y = dist_uncorrected)
  if (is.na(stats::var(df$x)) || is.na(stats::var(df$y))) return(list(slope = NA_real_, saturated = NA, reason = "variance_na"))
  fit <- tryCatch(stats::lm(y ~ x, data = df), error = function(e) NULL)
  if (is.null(fit)) return(list(slope = NA_real_, saturated = NA, reason = "regression_failed"))
  slope <- unname(stats::coef(fit)[2])
  list(slope = slope, saturated = isTRUE(slope < saturation_flag_cutoff), reason = "ok")
}

#' Fraction of outgroup k-mers present in the ingroup sequences of the same marker
#'
#' Two GenBank annotations can carry the same gene name and describe different regions: paralogues
#' of one family, a gene against the spacer beside it, or an amplicon that a second study placed
#' elsewhere. Aligners do not refuse such a pair. MAFFT will return a block, the block will look
#' like a marker with coverage on both sides of the root, and its columns will carry no positional
#' homology. The 2026-09-01 audit found this twice: `pepC` (0.000 of 20-mers shared) and
#' `trnT-psbD` (0.039), the second in a locus readmitted specifically to connect the outgroup.
#'
#' The measure is alignment-free, which is the point: it cannot be rescued by a better alignment,
#' so it separates "these sequences are hard to align" from "these sequences are not the same
#' region". Genuine counterparts in this dataset return 0.28 to 0.60 at `k = 20`; the two false
#' ones returned 0.000 and 0.039.
#'
#' Orientation is not guaranteed to agree between two independent phylotaR runs, so an outgroup
#' k-mer counts as shared when either it or its reverse complement occurs in the ingroup pool.
#' Only the outgroup side is reverse-complemented, that being the smaller of the two.
#'
#' @param ingroup_paths,outgroup_paths Character vectors of FASTA paths for one marker.
#' @param k Integer. k-mer length. 20 is strict enough that unrelated sequences share nothing;
#'   10 is permissive enough that a genuinely divergent homologue still registers.
#' @param max_seqs Integer. Sequences kept per side, to bound the cost on the large loci. The
#'   subsample is evenly spaced rather than random, so the figure is reproducible without a seed
#'   and without disturbing the caller's RNG state.
#' @return A list with `share` (the better of the two orientations, or `NA` if either side is
#'   empty), `share_forward`, `share_revcomp` and `n_kmers_outgroup`. The two directions are
#'   reported separately because a whole outgroup pool deposited on the opposite strand looks
#'   identical, in the combined figure, to one that agrees: both come back healthy. Only the split
#'   distinguishes "same region, other strand" from "same region, same strand".
#' @keywords internal
#' @noRd
.marker_homology_share <- function(ingroup_paths, outgroup_paths, k = 20L,
                                   max_seqs = 200L) {
  read_pool <- function(paths) {
    seqs <- unlist(lapply(paths, function(p) {
      as.character(Biostrings::readDNAStringSet(p))
    }), use.names = FALSE)
    seqs <- gsub("-", "", toupper(seqs), fixed = TRUE)
    seqs <- seqs[nchar(seqs) >= k]
    if (length(seqs) > max_seqs) {
      seqs <- seqs[unique(round(seq(1L, length(seqs), length.out = max_seqs)))]
    }
    seqs
  }
  kmer_set <- function(seqs) {
    out <- unlist(lapply(seqs, function(s) {
      n <- nchar(s) - k + 1L
      if (n < 1L) return(character(0))
      substring(s, seq_len(n), seq_len(n) + k - 1L)
    }), use.names = FALSE)
    out <- unique(out)
    out[!grepl("[^ACGT]", out)]
  }

  none <- list(share = NA_real_, share_forward = NA_real_, share_revcomp = NA_real_,
               n_kmers_outgroup = 0L)
  ing <- read_pool(ingroup_paths)
  outg <- read_pool(outgroup_paths)
  if (length(ing) == 0L || length(outg) == 0L) return(none)
  pool_in <- kmer_set(ing)
  ko <- kmer_set(outg)
  if (length(ko) == 0L || length(pool_in) == 0L) {
    none$n_kmers_outgroup <- length(ko)
    return(none)
  }
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(ko)))
  fwd <- mean(ko %in% pool_in)
  rev <- mean(rc %in% pool_in)
  list(share = max(fwd, rev), share_forward = fwd, share_revcomp = rev,
       n_kmers_outgroup = length(ko))
}

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
#' @param min_nseq_to_retain Integer. Minimum number of ingroup sequences required per locus, counted after outlier filtering. Defaults to `100L`.
#' @param outgroup_folder Character or `NULL`. Directory of aligned outgroup locus FASTA files, matched to the ingroup files by marker name. Used for **reporting only**: the summary table gains `n_outgroup` and `n_total` so a locus rejected here can be seen to carry outgroup data, but no retention decision depends on them. Outgroup markers with no ingroup counterpart simply leave those two columns empty; which ones they are, and what becomes of them, is reported by [integrate_and_clean_markers()], which owns that decision and applies `marker_aliases` before comparing the two sets. Defaults to `NULL`.
#'
#'   This module asks whether a locus resolves the ingroup radiation, and the answer cannot depend on how many outgroup accessions exist. `trnT-psbD` is the case that forced the distinction: 50 ingroup sequences against a threshold of 100, correctly rejected as an ingroup marker, while carrying 49 Portulaca accessions of 1347 bp that are the best outgroup coverage in the dataset. Bringing it back is a decision about connecting the two groups, which belongs to [integrate_and_clean_markers()] and its `readmit_markers` argument, not to a sequence count here.
#' @param max_marker_missing Numeric. Maximum allowable missing data fraction per locus. Defaults to `0.7`.
#' @param saturation_flag_cutoff Numeric. Uncorrected p-distance vs. raw distance slope threshold to flag substitution saturation. Defaults to `0.3`.
#' @param saturation_keep_cutoff Numeric. Saturation slope cutoff threshold below which saturated loci are excluded. Defaults to `0.5`.
#' @param iqr_multiplier Numeric. Interquartile range (IQR) multiplier for identifying site-length outlier bounds. Defaults to `1.5`.
#' @param saturation_method Character. Method for saturation test distance calculation. Defaults to `"corrected"`.
#' @param gap_handling Character. Method for handling ambiguous and gap characters. Defaults to `"iupac"`.
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
  outgroup_folder = NULL,
  max_marker_missing = 0.7,
  saturation_flag_cutoff = 0.3,
  saturation_keep_cutoff = 0.5,
  iqr_multiplier = 1.5,
  saturation_method = c("corrected", "legacy"),
  gap_handling = c("iupac", "legacy")
) {
  saturation_method <- match.arg(saturation_method)
  gap_handling <- match.arg(gap_handling)
  
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
    if (saturation_method == "legacy") {
      dist_1 <- ape::dist.dna(dna_bin, model = "K80", pairwise.deletion = TRUE)
      dist_2 <- ape::dist.dna(dna_bin, model = "K80", pairwise.deletion = TRUE, gamma = TRUE)
      xlab_txt <- "K80 pairwise distance"
      ylab_txt <- "K80 pairwise distance (gamma)"
    } else {
      dist_1 <- ape::dist.dna(dna_bin, model = "K80", pairwise.deletion = TRUE, gamma = TRUE)
      dist_2 <- ape::dist.dna(dna_bin, model = "raw", pairwise.deletion = TRUE)
      xlab_txt <- "Model-corrected distance (K80+Gamma)"
      ylab_txt <- "Uncorrected p-distance"
    }
    
    df <- data.frame(x = as.numeric(dist_1), y = as.numeric(dist_2))
    
    ggplot2::ggplot(df, ggplot2::aes(x, y)) +
      ggplot2::geom_point(size = 1.2) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::ggtitle(paste(marker, title_suffix)) +
      ggplot2::xlab(xlab_txt) + ggplot2::ylab(ylab_txt) +
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

  # Sequence counts of the outgroup alignments, keyed by marker name, to fill n_outgroup and
  # n_total for the loci that have a counterpart. Nothing here decides anything: this module asks
  # whether a locus resolves the ingroup radiation, which cannot depend on how many outgroup
  # accessions exist.
  #
  # Which outgroup markers have no ingroup counterpart is deliberately NOT reported here. That is
  # integrate_and_clean_markers() (Stage 4), which owns the decision, applies `marker_aliases`
  # before comparing, and writes tables/TABLE_outgroup_markers_dropped.csv. Reporting it here as
  # well duplicated that list without the aliases, so `trnL` was named as an orphan when Stage 4
  # renames it onto trnL_trnF and keeps it.
  outgroup_counts <- integer(0)
  if (!is.null(outgroup_folder)) {
    stop_if_missing_dir(outgroup_folder, "outgroup_folder")
    og_files <- list.files(outgroup_folder, pattern = "^ALN_masked_final_.*\\.fasta$", full.names = TRUE)
    if (length(og_files) == 0L) {
      warning("No aligned outgroup files found in '", outgroup_folder,
              "'; n_outgroup and n_total are left empty.", call. = FALSE)
    } else {
      outgroup_counts <- vapply(og_files, function(f) {
        length(Biostrings::readDNAStringSet(f))
      }, integer(1))
      names(outgroup_counts) <- vapply(og_files, extract_marker_name, character(1))
    }
  }

  # Statistical saturation test itself is defined at package level as .test_saturation_proxy()
  # (outside this function) for independent testability; bind it here to the already-validated
  # saturation_method so downstream call sites are unchanged.
  test_saturation_proxy <- function(dna_bin, saturation_flag_cutoff = 0.3) {
    .test_saturation_proxy(dna_bin, saturation_method = saturation_method, saturation_flag_cutoff = saturation_flag_cutoff)
  }

  alignment_stats <- function(dna_bin) {
    aln_mat <- toupper(safe_as_matrix(dna_bin))
    gaps_mask_orig <- aln_mat == "-"
    nseq <- nrow(aln_mat); total_cols <- ncol(aln_mat)
    aln_base_only <- aln_mat
    aln_base_only[!aln_base_only %in% c("A", "C", "G", "T")] <- NA
    informative_col <- apply(aln_base_only, 2, function(col) any(!is.na(col)))
    aln_len_effective <- sum(informative_col)
    missing_prop <- sum(aln_mat %in% c("N", "?", "X"), na.rm = TRUE) / (nseq * total_cols)
    gap_prop <- sum(gaps_mask_orig, na.rm = TRUE) / (nseq * total_cols)
    # IUPAC ambiguity codes preserved as literal characters under gap_handling = "iupac" are
    # neither "missing" (N/?/X) nor "gap" (-); tracked separately so the QC statistics fully
    # account for every non-diagnostic (non-ACGT) site in the alignment.
    iupac_ambiguous <- c("R", "Y", "S", "W", "K", "M", "B", "D", "H", "V")
    ambiguous_prop <- sum(aln_mat %in% iupac_ambiguous, na.rm = TRUE) / (nseq * total_cols)
    n_bases <- sum(aln_base_only %in% c("A", "C", "G", "T"), na.rm = TRUE)
    gc_prop <- if (n_bases == 0) NA_real_ else sum(aln_base_only == "G" | aln_base_only == "C", na.rm = TRUE) / n_bases
    var_sites <- sum(apply(aln_base_only, 2, function(col) length(unique(stats::na.omit(col))) > 1))
    pi_sites <- sum(apply(aln_base_only, 2, function(col) { tt <- table(stats::na.omit(col)); any(tt >= 2) && length(tt) >= 2 }))
    ent_vals <- apply(aln_base_only, 2, function(col) {
      tt <- table(stats::na.omit(col))
      if (length(tt) == 0) return(NA_real_)
      p <- tt / sum(tt); -sum(p * log(p))
    })
    list(aln_len = aln_len_effective, total_cols = total_cols, missing = missing_prop, gaps = gap_prop, ambiguous = ambiguous_prop, gc = gc_prop, var_sites = var_sites, pi_sites = pi_sites, entropy = mean(ent_vals, na.rm = TRUE))
  }
  
  make_summary_row <- function(marker, status, decision, decision_reason, n_original = NA_integer_, n_filtered = NA_integer_, n_removed = NA_integer_, n_ingroup = NA_integer_, n_outgroup = NA_integer_, n_total = NA_integer_, aln_len_before = NA_real_, aln_len_after = NA_real_, missing_before = NA_real_, missing_after = NA_real_, ambiguous_before = NA_real_, ambiguous_after = NA_real_, gc_before = NA_real_, gc_after = NA_real_, var_sites_before = NA_real_, var_sites_after = NA_real_, pi_sites_before = NA_real_, pi_sites_after = NA_real_, entropy_before = NA_real_, entropy_after = NA_real_, slope_before = NA_real_, slope_after = NA_real_, saturated_before = NA, saturated_after = NA, reason_before = NA_character_, reason_after = NA_character_, removed = "-") {
    data.frame(marker = marker, status = status, decision = decision, decision_reason = decision_reason, n_original = n_original, n_filtered = n_filtered, n_removed = n_removed, n_ingroup = n_ingroup, n_outgroup = n_outgroup, n_total = n_total, aln_len_before = aln_len_before, aln_len_after = aln_len_after, missing_before = missing_before, missing_after = missing_after, ambiguous_before = ambiguous_before, ambiguous_after = ambiguous_after, gc_before = gc_before, gc_after = gc_after, var_sites_before = var_sites_before, var_sites_after = var_sites_after, pi_sites_before = pi_sites_before, pi_sites_after = pi_sites_after, entropy_before = entropy_before, entropy_after = entropy_after, slope_before = slope_before, slope_after = slope_after, saturated_before = saturated_before, saturated_after = saturated_after, reason_before = reason_before, reason_after = reason_after, removed = removed, stringsAsFactors = FALSE)
  }
  
  results <- list(); filtered_alignments <- list(); pw_before_list <- list(); sat_before_list <- list(); pw_after_list <- list(); sat_after_list <- list()
  
  for (f in fasta_files) {
    marker <- extract_marker_name(f)
    log_message("Processing marker: ", marker)
    res <- tryCatch({
      aln <- ape::read.dna(f, format = "fasta")
      aln_char <- toupper(safe_as_matrix(aln))
      if (gap_handling == "legacy") {
        aln_char[!(aln_char %in% c("A", "C", "G", "T"))] <- "-"
      } else {
        valid_iupac <- c("A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N", "-")
        aln_char[aln_char %in% c("?", "X")] <- "N"
        aln_char[!aln_char %in% valid_iupac] <- "N"
      }
      aln <- ape::as.DNAbin(aln_char)
      seq_names <- rownames(aln)
      if (is.null(seq_names) || length(seq_names) != nrow(aln)) seq_names <- paste0("seq_", seq_len(nrow(aln))); rownames(aln) <- seq_names
      if (ncol(aln) < min_cols_to_evaluate) {
        log_message("Marker skipped for short alignment (<", min_cols_to_evaluate, " columns): ", marker)
        list(summary_row = make_summary_row(marker = marker, status = "skipped_short_alignment", decision = "NO", decision_reason = paste0("alignment_shorter_than_", min_cols_to_evaluate, "_columns"), n_original = nrow(aln), n_filtered = nrow(aln), n_removed = 0L, aln_len_before = ncol(aln), aln_len_after = ncol(aln), removed = "-"), pw_before = NULL, sat_before = NULL, pw_after = NULL, sat_after = NULL, filtered_alignment = NULL)
      } else {
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
      # Reported, never decisive. The retention criterion is the ingroup count, as it has always
      # been, because this module measures whether a locus resolves the ingroup radiation.
      n_ingroup  <- nrow(aln_filtered)
      n_outgroup <- if (marker %in% names(outgroup_counts)) unname(outgroup_counts[[marker]]) else 0L
      n_total    <- n_ingroup + n_outgroup

      keep_marker <- (!isTRUE(sat_after$saturated) && !is.na(sat_after$slope) && sat_after$slope > saturation_keep_cutoff && stats_after$missing < max_marker_missing && stats_after$var_sites > 0 && stats_after$aln_len >= min_aln_len_to_retain && n_ingroup >= min_nseq_to_retain)
      decision_reason <- if (keep_marker) "passed_all_thresholds" else paste(c(if (isTRUE(sat_after$saturated)) "flagged_as_saturated" else NULL, if (is.na(sat_after$slope)) "slope_na" else NULL, if (!is.na(sat_after$slope) && sat_after$slope <= saturation_keep_cutoff) paste0("slope_le_", saturation_keep_cutoff) else NULL, if (stats_after$missing >= max_marker_missing) paste0("missing_ge_", max_marker_missing) else NULL, if (stats_after$var_sites <= 0) "no_variable_sites" else NULL, if (stats_after$aln_len < min_aln_len_to_retain) paste0("aln_len_lt_", min_aln_len_to_retain) else NULL, if (n_ingroup < min_nseq_to_retain) paste0("nseq_lt_", min_nseq_to_retain) else NULL), collapse = ";")
      summary_row <- make_summary_row(marker = marker, status = if (keep_marker) "processed_retained" else "processed_rejected", decision = if (keep_marker) "YES" else "NO", decision_reason = decision_reason, n_original = nrow(aln), n_filtered = nrow(aln_filtered), n_removed = length(outliers), n_ingroup = n_ingroup, n_outgroup = n_outgroup, n_total = n_total, aln_len_before = stats_before$aln_len, aln_len_after = stats_after$aln_len, missing_before = stats_before$missing, missing_after = stats_after$missing, ambiguous_before = stats_before$ambiguous, ambiguous_after = stats_after$ambiguous, gc_before = stats_before$gc, gc_after = stats_after$gc, var_sites_before = stats_before$var_sites, var_sites_after = stats_after$var_sites, pi_sites_before = stats_before$pi_sites, pi_sites_after = stats_after$pi_sites, entropy_before = stats_before$entropy, entropy_after = stats_after$entropy, slope_before = sat_before$slope, slope_after = sat_after$slope, saturated_before = sat_before$saturated, saturated_after = sat_after$saturated, reason_before = sat_before$reason, reason_after = sat_after$reason, removed = if (length(outliers) == 0) "-" else paste(outliers, collapse = ";"))
      list(summary_row = summary_row, pw_before = pw_before, sat_before = sat_before_plot, pw_after = pw_after, sat_after = sat_after_plot, filtered_alignment = if (keep_marker) aln_filtered else NULL)
      }
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
  utils::write.csv(summary_table, out_summary_csv, row.names = FALSE)
  
  grDevices::cairo_pdf(filename = out_pdf, width = 11, height = 8.5, pointsize = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
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
  
  selected_markers <- dplyr::pull(dplyr::filter(summary_table, decision == "YES"), marker)
  for (marker in names(filtered_alignments)) {
    out_fasta <- file.path(dir_markers, paste0(marker, ".fasta"))
    ape::write.dna(filtered_alignments[[marker]], file = out_fasta, format = "fasta", nbcol = -1, colw = 9999)
  }
  
  saveRDS(list(pw_before = pw_before_list, pw_after = pw_after_list, sat_before = sat_before_list, sat_after = sat_after_list), file = file.path(out_base, "RDS_marker_screening_plots.rds"))
  writeLines(utils::capture.output(utils::sessionInfo()), con = out_sessioninfo_txt)
  
  message("\nMarker screening and quality check completed. ")
  return(summary_table)
}

#' Final Marker Integration, Taxonomic Cleaning, and Ingroup-Outgroup Partitioning
#'
#' Integrates independently curated ingroup and outgroup sequence datasets, standardizes species binomials
#' against the authoritative taxonomic checklist (Korotkova et al. 2021), and exports decoupled
#' FASTA sequence directories (`cleaned_markers_ingroup/`, `cleaned_markers_outgroup/`, `cleaned_markers_joint/`).
#' Computes isolated molecular informativeness metrics (tips, length, variable sites, parsimony informative sites,
#' GC content, and missingness) for the focal ingroup radiation to avoid outgroup-driven inflation artifacts.
#'
#' @param ingroup_dir Character. Path to directory containing filtered ingroup FASTA alignments.
#' @param outgroup_dir Character. Path to directory containing filtered outgroup FASTA alignments.
#' @param output_dir Character. Root destination directory for cleaned FASTA outputs and diagnostic tables.
#' @param accepted_list_file Character. Path to accepted botanical checklist CSV or Excel file (e.g., `CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx`).
#' @param metadata_in_file Character. Path to ingroup accession occupancy metadata CSV file.
#' @param metadata_out_file Character. Path to outgroup accession occupancy metadata CSV file.
#' @param marker_aliases Named character vector or `NULL`. Renames outgroup markers onto their ingroup counterpart before the two sets are joined, as `c(trnL = "trnL_trnF")`. Names and values are normalised the same way the filenames are, so either spelling works.
#'
#'   The two `phylotaR` runs cluster independently, so the same region can be named differently in each: the outgroup carries two `trnL` accessions of the trnL intron, which is the first half of the ingroup `trnL-trnF` amplicon, and they share 63% of their 20-mers with it. Aliasing is a homology claim and has to be evidenced, not assumed from the name: outgroup `ndhF` shares no 20-mer with ingroup `ndhF-rpl32` despite the obvious resemblance, because the two datasets cover different parts of the gene. Defaults to `NULL`.
#' @param readmit_markers Character vector or `NULL`. Markers that [run_marker_screening()] rejected and that are brought back into the joint dataset because they connect the outgroup to the ingroup, read from `readmit_dir`. This is a different question from the one Module 3 answers, and it is stated as such: `trnT-psbD` resolves the ingroup poorly, with 50 sequences over 1008 Cactaceae terminals, and carries the best outgroup coverage of the dataset, 49 *Portulaca* accessions of 1347 bp against the five or six terminals of overlap the retained loci provide. Retaining it for ingroup resolution would be wrong; retaining it for rooting is the point. Defaults to `NULL`.
#' @param readmit_dir Character or `NULL`. Directory holding the pre-screening ingroup alignments, normally `2_MAFFT_Cactaceae/alignments`. Required when `readmit_markers` is supplied. Defaults to `NULL`.
#' @param drop_report Logical. Write `tables/TABLE_outgroup_markers_dropped.csv` listing the outgroup markers that have no ingroup counterpart, with their sequence counts, and warn about them. The joint marker set is the ingroup's, so those markers leave the analysis at this point; a run takes hours and a console warning alone scrolls away. Defaults to `TRUE`.
#' @param homology_check Logical. For every marker present on both sides, measure the fraction of outgroup k-mers that occur in the ingroup sequences of the same name, write `tables/TABLE_marker_homology_check.csv`, and warn about the markers that share almost none. The measure is alignment-free, so it distinguishes sequences that are hard to align from sequences that are not the same region. Reported, never blocking. Defaults to `TRUE`.
#' @param homology_k20_min,homology_k10_min Numeric. A marker is flagged `no_detectable_homology` when it falls below **both**, at `k = 20` and `k = 10` respectively. Genuine counterparts in this dataset return 0.28 to 0.60 at `k = 20`; the two false pairs found on 2026-09-01 returned 0.000 (`pepC`, two PEPC paralogues) and 0.039 (`trnT-psbD`, ingroup median 598 bp against outgroup 1347 bp). Default to `0.05` and `0.15`.
#' @return A data frame containing the comprehensive marker summary with decoupled ingroup and joint metrics.
#' @references
#' Korotkova, N., Aquino, D., Arias, S., Eggli, U., Franck, A., Gómez-Hinostrosa, C., Guerrero, P. C.,
#' Hernández, H. M., Kohlbecker, A., Köhler, M., Luna, R., Machado, M., Merclinger, M., Nyffeler, R.,
#' Salvador-Montiel, S., Sánchez, D., Schlumpberger, B. O., & Berendsohn, W. G. (2021). Cactaceae at
#' Caryophyllales.org - a dynamic online species-level taxonomic backbone for the family.
#' *Willdenowia*, 51(2), 251–270. \doi{10.3372/wi.51.51208}
#' @export
integrate_and_clean_markers <- function(
  ingroup_dir,
  outgroup_dir,
  output_dir,
  accepted_list_file,
  metadata_in_file,
  metadata_out_file,
  marker_aliases = NULL,
  readmit_markers = NULL,
  readmit_dir = NULL,
  drop_report = TRUE,
  homology_check = TRUE,
  homology_k20_min = 0.05,
  homology_k10_min = 0.15
) {
  
  # Argument validation before any I/O, so an incoherent call fails on the argument rather than on
  # whichever directory scan happens to run first.
  if (!is.null(readmit_markers) && length(readmit_markers) > 0L && is.null(readmit_dir)) {
    stop("`readmit_markers` needs `readmit_dir`, the directory of pre-screening ingroup ",
         "alignments (normally \"2_MAFFT_Cactaceae/alignments\").", call. = FALSE)
  }

  out_tables_dir          <- file.path(output_dir, "tables")
  out_logs_dir            <- file.path(output_dir, "logs")
  out_marker_joint_dir    <- file.path(output_dir, "cleaned_markers_joint")
  out_marker_ingroup_dir  <- file.path(output_dir, "cleaned_markers_ingroup")
  out_marker_outgroup_dir <- file.path(output_dir, "cleaned_markers_outgroup")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_logs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_marker_joint_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_marker_ingroup_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_marker_outgroup_dir, recursive = TRUE, showWarnings = FALSE)
  
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
    x |>
      stringr::str_replace_all("\u00d7", "") |>
      stringr::str_replace_all("\\bx\\s+", "") |>
      stringr::str_trim() |>
      stringr::str_replace_all("[ -]+", "_")
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

  write_clean_table <- function(df, table_name) {
    table_path <- file.path(out_tables_dir, table_name)
    readr::write_csv(df, table_path, na = "")
  }

  collapse_unique <- function(x, sep = ";") {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    x <- sort(unique(x))
    if (length(x) == 0) return(NA_character_)
    paste(x, collapse = sep)
  }

  # `metrics_status` makes explicit why a row may carry empty site counts. The joint FASTA of a
  # marker is, before Stage 5, the union of two independently masked alignments whose widths do
  # not match, so it is a sequence set and not an alignment: per-column statistics (variable,
  # parsimony-informative and invariant sites) are undefined, and the composition percentages
  # (gaps, missing, GC) mix sequences of very different lengths and are therefore not comparable
  # with the ingroup-only or outgroup-only rows. Those non-comparable percentages are now
  # suppressed rather than published silently. Use TABLE_partition_informativeness_comparison.csv
  # (Stage 6) for the real joint informativeness of each partition.
  compute_alignment_metrics <- function(seqs_dna, marker_name, group_label) {
    if (length(seqs_dna) == 0) {
      return(dplyr::tibble(
        marker = marker_name,
        group = group_label,
        n_sequences = 0L,
        n_species = 0L,
        min_length = NA_real_,
        max_length = NA_real_,
        mean_length = NA_real_,
        missing_pct = NA_real_,
        gaps_pct = NA_real_,
        gc_pct = NA_real_,
        variable_sites = NA_integer_,
        parsimony_informative_sites = NA_integer_,
        invariant_sites = NA_integer_,
        metrics_status = "empty_no_sequences"
      ))
    }

    chars <- as.character(seqs_dna)
    lens <- nchar(chars)
    sp <- vapply(names(seqs_dna), extract_species, character(1))
    n_sp <- length(unique(sp[!is.na(sp) & nzchar(sp)]))
    
    is_aligned <- length(unique(lens)) == 1
    split_chars <- strsplit(toupper(chars), "", fixed = TRUE)
    all_cells <- unlist(split_chars)
    bases <- all_cells[all_cells %in% c("A", "C", "G", "T")]
    gc_pct <- if (length(bases) == 0) NA_real_ else round(sum(bases %in% c("G", "C")) / length(bases) * 100, 3)
    missing_pct <- round(mean(all_cells %in% c("N", "?")) * 100, 3)
    gaps_pct <- round(mean(all_cells == "-") * 100, 3)
    
    if (is_aligned && lens[1] > 0) {
      mat <- do.call(rbind, split_chars)
      var_sites <- sum(apply(mat, 2, function(col) {
        col_clean <- col[col %in% c("A", "C", "G", "T")]
        length(unique(col_clean)) > 1
      }))
      pi_sites <- sum(apply(mat, 2, function(col) {
        col_clean <- col[col %in% c("A", "C", "G", "T")]
        if (length(col_clean) == 0) return(FALSE)
        sum(table(col_clean) >= 2) >= 2
      }))
      inv_sites <- sum(apply(mat, 2, function(col) {
        col_clean <- col[col %in% c("A", "C", "G", "T")]
        length(col_clean) > 0 && length(unique(col_clean)) == 1
      }))
      metrics_status <- "ok_aligned"
    } else {
      var_sites <- NA_integer_
      pi_sites <- NA_integer_
      inv_sites <- NA_integer_
      # Composition percentages are computed over a ragged set of sequences and would not be
      # comparable with the aligned rows, so they are suppressed together with the site counts.
      gc_pct <- NA_real_
      missing_pct <- NA_real_
      gaps_pct <- NA_real_
      metrics_status <- sprintf(
        "unaligned_unequal_lengths_%d_to_%d_site_and_composition_metrics_not_computable",
        min(lens), max(lens)
      )
    }

    dplyr::tibble(
      marker = marker_name,
      group = group_label,
      n_sequences = length(seqs_dna),
      n_species = n_sp,
      min_length = min(lens),
      max_length = max(lens),
      mean_length = round(mean(lens), 1),
      missing_pct = missing_pct,
      gaps_pct = gaps_pct,
      gc_pct = gc_pct,
      variable_sites = var_sites,
      parsimony_informative_sites = pi_sites,
      invariant_sites = inv_sites,
      metrics_status = metrics_status
    )
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

  # Applied to the outgroup only, and before the join, because the join is by marker name and the
  # two phylotaR runs name their clusters independently.
  apply_aliases <- function(keys) {
    if (is.null(marker_aliases) || length(marker_aliases) == 0L) return(keys)
    from <- normalize_marker(names(marker_aliases))
    to   <- normalize_marker(unname(marker_aliases))
    hit  <- match(keys, from)
    ifelse(is.na(hit), keys, to[hit])
  }
  if (!is.null(marker_aliases) && length(marker_aliases) > 0L) {
    renamed <- outgroup_index$marker_key != apply_aliases(outgroup_index$marker_key)
    outgroup_index$marker_key <- apply_aliases(outgroup_index$marker_key)
    if (any(renamed)) {
      message("Outgroup markers renamed onto their ingroup counterpart: ",
              paste(sprintf("%s -> %s", normalize_marker(names(marker_aliases)),
                            normalize_marker(unname(marker_aliases))), collapse = ", "))
    }
  }

  # Readmission is an explicit, named decision, separate from the screening verdict it overrides.
  # A locus enters here because it connects the outgroup, not because it resolves the ingroup.
  if (!is.null(readmit_markers) && length(readmit_markers) > 0L) {
    stop_if_missing_dir(readmit_dir, "readmit_dir")
    pool <- list.files(readmit_dir, pattern = "\\.fasta$", full.names = TRUE)
    pool_index <- dplyr::tibble(path = pool, basename = basename(pool),
                                marker_key = normalize_marker(basename(pool)))
    wanted <- normalize_marker(readmit_markers)
    hit <- pool_index[pool_index$marker_key %in% wanted, , drop = FALSE]
    missing <- setdiff(wanted, hit$marker_key)
    if (length(missing) > 0L) {
      stop("Markers requested for readmission but absent from '", readmit_dir, "': ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    already <- intersect(wanted, ingroup_index$marker_key)
    if (length(already) > 0L) {
      warning("Markers requested for readmission that the screening already retained, ignored: ",
              paste(already, collapse = ", "), call. = FALSE)
      hit <- hit[!hit$marker_key %in% already, , drop = FALSE]
    }
    if (nrow(hit) > 0L) {
      ingroup_index <- dplyr::bind_rows(ingroup_index, hit)
      message("Markers readmitted for outgroup connection, bypassing the Module 3 verdict: ",
              paste(hit$marker_key, collapse = ", "))
    }
  }

  markers_to_process <- sort(unique(ingroup_index$marker_key))

  # The joint marker set is the ingroup's. An outgroup marker with no ingroup counterpart leaves
  # the analysis here, and used to leave without a word: `ndhA` alone carried 42 Portulaca
  # accessions that never reached the supermatrix.
  orphan_keys <- setdiff(unique(outgroup_index$marker_key), markers_to_process)
  if (length(orphan_keys) > 0L) {
    orphan_tbl <- do.call(rbind, lapply(orphan_keys, function(k) {
      f <- outgroup_index$path[outgroup_index$marker_key == k]
      n <- sum(vapply(f, function(x) length(Biostrings::readDNAStringSet(x)), integer(1)))
      data.frame(marker_key = k, files = length(f), n_sequences = n, stringsAsFactors = FALSE)
    }))
    orphan_tbl <- orphan_tbl[order(-orphan_tbl$n_sequences), , drop = FALSE]
    warning("Outgroup markers with no ingroup counterpart, dropped from the joint dataset: ",
            paste(sprintf("%s (%d seqs)", orphan_tbl$marker_key, orphan_tbl$n_sequences),
                  collapse = ", "),
            ". Renaming only helps when the regions are homologous; verify before using ",
            "marker_aliases.", call. = FALSE)
    if (isTRUE(drop_report)) {
      utils::write.csv(orphan_tbl,
                       file.path(out_tables_dir, "TABLE_outgroup_markers_dropped.csv"),
                       row.names = FALSE)
    }
  }

  # Sharing a marker name is not evidence of sharing a region. Both false pairs the 2026-09-01
  # audit found reached this point looking healthy: pepC with 506 aligned columns across the root,
  # trnT-psbD with the most balanced ingroup/outgroup counts in the dataset. Counts cannot see
  # this and neither can a coverage table; only the sequences can. Measured here because this is
  # where the two sets are joined, and reported rather than enforced, because a low share on a
  # fast-evolving locus is a reason to look, not a verdict.
  if (isTRUE(homology_check)) {
    shared_keys <- intersect(markers_to_process, unique(outgroup_index$marker_key))
    if (length(shared_keys) > 0L) {
      message("Checking ingroup/outgroup homology per marker (alignment-free)...")
      hom <- do.call(rbind, lapply(shared_keys, function(kk) {
        ip <- ingroup_index$path[ingroup_index$marker_key == kk]
        op <- outgroup_index$path[outgroup_index$marker_key == kk]
        s20 <- .marker_homology_share(ip, op, k = 20L)
        s10 <- .marker_homology_share(ip, op, k = 10L)
        data.frame(marker_key = kk,
                   share_k20 = round(s20$share, 4),
                   share_k10 = round(s10$share, 4),
                   share_k20_forward = round(s20$share_forward, 4),
                   share_k20_revcomp = round(s20$share_revcomp, 4),
                   n_kmers_outgroup = s20$n_kmers_outgroup,
                   stringsAsFactors = FALSE)
      }))
      hom$verdict <- ifelse(
        is.na(hom$share_k20), "not_evaluated",
        ifelse(hom$share_k20 < homology_k20_min & hom$share_k10 < homology_k10_min,
               "no_detectable_homology",
               # Same region, opposite strand. run_alignment_pipeline(fix_strand = TRUE) puts each
               # marker on one strand within its own pool, which cannot see that the two pools
               # disagree with each other: they are aligned in separate runs.
               ifelse(!is.na(hom$share_k20_revcomp) &
                        hom$share_k20_revcomp > 3 * pmax(hom$share_k20_forward, 1e-9),
                      "reversed_relative_to_ingroup",
                      ifelse(hom$share_k20 < 0.20, "low", "ok"))))
      hom <- hom[order(hom$share_k20), , drop = FALSE]
      utils::write.csv(hom, file.path(out_tables_dir, "TABLE_marker_homology_check.csv"),
                       row.names = FALSE)

      flipped <- hom[hom$verdict == "reversed_relative_to_ingroup", , drop = FALSE]
      if (nrow(flipped) > 0L) {
        warning("Markers whose outgroup sequences are the same region as the ingroup but on the ",
                "opposite strand: ",
                paste(sprintf("%s (forward %.3f, reverse %.3f)", flipped$marker_key,
                              flipped$share_k20_forward, flipped$share_k20_revcomp),
                      collapse = ", "),
                ". run_alignment_pipeline(fix_strand = TRUE) orients each marker within its own ",
                "pool and cannot see a disagreement between the two pools, which are aligned in ",
                "separate runs. Reverse-complement the outgroup file for these markers before ",
                "rerunning Stage 2.", call. = FALSE)
      }

      bad <- hom[hom$verdict == "no_detectable_homology", , drop = FALSE]
      if (nrow(bad) > 0L) {
        warning("Markers whose outgroup sequences share almost no subsequence with their ingroup ",
                "counterpart, so their aligned columns cannot carry positional homology: ",
                paste(sprintf("%s (k20 %.3f, k10 %.3f)", bad$marker_key, bad$share_k20,
                              bad$share_k10), collapse = ", "),
                ". Genuine counterparts in this dataset return 0.28 to 0.60 at k = 20. Such a ",
                "locus places non-homologous characters on the terminals that define the root, ",
                "which is where the outgroup branch is estimated. See ",
                "tables/TABLE_marker_homology_check.csv.", call. = FALSE)
      }
    }
  }

  message("Reading and validating metadata...")
  meta_in <- utils::read.csv(metadata_in_file, stringsAsFactors = FALSE)
  meta_out <- utils::read.csv(metadata_out_file, stringsAsFactors = FALSE)

  assert_cols(meta_in, c("species", "Marker_std", "sid"), "meta_in")
  assert_cols(meta_out, c("species", "Marker_std", "sid"), "meta_out")

  metadata_all <- dplyr::bind_rows(
    dplyr::mutate(meta_in, is_outgroup = FALSE),
    dplyr::mutate(meta_out, is_outgroup = TRUE)
  ) |>
    dplyr::mutate(species = normalize_species(species), marker_key = normalize_marker(as.character(Marker_std)), sid = as.character(sid)) |>
    # The alias has to reach the metadata too, or the renamed alignment finds no accession rows.
    dplyr::mutate(marker_key = ifelse(is_outgroup, apply_aliases(marker_key), marker_key)) |>
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
  write_clean_table(tabla_raw, "TABLE_raw_species_metadata.csv")
  if (nrow(unmatched_metadata) > 0) write_clean_table(unmatched_metadata, "TABLE_unmatched_metadata_rows.csv")

  if (grepl("\\.xlsx?$", accepted_list_file, ignore.case = TRUE)) {
    sheets <- readxl::excel_sheets(accepted_list_file)
    checklist_sheets <- sheets[!grepl("^facts", sheets, ignore.case = TRUE)]
    list_df <- lapply(checklist_sheets, function(sh) {
      df <- readxl::read_excel(accepted_list_file, sheet = sh)
      df$family <- sh
      df
    })
    full_checklist <- dplyr::bind_rows(list_df)
  } else {
    full_checklist <- utils::read.csv(accepted_list_file, stringsAsFactors = FALSE)
    if (!"family" %in% names(full_checklist)) full_checklist$family <- "Cactaceae"
  }

  if (!"pureName" %in% names(full_checklist)) {
    nm <- tolower(names(full_checklist))
    candidates <- intersect(nm, c("scientificname", "species", "name", "taxon", "purename"))
    if (length(candidates) == 0) stop("Checklist missing a scientific name column.", call. = FALSE)
    names(full_checklist)[nm == candidates[1]] <- "pureName"
  }

  clean_checklist_species <- function(df) {
    if ("RANK" %in% names(df)) df <- dplyr::filter(df, RANK == "Species")
    df |>
      dplyr::mutate(pureName = stringr::str_squish(as.character(pureName))) |>
      dplyr::filter(!stringr::str_detect(pureName, stringr::regex("\\b(subsp|ssp|var|forma|f\\.|subg|sect|ser|cf\\.|aff\\.|sp\\.|spp\\.|nr\\.)\\b", ignore_case = TRUE))) |>
      dplyr::mutate(species = normalize_species(pureName)) |>
      dplyr::filter(nzchar(species)) |>
      dplyr::pull(species) |> unique()
  }

  checklist_cactaceae <- clean_checklist_species(dplyr::filter(full_checklist, tolower(family) == "cactaceae"))
  checklist_anacampserotaceae <- clean_checklist_species(dplyr::filter(full_checklist, tolower(family) == "anacampserotaceae"))
  accepted_species <- unique(c(checklist_cactaceae, checklist_anacampserotaceae))

  outgroup_species <- unique(dplyr::pull(dplyr::filter(metadata_all, is_outgroup), species))

  tabla_raw <- tabla_raw |>
    dplyr::mutate(
      matched_to_metadata = !is.na(sid),
      accepted_in_checklist = dplyr::if_else(source_branch == "ingroup", species %in% checklist_cactaceae, species %in% accepted_species),
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
  write_clean_table(tabla_raw, "TABLE_raw_with_acceptance.csv")

  sequence_registry <- tabla_raw |>
    dplyr::select(marker_key, source_branch, species, fasta_name, sid, is_outgroup, matched_to_metadata, accepted_in_checklist, is_outgroup_species, accepted_or_outgroup, matched_to_checklist, record_status, is_duplicate_species_marker, n_records_total, n_records_ingroup, n_records_outgroup, duplicate_class, classification_conflict) |>
    dplyr::arrange(marker_key, species, source_branch, fasta_name)

  message("Exporting TABLE_sequence_registry_with_acceptance.csv ...")
  write_clean_table(sequence_registry, "TABLE_sequence_registry_with_acceptance.csv")

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

  message("Exporting TABLE_species_marker_sid_matrix.csv ...")
  write_clean_table(tabla_wide, "TABLE_species_marker_sid_matrix.csv")

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

  message("Exporting TABLE_duplicate_resolution_species_marker.csv ...")
  write_clean_table(duplicate_resolution, "TABLE_duplicate_resolution_species_marker.csv")

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
    write_clean_table(missing_map, "TABLE_missing_sid_header_map.csv")
  }

  message("Exporting decoupled cleaned FASTAs (Ingroup, Outgroup, Joint)...")
  exported_markers <- character(0)
  markers_without_sequences <- character(0)
  exported_registry_list <- vector("list", length(markers_to_process))
  
  metrics_ingroup_list  <- vector("list", length(markers_to_process))
  metrics_outgroup_list <- vector("list", length(markers_to_process))
  metrics_joint_list    <- vector("list", length(markers_to_process))

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
    
    # Decouple ingroup vs outgroup subsets
    sp_out_names <- vapply(names(aln_out), extract_species, character(1))
    is_out_taxon <- sp_out_names %in% outgroup_species
    
    aln_ingroup  <- aln_out[!is_out_taxon]
    aln_outgroup <- aln_out[is_out_taxon]
    
    # 1. Joint FASTA
    out_file_joint  <- file.path(out_marker_joint_dir, paste0(mk, ".fasta"))
    Biostrings::writeXStringSet(aln_out, out_file_joint)
    
    # 2. Ingroup only FASTA
    if (length(aln_ingroup) > 0) {
      out_file_in <- file.path(out_marker_ingroup_dir, paste0(mk, ".fasta"))
      Biostrings::writeXStringSet(aln_ingroup, out_file_in)
    }
    
    # 3. Outgroup only FASTA
    if (length(aln_outgroup) > 0) {
      out_file_out <- file.path(out_marker_outgroup_dir, paste0(mk, ".fasta"))
      Biostrings::writeXStringSet(aln_outgroup, out_file_out)
    }
    
    # Compute decoupled metrics
    metrics_ingroup_list[[i]]  <- compute_alignment_metrics(aln_ingroup, mk, "ingroup")
    metrics_outgroup_list[[i]] <- compute_alignment_metrics(aln_outgroup, mk, "outgroup")
    metrics_joint_list[[i]]    <- compute_alignment_metrics(aln_out, mk, "joint")
    
    exported_markers <- c(exported_markers, mk)
    message(sprintf("Exported %s: %d total (%d ingroup, %d outgroup)", mk, length(aln_out), length(aln_ingroup), length(aln_outgroup)))
    
    exported_registry_list[[i]] <- dplyr::tibble(
      marker_key = mk,
      fasta_name = names(aln_out),
      species = vapply(names(aln_out), extract_species, character(1)),
      source_branch = dplyr::if_else(names(aln_out) %in% names(aln_outgroup), "outgroup", "ingroup")
    )
  }

  exported_registry <- dplyr::bind_rows(exported_registry_list)
  
  # Export decoupled metrics tables
  table_metrics_ingroup  <- dplyr::bind_rows(metrics_ingroup_list)
  table_metrics_outgroup <- dplyr::bind_rows(metrics_outgroup_list)
  table_metrics_joint    <- dplyr::bind_rows(metrics_joint_list)
  
  write_clean_table(table_metrics_ingroup, "TABLE_marker_metrics_ingroup.csv")
  write_clean_table(table_metrics_outgroup, "TABLE_marker_metrics_outgroup.csv")
  write_clean_table(table_metrics_joint, "TABLE_marker_metrics_joint.csv")

  message("Exporting marker taxon composition table...")
  raw_counts_group <- tabla_raw |> dplyr::group_by(marker_key, source_branch) |> dplyr::summarise(n_species = dplyr::n_distinct(species[!is.na(species)]), n_records = dplyr::n_distinct(fasta_name[!is.na(fasta_name)]), .groups = "drop") |> dplyr::mutate(stage = "raw", group = source_branch) |> dplyr::select(marker_key, stage, group, n_species, n_records)
  raw_counts_total <- tabla_raw |> dplyr::group_by(marker_key) |> dplyr::summarise(n_species = dplyr::n_distinct(species[!is.na(species)]), n_records = dplyr::n_distinct(fasta_name[!is.na(fasta_name)]), .groups = "drop") |> dplyr::mutate(stage = "raw", group = "total") |> dplyr::select(marker_key, stage, group, n_species, n_records)
  selected_counts_group <- tabla_selected_final |> dplyr::mutate(group = dplyr::if_else(is_outgroup_species, "outgroup", "ingroup")) |> dplyr::group_by(marker_key, group) |> dplyr::summarise(n_species = dplyr::n_distinct(species[!is.na(species)]), n_records = dplyr::n_distinct(sid[!is.na(sid)]), .groups = "drop") |> dplyr::mutate(stage = "selected") |> dplyr::select(marker_key, stage, group, n_species, n_records)
  selected_counts_total <- tabla_selected_final |> dplyr::group_by(marker_key) |> dplyr::summarise(n_species = dplyr::n_distinct(species[!is.na(species)]), n_records = dplyr::n_distinct(sid[!is.na(sid)]), .groups = "drop") |> dplyr::mutate(stage = "selected", group = "total") |> dplyr::select(marker_key, stage, group, n_species, n_records)
  exported_counts_group <- exported_registry |> dplyr::group_by(marker_key, source_branch) |> dplyr::summarise(n_species = dplyr::n_distinct(species[!is.na(species)]), n_records = dplyr::n_distinct(fasta_name[!is.na(fasta_name)]), .groups = "drop") |> dplyr::mutate(stage = "exported", group = source_branch) |> dplyr::select(marker_key, stage, group, n_species, n_records)
  exported_counts_total <- exported_registry |> dplyr::group_by(marker_key) |> dplyr::summarise(n_species = dplyr::n_distinct(species[!is.na(species)]), n_records = dplyr::n_distinct(fasta_name[!is.na(fasta_name)]), .groups = "drop") |> dplyr::mutate(stage = "exported", group = "total") |> dplyr::select(marker_key, stage, group, n_species, n_records)

  marker_taxon_composition <- dplyr::bind_rows(raw_counts_group, raw_counts_total, selected_counts_group, selected_counts_total, exported_counts_group, exported_counts_total) |> dplyr::arrange(marker_key, factor(stage, levels = c("raw", "selected", "exported")), group)
  write_clean_table(marker_taxon_composition, "TABLE_marker_taxon_composition.csv")

  message("Exporting enriched marker summary table...")
  raw_summary <- tabla_raw |> dplyr::group_by(marker_key) |> dplyr::summarise(n_records_raw_total = dplyr::n_distinct(fasta_name[!is.na(fasta_name)]), n_records_raw_ingroup = dplyr::n_distinct(fasta_name[source_branch == "ingroup" & !is.na(fasta_name)]), n_records_raw_outgroup = dplyr::n_distinct(fasta_name[source_branch == "outgroup" & !is.na(fasta_name)]), n_species_raw_total = dplyr::n_distinct(species[!is.na(species)]), n_species_raw_ingroup = dplyr::n_distinct(species[source_branch == "ingroup" & !is.na(species)]), n_species_raw_outgroup = dplyr::n_distinct(species[source_branch == "outgroup" & !is.na(species)]), n_species_raw_accepted_ingroup = dplyr::n_distinct(species[source_branch == "ingroup" & accepted_in_checklist & !is.na(species)]), n_species_raw_rejected_ingroup = dplyr::n_distinct(species[source_branch == "ingroup" & !accepted_in_checklist & !is.na(species)]), .groups = "drop")
  duplicate_summary <- duplicate_groups |> dplyr::group_by(marker_key) |> dplyr::summarise(n_species_marker_duplicated = sum(n_records_total > 1, na.rm = TRUE), n_duplicate_records_excess_total = sum(pmax(n_records_total - 1L, 0L), na.rm = TRUE), n_duplicate_records_excess_ingroup = sum(pmax(n_records_ingroup - 1L, 0L), na.rm = TRUE), n_duplicate_records_excess_outgroup = sum(pmax(n_records_outgroup - 1L, 0L), na.rm = TRUE), n_duplicate_records_excess_mixed = sum(dplyr::if_else(n_records_ingroup > 0 & n_records_outgroup > 0, n_records_total - 1L, 0L), na.rm = TRUE), n_classification_conflicts = sum(classification_conflict, na.rm = TRUE), .groups = "drop")
  selected_summary <- tabla_selected_final |> dplyr::group_by(marker_key) |> dplyr::summarise(n_species_selected_total = dplyr::n_distinct(species[!is.na(species)]), n_species_selected_ingroup = dplyr::n_distinct(species[!is.na(species) & !is_outgroup_species]), n_species_selected_outgroup = dplyr::n_distinct(species[!is.na(species) & is_outgroup_species]), n_records_selected_total = dplyr::n_distinct(sid[!is.na(sid)]), .groups = "drop")
  exported_summary <- exported_registry |> dplyr::group_by(marker_key) |> dplyr::summarise(n_sequences_exported = dplyr::n_distinct(fasta_name[!is.na(fasta_name)]), n_species_exported_total = dplyr::n_distinct(species[!is.na(species)]), n_species_exported_ingroup = dplyr::n_distinct(species[!is.na(species) & source_branch == "ingroup"]), n_species_exported_outgroup = dplyr::n_distinct(species[!is.na(species) & source_branch == "outgroup"]), .groups = "drop")

  # Join metrics
  metrics_comparison <- table_metrics_ingroup |>
    dplyr::select(marker, n_seq_ingroup = n_sequences, n_sp_ingroup = n_species, pis_ingroup = parsimony_informative_sites, var_sites_ingroup = variable_sites, gc_pct_ingroup = gc_pct) |>
    dplyr::left_join(
      table_metrics_joint |> dplyr::select(marker, n_seq_joint = n_sequences, n_sp_joint = n_species, pis_joint = parsimony_informative_sites, var_sites_joint = variable_sites, gc_pct_joint = gc_pct, joint_metrics_status = metrics_status),
      by = "marker"
    ) |>
    # joint_metrics_status is carried through so that the blank *_joint cells of
    # TABLE_marker_summary.csv explain themselves. Any marker that mixes ingroup and outgroup
    # sequences is unaligned at this stage (Stage 5 realigns it), so its joint site counts and
    # the inflation ratio are legitimately not computable here.
    dplyr::mutate(
      pis_inflation_ratio = dplyr::if_else(!is.na(pis_ingroup) & pis_ingroup > 0 & !is.na(pis_joint), round(pis_joint / pis_ingroup, 2), NA_real_)
    )

  tabla_summary <- raw_summary |>
    dplyr::left_join(duplicate_summary, by = "marker_key") |>
    dplyr::left_join(selected_summary, by = "marker_key") |>
    dplyr::left_join(exported_summary, by = "marker_key") |>
    dplyr::left_join(metrics_comparison, by = c("marker_key" = "marker")) |>
    dplyr::mutate(
      n_species_marker_duplicated = dplyr::coalesce(n_species_marker_duplicated, 0L),
      n_duplicate_records_excess_total = dplyr::coalesce(n_duplicate_records_excess_total, 0L),
      n_duplicate_records_excess_ingroup = dplyr::coalesce(n_duplicate_records_excess_ingroup, 0L),
      n_duplicate_records_excess_outgroup = dplyr::coalesce(n_duplicate_records_excess_outgroup, 0L),
      n_duplicate_records_excess_mixed = dplyr::coalesce(n_duplicate_records_excess_mixed, 0L),
      n_classification_conflicts = dplyr::coalesce(n_classification_conflicts, 0L),
      n_species_selected_total = dplyr::coalesce(n_species_selected_total, 0L),
      n_species_selected_ingroup = dplyr::coalesce(n_species_selected_ingroup, 0L),
      n_species_selected_outgroup = dplyr::coalesce(n_species_selected_outgroup, 0L),
      n_records_selected_total = dplyr::coalesce(n_records_selected_total, 0L),
      n_sequences_exported = dplyr::coalesce(n_sequences_exported, 0L),
      n_species_exported_total = dplyr::coalesce(n_species_exported_total, 0L),
      n_species_exported_ingroup = dplyr::coalesce(n_species_exported_ingroup, 0L),
      n_species_exported_outgroup = dplyr::coalesce(n_species_exported_outgroup, 0L),
      pct_species_raw_accepted_ingroup = dplyr::if_else(n_species_raw_ingroup > 0, round(100 * n_species_raw_accepted_ingroup / n_species_raw_ingroup, 2), NA_real_),
      pct_species_retained_from_raw = dplyr::if_else(n_species_raw_total > 0, round(100 * n_species_exported_total / n_species_raw_total, 2), NA_real_),
      pct_records_retained_from_raw = dplyr::if_else(n_records_raw_total > 0, round(100 * n_sequences_exported / n_records_raw_total, 2), NA_real_)
    ) |>
    dplyr::arrange(dplyr::desc(n_species_raw_total), marker_key)

  write_clean_table(tabla_summary, "TABLE_marker_summary.csv")

  message("Exporting dataset species summary and species curation status tables...")
  n_checklist_cactaceae <- length(checklist_cactaceae)
  n_checklist_anacampserotaceae <- length(checklist_anacampserotaceae)
  
  n_raw_ingroup_species <- dplyr::n_distinct(tabla_raw$species[tabla_raw$source_branch == "ingroup" & !is.na(tabla_raw$species)])
  n_accepted_ingroup_species <- dplyr::n_distinct(tabla_selected_final$species[!tabla_selected_final$is_outgroup_species & !is.na(tabla_selected_final$species)])
  n_rejected_ingroup_species <- dplyr::n_distinct(tabla_raw$species[tabla_raw$source_branch == "ingroup" & !tabla_raw$accepted_in_checklist & !is.na(tabla_raw$species)])
  
  pct_cactaceae_recovered <- if (n_checklist_cactaceae > 0) round(100 * n_accepted_ingroup_species / n_checklist_cactaceae, 2) else NA_real_

  outgroup_final_species <- unique(tabla_selected_final$species[tabla_selected_final$is_outgroup_species & !is.na(tabla_selected_final$species)])
  # The trailing underscore is required: without it "^Portulaca" also matches Portulacaria
  # (Didiereaceae). Species labels in the registry use underscore as the binomial separator.
  n_selected_anacampserotaceae <- sum(outgroup_final_species %in% checklist_anacampserotaceae | grepl("^(Anacampseros|Grahamia|Talinopsis)_", outgroup_final_species))
  n_selected_portulacaceae <- sum(grepl("^Portulaca_", outgroup_final_species))
  n_selected_talinaceae <- sum(grepl("^(Talinum|Talinella)_", outgroup_final_species))
  n_selected_outgroup_species <- length(outgroup_final_species)
  n_final_joint_species <- dplyr::n_distinct(tabla_selected_final$species[!is.na(tabla_selected_final$species)])

  dataset_species_summary <- tibble::tribble(
    ~metric, ~value, ~details,
    "Total accepted species in Cactaceae checklist (Focal Ingroup)", as.character(n_checklist_cactaceae), "Accepted species in Cactaceae sheet (Caryophyllales.org; Korotkova et al. 2021)",
    "Total accepted species in Anacampserotaceae checklist (Outgroup)", as.character(n_checklist_anacampserotaceae), "Accepted species in Anacampserotaceae sheet (Caryophyllales.org)",
    "Total unique ingroup species with raw sequences (Cactaceae)", as.character(n_raw_ingroup_species), "Distinct Cactaceae binomials retrieved from GenBank ingroup mining",
    "Accepted unique ingroup species retained (Cactaceae)", as.character(n_accepted_ingroup_species), "Cactaceae species validated against checklist and retained in final dataset",
    "Rejected unique ingroup species excluded (Cactaceae)", as.character(n_rejected_ingroup_species), "Synonyms, unresolved names, or invalid binomials excluded from ingroup",
    "Cactaceae focal species recovery rate (%)", paste0(pct_cactaceae_recovered, "%"), "Percentage of accepted Cactaceae checklist species represented in molecular dataset",
    "Unique Anacampserotaceae outgroup species retained", as.character(n_selected_anacampserotaceae), "Provisional (Stage 4) count of Anacampserotaceae retained for rooting; updated in place by run_concatenation_pipeline() (Stage 6)",
    "Unique Portulacaceae outgroup species retained", as.character(n_selected_portulacaceae), "Provisional (Stage 4) count of Portulaca retained for rooting; updated in place by run_concatenation_pipeline() (Stage 6)",
    "Unique Talinaceae outgroup species retained", as.character(n_selected_talinaceae), "Provisional (Stage 4) count of Talinaceae retained for rooting; updated in place by run_concatenation_pipeline() (Stage 6)",
    "Total unique outgroup species retained", as.character(n_selected_outgroup_species), "Provisional (Stage 4) count of outgroups (Anacampserotaceae + Portulacaceae + Talinaceae); updated in place by run_concatenation_pipeline() (Stage 6)",
    "Total unique species in final dataset (Joint)", as.character(n_final_joint_species), "Provisional (Stage 4) species total (Cactaceae + Outgroups); updated in place by run_concatenation_pipeline() (Stage 6)",
    "Total curated sequence records (Cactaceae Ingroup)", as.character(sum(!tabla_selected_final$is_outgroup_species)), "Total single-locus accessions in cleaned_markers_ingroup",
    "Total curated sequence records (Outgroups)", as.character(sum(tabla_selected_final$is_outgroup_species)), "Total single-locus accessions in cleaned_markers_outgroup",
    "Total curated sequence records (Joint)", as.character(nrow(tabla_selected_final)), "Total single-locus accessions in cleaned_markers_joint"
  )

  # --- Loci accounting, added per user request (initial mined loci per group, plus a
  # provisional final-loci count per group). integrate_and_clean_markers() (Stage 4) runs
  # before run_joint_realignment() (Stage 5), which can drop individual taxa from a locus even
  # when the marker file itself survives -- so the "Final loci retained" values below are
  # provisional. run_concatenation_pipeline() (Stage 6) overwrites these same rows in place,
  # in this same file, once the true post-realignment supermatrix exists.
  # Trailing underscore required, as in the species-level split above: species labels in the
  # registry use underscore as the binomial separator, and a bare "^Portulaca" would also match
  # Portulacaria (Didiereaceae), which is neither Portulacaceae nor part of the rooting sample.
  is_anacampserotaceae_sp <- function(sp) sp %in% checklist_anacampserotaceae | grepl("^(Anacampseros|Grahamia|Talinopsis)_", sp)
  is_portulacaceae_sp <- function(sp) grepl("^Portulaca_", sp)
  is_talinaceae_sp <- function(sp) grepl("^(Talinum|Talinella)_", sp)

  outgroup_meta_raw <- dplyr::filter(metadata_all, is_outgroup, !is.na(marker_key), !is.na(species))
  n_initial_loci_cactaceae <- dplyr::n_distinct(metadata_all$marker_key[!metadata_all$is_outgroup & !is.na(metadata_all$marker_key)])
  n_initial_loci_anacampserotaceae <- dplyr::n_distinct(outgroup_meta_raw$marker_key[is_anacampserotaceae_sp(outgroup_meta_raw$species)])
  n_initial_loci_portulacaceae <- dplyr::n_distinct(outgroup_meta_raw$marker_key[is_portulacaceae_sp(outgroup_meta_raw$species)])
  n_initial_loci_talinaceae <- dplyr::n_distinct(outgroup_meta_raw$marker_key[is_talinaceae_sp(outgroup_meta_raw$species)])

  outgroup_exported <- dplyr::filter(exported_registry, source_branch == "outgroup", !is.na(marker_key), !is.na(species))
  n_stage4_loci_cactaceae <- dplyr::n_distinct(exported_registry$marker_key[exported_registry$source_branch == "ingroup" & !is.na(exported_registry$marker_key)])
  n_stage4_loci_anacampserotaceae <- dplyr::n_distinct(outgroup_exported$marker_key[is_anacampserotaceae_sp(outgroup_exported$species)])
  n_stage4_loci_portulacaceae <- dplyr::n_distinct(outgroup_exported$marker_key[is_portulacaceae_sp(outgroup_exported$species)])
  n_stage4_loci_talinaceae <- dplyr::n_distinct(outgroup_exported$marker_key[is_talinaceae_sp(outgroup_exported$species)])
  n_stage4_loci_joint <- length(exported_markers)

  dataset_species_summary <- dplyr::bind_rows(
    dataset_species_summary,
    tibble::tribble(
      ~metric, ~value, ~details,
      "Initial loci mined (Cactaceae Ingroup, Stage 1)", as.character(n_initial_loci_cactaceae), "Distinct standardized markers recovered by assemble_ingroup_phylotar() before saturation screening",
      "Initial loci mined (Anacampserotaceae Outgroup, Stage 1)", as.character(n_initial_loci_anacampserotaceae), "Distinct standardized markers recovered by assemble_outgroup_phylotar() for Anacampseros/Grahamia/Talinopsis",
      "Initial loci mined (Portulacaceae Outgroup, Stage 1)", as.character(n_initial_loci_portulacaceae), "Distinct standardized markers recovered by assemble_outgroup_phylotar() for Portulaca",
      "Initial loci mined (Talinaceae Outgroup, Stage 1)", as.character(n_initial_loci_talinaceae), "Distinct standardized markers recovered by assemble_outgroup_phylotar() for Talinum/Talinella",
      "Final loci retained (Cactaceae Ingroup)", as.character(n_stage4_loci_cactaceae), "Provisional (Stage 3-4) count; updated in place by run_concatenation_pipeline() (Stage 6) once the joint realignment is complete",
      "Final loci retained (Anacampserotaceae Outgroup)", as.character(n_stage4_loci_anacampserotaceae), "Provisional (Stage 3-4) count; updated in place by run_concatenation_pipeline() (Stage 6) once the joint realignment is complete",
      "Final loci retained (Portulacaceae Outgroup)", as.character(n_stage4_loci_portulacaceae), "Provisional (Stage 3-4) count; updated in place by run_concatenation_pipeline() (Stage 6) once the joint realignment is complete",
      "Final loci retained (Talinaceae Outgroup)", as.character(n_stage4_loci_talinaceae), "Provisional (Stage 3-4) count; updated in place by run_concatenation_pipeline() (Stage 6) once the joint realignment is complete",
      "Total unique loci in final joint dataset", as.character(n_stage4_loci_joint), "Provisional (Stage 3-4) count; updated in place by run_concatenation_pipeline() (Stage 6) once the joint realignment is complete"
    )
  )
  write_clean_table(dataset_species_summary, "TABLE_dataset_species_summary.csv")

  species_curation_status <- tabla_raw |>
    dplyr::filter(!is.na(species)) |>
    dplyr::group_by(species) |>
    dplyr::summarise(
      species_clean = dplyr::first(stringr::str_replace_all(species, "_", " ")),
      group = dplyr::case_when(
        any(source_branch == "outgroup") ~ "outgroup",
        TRUE ~ "ingroup"
      ),
      accepted_in_checklist = any(accepted_in_checklist, na.rm = TRUE),
      record_status = dplyr::case_when(
        group == "outgroup" ~ "accepted_outgroup",
        accepted_in_checklist ~ "accepted_ingroup",
        TRUE ~ "rejected_ingroup"
      ),
      retained_in_final_dataset = dplyr::first(species) %in% tabla_selected_final$species,
      n_raw_records = dplyr::n(),
      n_markers_raw = dplyr::n_distinct(marker_key),
      raw_markers = collapse_unique(marker_key),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      tabla_selected_final |>
        dplyr::group_by(species) |>
        dplyr::summarise(
          n_markers_retained = dplyr::n_distinct(marker_key),
          retained_markers = collapse_unique(marker_key),
          .groups = "drop"
        ),
      by = "species"
    ) |>
    dplyr::mutate(
      n_markers_retained = dplyr::coalesce(n_markers_retained, 0L),
      retained_markers = dplyr::coalesce(retained_markers, "none")
    ) |>
    dplyr::arrange(group, dplyr::desc(retained_in_final_dataset), dplyr::desc(n_markers_retained), species)

  write_clean_table(species_curation_status, "TABLE_species_curation_status.csv")

  message("Exporting comprehensive taxonomic backbone molecular matrix...")
  accepted_checklist_species_df <- full_checklist |>
    dplyr::filter(if ("RANK" %in% names(full_checklist)) RANK == "Species" else TRUE) |>
    dplyr::mutate(pureName = stringr::str_squish(as.character(pureName))) |>
    dplyr::filter(!stringr::str_detect(pureName, stringr::regex("\\b(subsp|ssp|var|forma|f\\.|subg|sect|ser|cf\\.|aff\\.|sp\\.|spp\\.|nr\\.)\\b", ignore_case = TRUE))) |>
    dplyr::mutate(species = normalize_species(pureName), species_clean = stringr::str_replace_all(species, "_", " ")) |>
    dplyr::filter(nzchar(species)) |>
    dplyr::distinct(species, .keep_all = TRUE)

  outgroup_df <- dplyr::tibble(
    species = outgroup_species,
    species_clean = stringr::str_replace_all(outgroup_species, "_", " "),
    # Family is assigned from the genus rather than by a Portulacaceae/Anacampserotaceae binary,
    # which silently mislabelled any third outgroup family as Anacampserotaceae.
    family = dplyr::case_when(
      grepl("^Portulaca_", outgroup_species) ~ "Portulacaceae",
      grepl("^(Talinum|Talinella|Amphipetalum)_", outgroup_species) ~ "Talinaceae",
      grepl("^(Anacampseros|Grahamia|Talinopsis)_", outgroup_species) ~ "Anacampserotaceae",
      TRUE ~ NA_character_
    ),
    type = "Outgroup",
    uuid = NA_character_,
    fullName = species_clean,
    pureName = species_clean,
    author = NA_character_,
    RANK = "Species"
  )

  backbone_df <- dplyr::bind_rows(
    accepted_checklist_species_df |> dplyr::mutate(is_outgroup = FALSE),
    outgroup_df |> dplyr::filter(!species %in% accepted_checklist_species_df$species) |> dplyr::mutate(is_outgroup = TRUE)
  ) |>
    dplyr::mutate(
      group = dplyr::if_else(is_outgroup | family %in% c("Portulacaceae", "Anacampserotaceae", "Talinaceae"), "outgroup", "ingroup")
    )

  marker_sids_wide <- tabla_selected_final |>
    dplyr::select(species, marker_key, sid) |>
    dplyr::distinct(species, marker_key, .keep_all = TRUE) |>
    tidyr::pivot_wider(
      names_from = marker_key,
      values_from = sid,
      names_prefix = "sid_"
    )

  taxonomic_backbone_matrix <- backbone_df |>
    dplyr::select(
      family,
      species,
      species_clean,
      fullName,
      uuid,
      author,
      group
    ) |>
    dplyr::left_join(marker_sids_wide, by = "species") |>
    dplyr::left_join(
      tabla_selected_final |>
        dplyr::group_by(species) |>
        dplyr::summarise(
          n_markers_retained = dplyr::n_distinct(marker_key),
          retained_markers = collapse_unique(marker_key),
          .groups = "drop"
        ),
      by = "species"
    ) |>
    dplyr::mutate(
      n_markers_retained = dplyr::coalesce(n_markers_retained, 0L),
      has_molecular_data = n_markers_retained > 0,
      retained_in_curated_dataset = n_markers_retained > 0,
      retained_in_phylogeny = NA,
      iucn_redlist_category = NA_character_
    ) |>
    dplyr::arrange(family, group, dplyr::desc(has_molecular_data), dplyr::desc(n_markers_retained), species)

  write_clean_table(taxonomic_backbone_matrix, "TABLE_taxonomic_backbone_molecular_matrix.csv")

  log_lines <- c(
    "FINAL MARKER INTEGRATION LOG",
    sprintf("Input ingroup directory: %s", ingroup_dir),
    sprintf("Input outgroup directory: %s", outgroup_dir),
    sprintf("Output directory: %s", output_dir),
    sprintf("Markers in ingroup index: %d", nrow(ingroup_index)),
    sprintf("Markers in outgroup index: %d", nrow(outgroup_index)),
    sprintf("Unique ingroup-driven markers processed: %d", length(markers_to_process)),
    sprintf("Total accepted species in Cactaceae checklist: %d", n_checklist_cactaceae),
    sprintf("Total accepted species in Anacampserotaceae checklist: %d", n_checklist_anacampserotaceae),
    sprintf("Unique Cactaceae ingroup species (raw): %d | Accepted: %d | Rejected: %d", n_raw_ingroup_species, n_accepted_ingroup_species, n_rejected_ingroup_species),
    sprintf("Cactaceae checklist recovery rate: %.2f%%", pct_cactaceae_recovered),
    sprintf("Unique outgroup species retained (Anacampserotaceae: %d, Portulacaceae: %d, Total: %d)", n_selected_anacampserotaceae, n_selected_portulacaceae, n_selected_outgroup_species),
    sprintf("Total unique species in final joint dataset: %d", n_final_joint_species),
    sprintf("Rows in raw integration table: %d", nrow(tabla_raw)),
    sprintf("Unmatched metadata rows: %d", nrow(unmatched_metadata)),
    sprintf("Duplicate species \U00d7 marker combinations: %d", nrow(duplicates)),
    sprintf("Selected unique species \U00d7 marker rows (accepted final): %d", nrow(tabla_selected_final)),
    sprintf("Missing SID-to-header mappings: %d", nrow(missing_map)),
    sprintf("Exported FASTA markers: %d", length(exported_markers)),
    sprintf("Decoupled directories populated: cleaned_markers_ingroup, cleaned_markers_outgroup, cleaned_markers_joint"),
    sprintf("Markers with zero combined sequences: %d", length(markers_without_sequences))
  )
  writeLines(log_lines, con = file.path(out_logs_dir, "LOG_clean_integration_summary.txt"))
  
  message("\nMarker integration and cleaning pipeline completed. \U0001f335")
  return(invisible(tabla_summary))
}
