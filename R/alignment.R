#' Infer Positional Homology via MAFFT Alignment
#'
#' Establishes hypotheses of positional homology across unaligned orthologous nucleotide sequence clusters.
#' Positional homology alignment is a crucial prerequisite for maximum-likelihood phylogenetic inference,
#' ensuring that corresponding nucleotide sites derived from common evolutionary ancestry are aligned
#' prior to substitution model evaluation.
#'
#' @param input_fasta Character. Path to unaligned input FASTA file.
#' @param output_fasta Character. Path to destination aligned FASTA output file.
#' @param mafft_exec Character. System command or full path to the executable `MAFFT` binary. Defaults to `"mafft"`.
#' @param mafft_opts Character. Command-line parameters passed directly to `MAFFT`. Defaults to `"--auto"`.
#' @return Invisible numeric exit status code (0 for successful alignment completion).
#' @references
#' Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7:
#' Improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772–780.
#' \doi{10.1093/molbev/mst010}
#' @examples
#' \dontrun{
#' run_mafft(
#'   input_fasta = "raw_cluster.fasta",
#'   output_fasta = "aligned_cluster.fasta",
#'   mafft_opts = "--auto"
#' )
#' }
#' @export
run_mafft <- function(input_fasta, output_fasta, mafft_exec = "mafft", mafft_opts = "--auto") {
  args <- c(strsplit(mafft_opts, "\\s+")[[1]], shQuote(input_fasta))
  stderr_fasta <- paste0(output_fasta, ".stderr.log")
  
  status <- system2(
    mafft_exec,
    args = args,
    stdout = output_fasta,
    stderr = stderr_fasta
  )
  
  if (!file.exists(output_fasta) || file.info(output_fasta)$size == 0) {
    stop("MAFFT did not produce a readable output FASTA. Check logs: ", stderr_fasta, call. = FALSE)
  }
  
  if (!is.null(status) && status != 0) {
    stop(
      sprintf("MAFFT returned non-zero exit status (%s). Check: %s", status, stderr_fasta),
      call. = FALSE
    )
  }
  invisible(status)
}

stop_if_missing_dir <- function(path_dir, label) {
  if (!dir.exists(path_dir)) {
    stop(sprintf("%s does not exist: %s", label, path_dir), call. = FALSE)
  }
}

check_mafft_available <- function(mafft_exec) {
  test <- suppressWarnings(
    system2(mafft_exec, args = "--version", stdout = TRUE, stderr = TRUE)
  )
  status <- attr(test, "status")
  if (!is.null(status) && status != 0) {
    stop("MAFFT is not available or failed to run with '--version'.", call. = FALSE)
  }
  test
}

safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_string_counts <- function(s, pattern) {
  if (length(s) == 0 || is.na(s) || nchar(s) == 0) return(0L)
  sum(strsplit(s, "", fixed = TRUE)[[1]] == pattern)
}

# The ten IUPAC nucleotide ambiguity codes. Kept as a single definition so that the counting
# performed by get_alignment_stats() and the retention performed by clean_ambiguous() can never
# fall out of step.
IUPAC_AMBIGUITY_CODES <- c("R", "Y", "S", "W", "K", "M", "B", "D", "H", "V")

safe_ambiguity_counts <- function(s) {
  if (length(s) == 0 || is.na(s) || nchar(s) == 0) return(0L)
  sum(strsplit(toupper(s), "", fixed = TRUE)[[1]] %in% IUPAC_AMBIGUITY_CODES)
}

clean_ambiguous <- function(x, preserve_iupac = TRUE) {
  seqs <- as.character(x)
  seqs_up <- toupper(seqs)
  # With preserve_iupac = TRUE (the default) the IUPAC ambiguity codes R/Y/S/W/K/M/B/D/H/V are
  # retained. This is the correct default for a phylogenetic pipeline: RAxML-NG and ModelTest-NG
  # incorporate ambiguity into the likelihood as a partial constraint, so an "R" site restricts
  # the state to A or G, whereas an "N" restricts nothing. Collapsing R to N therefore discards
  # real information for no analytical gain, and it also preserves genuine heterozygous signal in
  # multicopy nuclear markers such as nrITS. MAFFT, DECIPHER::MaskAlignment and ape::dist.dna all
  # accept the codes.
  #
  # Setting preserve_iupac = FALSE restores the historical behaviour of collapsing every non-ACGT
  # character to "N", which is available for anyone who explicitly needs an alignment free of
  # ambiguity. Note that doing so makes the gap_handling = "iupac" branch of
  # run_marker_screening() unreachable, because the alignments it reads would already be N-ified
  # here, and the mean_fraction_ambiguous_* columns of the marker summary would report the
  # ambiguity present in the raw input rather than in the exported alignment.
  drop_pattern <- if (isTRUE(preserve_iupac)) "[^ACGTRYSWKMBDHV-]" else "[^ACGT-]"
  seqs_clean <- gsub(drop_pattern, "N", seqs_up, perl = TRUE)
  out <- Biostrings::DNAStringSet(seqs_clean)
  names(out) <- names(x)
  out
}

get_alignment_stats <- function(x) {
  seqs <- as.character(x)
  
  if (length(seqs) == 0) {
    return(data.frame(
      n_sequences = 0L,
      min_sequence_length = NA_integer_,
      max_sequence_length = NA_integer_,
      mean_sequence_length = NA_real_,
      mean_fraction_missing = NA_real_,
      mean_fraction_gaps = NA_real_,
      mean_fraction_ambiguous = NA_real_,
      n_sites_ambiguous = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  seqlen <- nchar(seqs)
  n_missing <- vapply(seqs, safe_string_counts, integer(1), pattern = "N")
  n_gaps <- vapply(seqs, safe_string_counts, integer(1), pattern = "-")
  # IUPAC ambiguity is neither "missing" (N) nor "gap" (-) and was previously invisible to every
  # quality-control table. Counting it here, at each pipeline stage, is what makes a per-locus
  # decision about ambiguity possible: the raw_input stage is measured before clean_ambiguous()
  # runs, so the count reflects what the source records actually contained.
  n_ambig <- vapply(seqs, safe_ambiguity_counts, integer(1))

  pct_missing <- ifelse(seqlen > 0, n_missing / seqlen, NA_real_)
  pct_gaps <- ifelse(seqlen > 0, n_gaps / seqlen, NA_real_)
  pct_ambig <- ifelse(seqlen > 0, n_ambig / seqlen, NA_real_)

  data.frame(
    n_sequences = length(seqs),
    min_sequence_length = min(seqlen),
    max_sequence_length = max(seqlen),
    mean_sequence_length = safe_mean(seqlen),
    mean_fraction_missing = safe_mean(pct_missing),
    mean_fraction_gaps = safe_mean(pct_gaps),
    mean_fraction_ambiguous = safe_mean(pct_ambig),
    n_sites_ambiguous = sum(n_ambig),
    stringsAsFactors = FALSE
  )
}

write_run_log <- function(log_file, input_folder, output_root, mafft_exec, mafft_opts,
                          mask_alignment_regions = NA, min_non_gap_fraction = NA_real_,
                          max_missing_fraction = NA_real_, min_masked_alignment_length = NA_integer_,
                          preserve_iupac = NA) {
  mafft_version <- tryCatch(
    paste(check_mafft_available(mafft_exec), collapse = " "),
    error = function(e) paste("Unavailable:", conditionMessage(e))
  )

  lines <- c(
    "Alignment run log",
    paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("Input folder:", input_folder),
    paste("Output root:", output_root),
    paste("MAFFT executable:", mafft_exec),
    paste("MAFFT options:", mafft_opts),
    paste("MAFFT version:", mafft_version),
    "",
    # Every argument that changes the exported alignments is recorded here, so the call that
    # produced a given 2_MAFFT_* directory can be reconstructed from its own log alone. Without
    # this the only trace of mask_alignment_regions was the masking_applied column of the
    # per-marker manifest, which is easy to miss.
    "Pipeline parameters:",
    paste("  mask_alignment_regions:", mask_alignment_regions),
    paste("  min_non_gap_fraction:", min_non_gap_fraction),
    paste("  max_missing_fraction:", max_missing_fraction),
    paste("  min_masked_alignment_length:", min_masked_alignment_length),
    paste("  preserve_iupac:", preserve_iupac),
    if (isFALSE(mask_alignment_regions)) {
      paste("  NOTE: masking disabled, so min_non_gap_fraction, max_missing_fraction and",
            "min_masked_alignment_length are NOT applied in this stage; occupancy filtering",
            "is deferred to run_joint_realignment().")
    } else {
      NULL
    },
    "",
    "R session info:",
    paste(utils::capture.output(utils::sessionInfo()), collapse = "\n")
  )
  writeLines(lines, con = log_file)
}

build_marker_consolidated_row <- function(marker_id, summary_all, filter_info, final_alignment_written, masking_applied, status) {
  get_stage_row <- function(stage_name) {
    out <- summary_all[summary_all$stage == stage_name, , drop = FALSE]
    if (nrow(out) == 0) {
      return(data.frame(
        stage = stage_name,
        n_sequences = NA_integer_,
        min_sequence_length = NA_integer_,
        max_sequence_length = NA_integer_,
        mean_sequence_length = NA_real_,
        mean_fraction_missing = NA_real_,
        mean_fraction_gaps = NA_real_,
        mean_fraction_ambiguous = NA_real_,
        n_sites_ambiguous = NA_integer_,
        stringsAsFactors = FALSE
      ))
    }
    out[1, , drop = FALSE]
  }
  
  raw_row    <- get_stage_row("raw_input")
  mafft_row  <- get_stage_row("mafft_aligned")
  common_row <- get_stage_row("common_gap_removed")
  masked_row <- get_stage_row("masked_alignment")
  final_row  <- get_stage_row("final_filtered")
  
  n_removed_low_coverage <- if (!is.null(filter_info) && nrow(filter_info) > 0) sum(filter_info$removed_low_coverage, na.rm = TRUE) else 0L
  n_removed_high_missingness <- if (!is.null(filter_info) && nrow(filter_info) > 0) sum(filter_info$removed_high_missingness, na.rm = TRUE) else 0L
  n_removed_total_filter <- if (!is.null(filter_info) && nrow(filter_info) > 0) sum(filter_info$removed_final, na.rm = TRUE) else 0L
  
  raw_n    <- raw_row$n_sequences
  mafft_n  <- mafft_row$n_sequences
  common_n <- common_row$n_sequences
  masked_n <- masked_row$n_sequences
  final_n  <- final_row$n_sequences
  
  raw_len    <- raw_row$max_sequence_length
  mafft_len  <- mafft_row$max_sequence_length
  common_len <- common_row$max_sequence_length
  masked_len <- masked_row$max_sequence_length
  final_len  <- final_row$max_sequence_length
  
  data.frame(
    marker = marker_id,
    
    n_sequences_raw_input = raw_n,
    n_sequences_mafft_aligned = mafft_n,
    n_sequences_common_gap_removed = common_n,
    n_sequences_masked_alignment = masked_n,
    n_sequences_final_filtered = final_n,
    
    alignment_length_raw_input = raw_len,
    alignment_length_mafft_aligned = mafft_len,
    alignment_length_common_gap_removed = common_len,
    alignment_length_masked_alignment = masked_len,
    alignment_length_final_filtered = final_len,
    
    mean_fraction_missing_raw_input = raw_row$mean_fraction_missing,
    mean_fraction_missing_mafft_aligned = mafft_row$mean_fraction_missing,
    mean_fraction_missing_common_gap_removed = common_row$mean_fraction_missing,
    mean_fraction_missing_masked_alignment = masked_row$mean_fraction_missing,
    mean_fraction_missing_final_filtered = final_row$mean_fraction_missing,
    
    mean_fraction_gaps_raw_input = raw_row$mean_fraction_gaps,
    mean_fraction_gaps_mafft_aligned = mafft_row$mean_fraction_gaps,
    mean_fraction_gaps_common_gap_removed = common_row$mean_fraction_gaps,
    mean_fraction_gaps_masked_alignment = masked_row$mean_fraction_gaps,
    mean_fraction_gaps_final_filtered = final_row$mean_fraction_gaps,

    mean_fraction_ambiguous_raw_input = raw_row$mean_fraction_ambiguous,
    mean_fraction_ambiguous_mafft_aligned = mafft_row$mean_fraction_ambiguous,
    mean_fraction_ambiguous_common_gap_removed = common_row$mean_fraction_ambiguous,
    mean_fraction_ambiguous_masked_alignment = masked_row$mean_fraction_ambiguous,
    mean_fraction_ambiguous_final_filtered = final_row$mean_fraction_ambiguous,

    n_sites_ambiguous_raw_input = raw_row$n_sites_ambiguous,
    n_sites_ambiguous_final_filtered = final_row$n_sites_ambiguous,
    
    n_sequences_removed_total = if (!is.na(raw_n) && !is.na(final_n)) raw_n - final_n else NA_integer_,
    n_sequences_removed_low_coverage = n_removed_low_coverage,
    n_sequences_removed_high_missingness = n_removed_high_missingness,
    n_sequences_removed_by_filtering = n_removed_total_filter,
    
    n_sites_removed_by_common_gap_removal = if (!is.na(mafft_len) && !is.na(common_len)) mafft_len - common_len else NA_integer_,
    n_sites_removed_by_masking = if (!is.na(common_len) && !is.na(masked_len)) common_len - masked_len else NA_integer_,
    n_sites_removed_total_from_mafft = if (!is.na(mafft_len) && !is.na(final_len)) mafft_len - final_len else NA_integer_,
    
    pct_sequences_retained_after_mafft = if (!is.na(raw_n) && raw_n > 0) mafft_n / raw_n else NA_real_,
    pct_sequences_retained_after_common_gap_removal = if (!is.na(raw_n) && raw_n > 0) common_n / raw_n else NA_real_,
    pct_sequences_retained_after_masking = if (!is.na(raw_n) && raw_n > 0) masked_n / raw_n else NA_real_,
    pct_sequences_retained_final = if (!is.na(raw_n) && raw_n > 0) final_n / raw_n else NA_real_,
    
    pct_sites_retained_after_common_gap_removal = if (!is.na(mafft_len) && mafft_len > 0) common_len / mafft_len else NA_real_,
    pct_sites_retained_after_masking = if (!is.na(mafft_len) && mafft_len > 0) masked_len / mafft_len else NA_real_,
    pct_sites_retained_final_vs_mafft = if (!is.na(mafft_len) && mafft_len > 0) final_len / mafft_len else NA_real_,
    
    masking_applied = masking_applied,
    final_alignment_written = final_alignment_written,
    status = status,
    stringsAsFactors = FALSE
  )
}

# Parameter stamp appended to every per-marker summary CSV. It is what makes the cache in
# process_marker_file() parameter-aware: reusing an alignment produced under different masking
# or occupancy settings silently invalidates the whole stage, which is exactly how the
# mask_alignment_regions fix could appear to have "no effect" on a re-run over a populated
# 2_MAFFT_* directory.
alignment_param_stamp <- function(mask_alignment_regions, min_non_gap_fraction, max_missing_fraction,
                                  min_masked_alignment_length, preserve_iupac, fix_strand = TRUE) {
  data.frame(
    param_mask_alignment_regions = as.logical(mask_alignment_regions),
    param_min_non_gap_fraction = as.numeric(min_non_gap_fraction),
    param_max_missing_fraction = as.numeric(max_missing_fraction),
    param_min_masked_alignment_length = as.integer(min_masked_alignment_length),
    param_preserve_iupac = as.logical(preserve_iupac),
    # In the stamp so that switching strand correction on or off invalidates the cache. Without
    # this a rerun would silently reuse alignments built under the other setting.
    param_fix_strand = as.logical(fix_strand),
    stringsAsFactors = FALSE
  )
}

cached_params_match <- function(cached_summary, stamp) {
  needed <- names(stamp)
  if (!all(needed %in% names(cached_summary))) return(FALSE)
  if (nrow(cached_summary) == 0) return(FALSE)
  cached <- cached_summary[1, needed, drop = FALSE]
  same <- vapply(needed, function(k) {
    a <- cached[[k]]; b <- stamp[[k]]
    if (is.na(a) && is.na(b)) return(TRUE)
    if (is.na(a) || is.na(b)) return(FALSE)
    if (is.numeric(b)) return(isTRUE(all.equal(as.numeric(a), as.numeric(b))))
    identical(as.character(a), as.character(b))
  }, logical(1))
  all(same)
}

#' Put every sequence of a marker on the same strand before alignment
#'
#' GenBank stores a record on whichever strand the submitter deposited. Two accessions of the same
#' gene can therefore be reverse complements of each other, and MAFFT, which compares only the
#' orientation it is given, will align them anyway: it returns a block, the block looks like a
#' marker, and the reversed rows carry no positional homology to the rest.
#'
#' Found on 2026-09-01 in `Portulaca oleracea` and `P. pilosa`, whose `rbcL` and `matK` accessions
#' are deposited reversed. Their aligned `rbcL` sat at 0.51 observed divergence from Cactaceae
#' where `P. grandiflora`, the same genus, sits at 0.030. `P. oleracea` was at that moment the
#' terminal the acceptance table nominated for rooting the tree.
#'
#' Orientation is decided against the marker's own majority, not against an external reference,
#' because the pipeline mines whatever GenBank holds and no reference is guaranteed. The seed is
#' the longest sequence; every sequence agreeing with the seed extends the reference pool; each
#' remaining sequence is then compared in both directions and flipped when the reverse complement
#' matches better. A sequence matching neither direction is left untouched and reported, since it
#' is a homology problem rather than an orientation one.
#'
#' `mafft --adjustdirection` solves the same problem, and is not used here because it renames the
#' sequences it flips with an `_R_` prefix, which would then have to be undone in the sequence
#' filter log, the name crosswalk and the concatenation. Flipping before MAFFT leaves every
#' downstream stage untouched.
#'
#' @param dna A `DNAStringSet` for one marker.
#' @param k Integer. k-mer length used to compare orientations.
#' @param min_share Numeric. Share of k-mers a sequence must match, in whichever direction, to be
#'   considered part of the marker at all.
#' @param ratio Numeric. How many times better the reverse complement must match before a sequence
#'   is flipped. Above 1 so that a marginal difference never flips anything.
#' @return A list with `dna` (the corrected set) and `log`, one row per sequence.
#' @keywords internal
#' @noRd
.normalise_strand <- function(dna, k = 20L, min_share = 0.15, ratio = 3) {
  n <- length(dna)
  empty_log <- data.frame(Seq = character(0), Length = integer(0), ShareForward = numeric(0),
                          ShareReverse = numeric(0), Action = character(0),
                          stringsAsFactors = FALSE)
  if (n < 3L) return(list(dna = dna, log = empty_log))

  seqs <- gsub("-", "", as.character(dna), fixed = TRUE)
  kset <- function(s) {
    m <- nchar(s) - k + 1L
    if (m < 1L) return(character(0))
    w <- substring(s, seq_len(m), seq_len(m) + k - 1L)
    unique(w[!grepl("[^ACGTacgt]", w)])
  }
  sets <- lapply(toupper(seqs), kset)
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seqs)))
  sets_rc <- lapply(toupper(rc), kset)

  # Seed on the longest sequence, then grow the reference with everything that already agrees with
  # it. Growing first matters: a single seed shares few k-mers with a divergent congener, and the
  # comparison has to be against the marker as a whole.
  usable <- vapply(sets, length, integer(1)) >= 5L
  if (!any(usable)) return(list(dna = dna, log = empty_log))
  seed <- which.max(ifelse(usable, nchar(seqs), -1L))
  pool <- sets[[seed]]
  share <- function(s, p) if (length(s) == 0L) 0 else length(intersect(s, p)) / length(s)
  agree <- vapply(seq_len(n), function(i) usable[i] && share(sets[[i]], pool) >= min_share, logical(1))
  if (sum(agree) > 1L) pool <- unique(unlist(sets[agree], use.names = FALSE))

  fwd <- vapply(seq_len(n), function(i) share(sets[[i]], pool), numeric(1))
  rev <- vapply(seq_len(n), function(i) share(sets_rc[[i]], pool), numeric(1))

  flip <- usable & rev > min_share & rev > ratio * fwd
  action <- ifelse(!usable, "too_short",
            ifelse(flip, "reverse_complemented",
            ifelse(pmax(fwd, rev) < min_share, "no_match_either_direction", "kept")))

  out <- dna
  if (any(flip)) out[flip] <- Biostrings::reverseComplement(dna[flip])

  list(dna = out,
       log = data.frame(Seq = names(dna), Length = nchar(seqs),
                        ShareForward = round(fwd, 4), ShareReverse = round(rev, 4),
                        Action = action, stringsAsFactors = FALSE))
}

process_marker_file <- function(
    fasta_file,
    output_aln_dir,
    output_table_dir,
    mafft_exec,
    mafft_opts,
    mask_alignment_regions = TRUE,
    min_non_gap_fraction = 0.30,
    max_missing_fraction = 0.30,
    min_masked_alignment_length = 100L,
    preserve_iupac = TRUE,
    fix_strand = TRUE
) {
  marker_file <- basename(fasta_file)
  marker_id <- tools::file_path_sans_ext(marker_file)

  clean_fasta_path <- file.path(output_aln_dir, paste0("ALN_clean_input_", marker_id, ".fasta"))
  mafft_fasta_path <- file.path(output_aln_dir, paste0("ALN_mafft_raw_", marker_id, ".fasta"))
  final_fasta_path <- file.path(output_aln_dir, paste0("ALN_masked_final_", marker_id, ".fasta"))

  summary_csv <- file.path(output_table_dir, paste0("TABLE_marker_alignment_summary_", marker_id, ".csv"))
  filter_csv  <- file.path(output_table_dir, paste0("TABLE_marker_sequence_filter_log_", marker_id, ".csv"))

  param_stamp <- alignment_param_stamp(mask_alignment_regions, min_non_gap_fraction,
                                       max_missing_fraction, min_masked_alignment_length,
                                       preserve_iupac, fix_strand)

  if (file.exists(final_fasta_path) && file.exists(summary_csv) && file.exists(filter_csv)) {
    cached_summary <- utils::read.csv(summary_csv, stringsAsFactors = FALSE)
    if (cached_params_match(cached_summary, param_stamp)) {
      message("CACHE: Marker already processed with identical parameters: ", marker_id)
      cached_filter <- utils::read.csv(filter_csv, stringsAsFactors = FALSE)
      # Rebuilt through build_marker_consolidated_row() rather than by returning the last row of
      # the stage table: those two objects have different schemas, and returning the stage row
      # made the do.call(rbind, ...) of run_alignment_pipeline() fail on any partially cached run.
      return(build_marker_consolidated_row(marker_id, cached_summary, cached_filter,
                                           TRUE, isTRUE(mask_alignment_regions), "CACHED_OK"))
    }
    message("CACHE INVALIDATED (alignment parameters changed), reprocessing marker: ", marker_id)
  }

  raw_in <- Biostrings::readDNAStringSet(fasta_file)
  if (length(raw_in) == 0) stop("Input FASTA contains zero sequences.", call. = FALSE)
  if (anyDuplicated(names(raw_in)) > 0) warning(sprintf("Duplicated sequence names detected in %s", marker_file), call. = FALSE)
  
  stats_before <- get_alignment_stats(raw_in)
  
  raw_in_clean <- clean_ambiguous(raw_in, preserve_iupac = preserve_iupac)

  # Before MAFFT, not after: a reversed sequence aligns to nothing, ends with almost no occupancy,
  # and is then dropped by the occupancy filter downstream without ever saying why.
  if (isTRUE(fix_strand)) {
    strand <- .normalise_strand(raw_in_clean)
    if (nrow(strand$log) > 0L) {
      utils::write.csv(strand$log,
                       file.path(output_table_dir, paste0("LOG_STRAND_", marker_id, ".csv")),
                       row.names = FALSE)
      flipped <- strand$log$Seq[strand$log$Action == "reverse_complemented"]
      if (length(flipped) > 0L) {
        message(sprintf("  Marker '%s': %d sequence(s) reverse-complemented onto the marker's strand: %s",
                        marker_id, length(flipped),
                        paste(utils::head(flipped, 4L), collapse = ", ")))
      }
      unmatched <- strand$log$Seq[strand$log$Action == "no_match_either_direction"]
      if (length(unmatched) > 0L) {
        warning(sprintf("Marker '%s': %d sequence(s) match the marker in neither direction: %s. That is a homology problem, not an orientation one; see LOG_STRAND_%s.csv.",
                        marker_id, length(unmatched),
                        paste(utils::head(unmatched, 4L), collapse = ", "), marker_id),
                call. = FALSE)
      }
    }
    raw_in_clean <- strand$dna
  }

  Biostrings::writeXStringSet(raw_in_clean, clean_fasta_path)

  run_mafft(input_fasta = clean_fasta_path, output_fasta = mafft_fasta_path, mafft_exec = mafft_exec, mafft_opts = mafft_opts)
  
  aln <- Biostrings::readDNAStringSet(mafft_fasta_path)
  if (length(aln) == 0) stop("Aligned FASTA contains zero sequences after MAFFT.", call. = FALSE)
  
  stats_after_mafft <- get_alignment_stats(aln)
  
  if (!mask_alignment_regions) {
    # No masking: the MAFFT alignment is exported as-is and the per-sequence occupancy filters
    # (min_non_gap_fraction, max_missing_fraction) and the min_masked_alignment_length floor are
    # deliberately NOT applied here. They are relative to a post-masking column set that does not
    # exist in this branch; occupancy filtering happens once, on the joint ingroup + outgroup
    # alignment, in run_joint_realignment() (Module 5). This is recorded in the parameter stamp
    # and flagged in LOG_alignment_run_info.txt so the omission is never silent.
    Biostrings::writeXStringSet(aln, final_fasta_path)
    summary_all <- rbind(
      cbind(stage = "raw_input", stats_before),
      cbind(stage = "mafft_aligned", stats_after_mafft),
      cbind(stage = "common_gap_removed", stats_after_mafft),
      cbind(stage = "masked_alignment", stats_after_mafft),
      cbind(stage = "final_filtered", stats_after_mafft)
    )
    for (nm in names(param_stamp)) summary_all[[nm]] <- param_stamp[[nm]]
    utils::write.csv(summary_all, summary_csv, row.names = FALSE)

    filter_info <- data.frame(
      sequence_id = names(aln),
      non_gap_sites_after_masking = NA_integer_,
      fraction_missing_after_masking = NA_real_,
      removed_low_coverage = FALSE,
      removed_high_missingness = FALSE,
      removed_final = FALSE,
      stringsAsFactors = FALSE
    )
    utils::write.csv(filter_info, filter_csv, row.names = FALSE)
    
    return(build_marker_consolidated_row(marker_id, summary_all, filter_info, TRUE, FALSE, "OK_NO_MASK"))
  }
  
  aln_no_gaps <- DECIPHER::RemoveGaps(aln, removeGaps = "common")
  stats_after_remove_gaps <- get_alignment_stats(aln_no_gaps)
  
  aln_masked <- methods::as(DECIPHER::MaskAlignment(aln_no_gaps, correction = (length(aln_no_gaps) < 200)), "DNAStringSet")
  stats_after_mask <- get_alignment_stats(aln_masked)
  
  masked_strings <- as.character(aln_masked)
  masked_widths <- nchar(masked_strings)
  
  # Absolute floor on the masked alignment length. min_non_gap_fraction is relative to the
  # POST-masking width, so without this guard an alignment that DECIPHER collapses to a handful
  # of columns still passes: with a 1-column alignment the threshold is 0.3 sites and any
  # sequence carrying one base survives, producing a 1 bp "OK" alignment that propagates to
  # Stage 4. Below min_masked_alignment_length the marker is rejected as ZERO_RETAINED instead.
  if (length(masked_widths) == 0 ||
      max(masked_widths, na.rm = TRUE) < min_masked_alignment_length) {
    remove_gap_filter <- rep(TRUE, length(masked_strings))
    non_gaps <- rep(0L, length(masked_strings))
    pct_missing_masked <- rep(NA_real_, length(masked_strings))
  } else {
    aln_length_masked <- max(masked_widths, na.rm = TRUE)
    non_gaps <- nchar(gsub("-", "", masked_strings))
    min_non_gap_sites <- aln_length_masked * min_non_gap_fraction
    remove_gap_filter <- non_gaps < min_non_gap_sites
    
    pct_missing_masked <- vapply(masked_strings, function(s) {
      if (is.na(s) || nchar(s) == 0) return(NA_real_)
      safe_string_counts(s, "N") / nchar(s)
    }, numeric(1))
  }
  
  remove_missing_filter <- pct_missing_masked > max_missing_fraction
  remove_missing_filter[is.na(remove_missing_filter)] <- TRUE
  remove_final <- remove_gap_filter | remove_missing_filter
  
  aln_final <- aln_masked[!remove_final]
  stats_final <- get_alignment_stats(aln_final)
  final_written <- FALSE
  if (length(aln_final) > 0) {
    Biostrings::writeXStringSet(aln_final, final_fasta_path)
    final_written <- TRUE
  }
  
  filter_info <- data.frame(
    sequence_id = names(aln_masked),
    non_gap_sites_after_masking = non_gaps,
    fraction_missing_after_masking = pct_missing_masked,
    removed_low_coverage = remove_gap_filter,
    removed_high_missingness = remove_missing_filter,
    removed_final = remove_final,
    stringsAsFactors = FALSE
  )
  utils::write.csv(filter_info, filter_csv, row.names = FALSE)
  
  summary_all <- rbind(
    cbind(stage = "raw_input", stats_before),
    cbind(stage = "mafft_aligned", stats_after_mafft),
    cbind(stage = "common_gap_removed", stats_after_remove_gaps),
    cbind(stage = "masked_alignment", stats_after_mask),
    cbind(stage = "final_filtered", stats_final)
  )
  for (nm in names(param_stamp)) summary_all[[nm]] <- param_stamp[[nm]]
  utils::write.csv(summary_all, summary_csv, row.names = FALSE)

  return(build_marker_consolidated_row(marker_id, summary_all, filter_info, final_written, TRUE, if (final_written) "OK" else "ZERO_RETAINED"))
}

#' Execute Complete Alignment and Gap-Masking Pipeline
#'
#' Orchestrates multiple sequence alignment (MSA) and automated quality-control masking across orthologous sequence clusters.
#' Primary alignment hypotheses are inferred using `MAFFT` (Katoh & Standley, 2013). Subsequently, ambiguous sites,
#' poorly aligned terminal fragments, and non-homologous insertions are masked using the `DECIPHER` framework (Wright, 2024),
#' eliminating systematic noise while retaining phylogenetically informative nucleotide positions for downstream supermatrix assembly.
#'
#' @param input_folder Character. Directory path containing raw unaligned orthologous FASTA files.
#' @param output_dir Character. Root destination directory for output subfolders (`alignments/`, `tables/`, `logs/`).
#' @param fasta_pattern Character. Regular expression pattern matching target FASTA files. Defaults to `"\\.fasta$"`.
#' @param mask_alignment_regions Logical. Apply automated alignment masking via `DECIPHER`? Defaults to `TRUE`.
#' @param min_non_gap_fraction Numeric. Minimum allowable proportion of non-gap characters required to retain a site column. Defaults to `0.30`. Applied after masking and relative to the post-masking alignment width.
#' @param max_missing_fraction Numeric. Maximum allowable proportion of missing or ambiguous characters (`N`) allowed per sequence. Defaults to `0.30`.
#' @param min_masked_alignment_length Integer. Absolute minimum number of alignment columns that must survive masking for the locus to be retained. Loci falling below this floor are reported as `ZERO_RETAINED` rather than exported as near-empty alignments. **Defaults to `100L`.** The floor exists to catch alignments that masking has degraded to the point of being uninformative, not to arbitrate between loci of different lengths: no marker in the reference Cactaceae dataset is shorter than 100 columns, so any masked alignment falling below that value is degenerate rather than merely short. Raising the floor can only reject markers, never admit them; a marker rejected by it is reported with `decision_reason` naming the threshold, so the effect is always visible in the screening table.
#'
#'   **Scope.** This floor is evaluated only in the masking branch, that is when `mask_alignment_regions = TRUE`. With masking deferred (`FALSE`, the configuration used for the outgroup) it is deliberately not applied, because it is defined against a post-masking column set that does not exist in that branch. Outgroup terminals are filtered instead by per-sequence occupancy in `run_joint_realignment()` (Module 5). Do not assume this parameter protects both branches.
#' @param preserve_iupac Logical. Retain IUPAC ambiguity codes (`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`) instead of collapsing them to `N` before alignment. **Defaults to `TRUE`.** `RAxML-NG` and `ModelTest-NG` incorporate ambiguity into the likelihood as a partial constraint, so an `R` site restricts the state to A or G whereas an `N` restricts nothing: collapsing the codes discards real information for no analytical gain. Set to `FALSE` only when an alignment free of ambiguity is explicitly required. The amount of ambiguity present at every stage is reported in the `mean_fraction_ambiguous_*` and `n_sites_ambiguous_*` columns of the marker summary, so the decision can be revisited per locus with data. See `@details`.
#' @param fix_strand Logical. Put every sequence of a marker on the same strand before alignment, deciding orientation against the marker's own majority. GenBank stores each record on whichever strand the submitter deposited, and `MAFFT` compares only the orientation it is given: a reverse-complemented accession is aligned anyway and contributes columns with no positional homology. Found on 2026-09-01 in the `rbcL` and `matK` accessions of `Portulaca oleracea` and `P. pilosa`, which sat at 0.51 and 0.40 observed divergence from Cactaceae where `P. grandiflora` sits at 0.030 and 0.066. Each sequence is logged with its match in both directions in `tables/LOG_STRAND_<marker>.csv`. Defaults to `TRUE`.
#' @param mafft_exec Character. System command or full path to the executable `MAFFT` binary. Defaults to `"mafft"`.
#' @param mafft_opts Character. Command-line parameters passed directly to `MAFFT`. Defaults to `"--auto"`.
#' @return A data frame containing site length, missingness, and sequence retention statistics across processed loci.
#' @details
#' **Masking policy.** With `mask_alignment_regions = TRUE` (the default) the `MAFFT` alignment is
#' passed through `DECIPHER::RemoveGaps()` and `DECIPHER::MaskAlignment()`, and the resulting
#' post-masking column set is what `min_non_gap_fraction`, `max_missing_fraction` and
#' `min_masked_alignment_length` are evaluated against. With `mask_alignment_regions = FALSE` the
#' `MAFFT` alignment is exported unmodified and **those three thresholds are not applied at this
#' stage**: they are relative to a masked column set that does not exist in that branch. Occupancy
#' filtering then happens once, on the joint ingroup plus outgroup alignment, in
#' `run_joint_realignment()`. This is the recommended setting for a sparsely sampled outgroup:
#' masking two to seven highly divergent accessions on their own defines a column set the ingroup
#' does not share and can erode an outgroup alignment to a few base pairs, which then fails the
#' occupancy filter of the joint realignment and removes the outgroup precisely from the loci
#' needed to root the tree. Both the chosen policy and the parameter values are written to
#' `logs/LOG_alignment_run_info.txt` and stamped into every per-marker summary table.
#'
#' **Caching.** A marker is reused from a previous run only when its exported alignment, summary
#' table and filter log all exist *and* the five parameters stamped into the cached summary table
#' match the current call. Any parameter change invalidates the cache and the marker is
#' reprocessed, so changing `mask_alignment_regions` on an already populated output directory
#' takes effect instead of silently returning the previous alignments.
#'
#' **Ambiguity policy.** IUPAC ambiguity codes are retained by default (`preserve_iupac = TRUE`).
#' The downstream tools all accept them, and they carry information that `N` does not: an `R`
#' site constrains the state to A or G, an `N` constrains nothing. There is no analytical reason
#' to discard that constraint, so the pipeline does not.
#'
#' Ambiguity is nevertheless measured rather than assumed away. The marker summary reports
#' `mean_fraction_ambiguous_*` for each of the five processing stages, and `n_sites_ambiguous_*`
#' for the raw input and the final alignment. The `raw_input` figures are computed before
#' `clean_ambiguous()` is applied, so they record what the source records actually contained
#' regardless of the policy in force. That is what makes a per-locus decision possible: a marker
#' whose ambiguity is concentrated rather than diffuse can be examined on its own evidence
#' instead of being subjected to a global rule.
#' @references
#' Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7:
#' Improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772–780.
#' \doi{10.1093/molbev/mst010}
#'
#' Wright, E. S. (2024). Fast and Flexible Search for Homologous Biological Sequences with DECIPHER v3.
#' *The R Journal*, 16(2), 191-200. \doi{10.18129/B9.bioc.DECIPHER}
#' @examples
#' \dontrun{
#' run_alignment_pipeline(
#'   input_folder = "1_phylotaR_out_ingroup",
#'   output_dir = "2_MAFFT_Cactaceae",
#'   min_non_gap_fraction = 0.30,
#'   max_missing_fraction = 0.30
#' )
#' }
#' @export
run_alignment_pipeline <- function(
    input_folder,
    output_dir,
    fasta_pattern = "\\.fasta$",
    mask_alignment_regions = TRUE,
    min_non_gap_fraction = 0.30,
    max_missing_fraction = 0.30,
    min_masked_alignment_length = 100L,
    preserve_iupac = TRUE,
    fix_strand = TRUE,
    mafft_exec = "mafft",
    mafft_opts = "--auto"
) {
  stop_if_missing_dir(input_folder, "Input folder")
  
  output_aln_dir <- file.path(output_dir, "alignments")
  output_table_dir <- file.path(output_dir, "tables")
  output_log_dir <- file.path(output_dir, "logs")

  dir.create(output_aln_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_log_dir, recursive = TRUE, showWarnings = FALSE)
  
  check_mafft_available(mafft_exec)
  
  log_file <- file.path(output_log_dir, "LOG_alignment_run_info.txt")
  write_run_log(log_file, input_folder, output_dir, mafft_exec, mafft_opts,
                mask_alignment_regions = mask_alignment_regions,
                min_non_gap_fraction = min_non_gap_fraction,
                max_missing_fraction = max_missing_fraction,
                min_masked_alignment_length = min_masked_alignment_length,
                preserve_iupac = preserve_iupac)


  fasta_files <- sort(list.files(input_folder, pattern = fasta_pattern, full.names = TRUE))
  if (length(fasta_files) == 0) stop("No FASTA files found in input folder.", call. = FALSE)
  
  manifest_list <- vector("list", length(fasta_files))
  for (i in seq_along(fasta_files)) {
    f <- fasta_files[i]
    message("\nProcessing: ", basename(f))
    
    manifest_list[[i]] <- tryCatch(
      process_marker_file(f, output_aln_dir, output_table_dir, mafft_exec, mafft_opts,
                          mask_alignment_regions, min_non_gap_fraction, max_missing_fraction,
                          min_masked_alignment_length, preserve_iupac, fix_strand),
      error = function(e) {
        data.frame(
          marker = tools::file_path_sans_ext(basename(f)),
          n_sequences_raw_input = NA_integer_, n_sequences_mafft_aligned = NA_integer_, n_sequences_common_gap_removed = NA_integer_,
          n_sequences_masked_alignment = NA_integer_, n_sequences_final_filtered = NA_integer_,
          alignment_length_raw_input = NA_integer_, alignment_length_mafft_aligned = NA_integer_, alignment_length_common_gap_removed = NA_integer_,
          alignment_length_masked_alignment = NA_integer_, alignment_length_final_filtered = NA_integer_,
          mean_fraction_missing_raw_input = NA_real_, mean_fraction_missing_mafft_aligned = NA_real_, mean_fraction_missing_common_gap_removed = NA_real_,
          mean_fraction_missing_masked_alignment = NA_real_, mean_fraction_missing_final_filtered = NA_real_,
          mean_fraction_gaps_raw_input = NA_real_, mean_fraction_gaps_mafft_aligned = NA_real_, mean_fraction_gaps_common_gap_removed = NA_real_,
          mean_fraction_gaps_masked_alignment = NA_real_, mean_fraction_gaps_final_filtered = NA_real_,
          n_sequences_removed_total = NA_integer_, n_sequences_removed_low_coverage = NA_integer_, n_sequences_removed_high_missingness = NA_integer_,
          n_sequences_removed_by_filtering = NA_integer_, n_sites_removed_by_common_gap_removal = NA_integer_, n_sites_removed_by_masking = NA_integer_,
          n_sites_removed_total_from_mafft = NA_integer_, pct_sequences_retained_after_mafft = NA_real_, pct_sequences_retained_after_common_gap_removal = NA_real_,
          pct_sequences_retained_after_masking = NA_real_, pct_sequences_retained_final = NA_real_, pct_sites_retained_after_common_gap_removal = NA_real_,
          pct_sites_retained_after_masking = NA_real_, pct_sites_retained_final_vs_mafft = NA_real_,
          masking_applied = mask_alignment_regions, final_alignment_written = FALSE, status = paste("ERROR:", conditionMessage(e)), stringsAsFactors = FALSE
        )
      }
    )
  }
  
  valid_manifests <- Filter(function(x) !is.null(x), manifest_list)
  manifest <- do.call(rbind, valid_manifests)
  
  manifest_path <- file.path(output_log_dir, "LOG_marker_processing_manifest.csv")
  utils::write.csv(manifest, manifest_path, row.names = FALSE)
  
  final_summary_path <- file.path(output_table_dir, "TABLE_alignment_summary_all_markers.csv")
  utils::write.csv(manifest, final_summary_path, row.names = FALSE)
  
  message("\nRun completed. \U0001f335")
  return(manifest)
}

#' Reconcile and Validate Taxonomic Nomenclature
#'
#' Reconciles sequence tip labels against authoritative botanical checklists (e.g., Caryophyllales.org checklist; Korotkova et al., 2021).
#' Resolves taxonomic synonymies, infraspecific variants, and orthographic errors, guaranteeing nomenclatural stability
#' across public GenBank sequence downloads.
#'
#' @param raw_input_fasta Character. Path to input FASTA file containing raw GenBank sequence accessions.
#' @param checklist_path Character. Path to accepted taxonomic checklist CSV or Excel file.
#' @param output_clean_dir Character. Directory path where standardized FASTA sequence output will be saved.
#' @param force_process Logical. Force reprocessing if output cached file exists? Defaults to `FALSE`.
#' @return Character string path to the generated FASTA file with standardized species binomials.
#' @examples
#' \dontrun{
#' clean_taxonomic_names(
#'   raw_input_fasta = "marker_rbcl.fasta",
#'   checklist_path = "CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx",
#'   output_clean_dir = "cleaned_names"
#' )
#' }
#' @export
clean_taxonomic_names <- function(raw_input_fasta, checklist_path, output_clean_dir, force_process = FALSE) {
  
  dir.create(output_clean_dir, recursive = TRUE, showWarnings = FALSE)
  marker_id <- tools::file_path_sans_ext(basename(raw_input_fasta))
  target_out <- file.path(output_clean_dir, paste0(marker_id, "_clean_names.fasta"))
  
  if (!force_process && file.exists(target_out)) {
    message("CACHE: Target clean names FASTA already calculated: ", target_out)
    return(target_out)
  }
  
  dna_in <- Biostrings::readDNAStringSet(raw_input_fasta)
  
  if (grepl("\\.xlsx?$", checklist_path, ignore.case = TRUE)) {
    sheets <- readxl::excel_sheets(checklist_path)
    checklist_sheets <- sheets[!grepl("^facts", sheets, ignore.case = TRUE)]
    list_df <- lapply(checklist_sheets, function(sh) {
      readxl::read_excel(checklist_path, sheet = sh)
    })
    checklist_df <- dplyr::bind_rows(list_df)
  } else {
    checklist_df <- utils::read.csv(checklist_path, stringsAsFactors = FALSE)
  }
  
  accepted_names <- unique(trimws(checklist_df$pureName))
  
  names_seqs <- names(dna_in)
  binomial_leaves <- sub("^([^ ]+_[^ ]+).*", "\\1", names_seqs)
  keep_idx <- binomial_leaves %in% accepted_names
  dna_clean <- dna_in[keep_idx]
  names(dna_clean) <- sub("_", " ", binomial_leaves[keep_idx])
  
  dna_final <- dna_clean[!duplicated(names(dna_clean))]
  Biostrings::writeXStringSet(dna_final, target_out)
  return(target_out)
}

#' Perform Joint Realignment Across Integrated Ingroup and Outgroup Sequences
#'
#' Perform Positional Homology Realignment Across Locus Sequence Alignments
#'
#' Re-estimates positional homology alignments (`MAFFT`) across curated locus FASTA files,
#' performs alignment quality masking with `DECIPHER`, filters low-occupancy sequences,
#' and generates comprehensive alignment statistics.
#'
#' @param input_dir Character. Directory containing curated locus FASTA files (e.g., `4_Cleaned/cleaned_markers_joint` or `4_Cleaned/cleaned_markers_ingroup`).
#' @param output_fasta_dir Character. Directory path to save output realigned FASTA sequence files and QC logs.
#' @param output_aln_dir Character. Directory path to store final masked aligned FASTA files.
#' @param min_non_gap_fraction Numeric. Minimum proportion of non-gap characters, relative to the masked alignment width, required to retain an individual sequence in a locus. Defaults to `0.30`.
#' @param max_missing_fraction Numeric. Maximum proportion of missing characters (`N`) tolerated per sequence. Defaults to `0.30`.
#' @param preserve_iupac Logical. Retain IUPAC ambiguity codes instead of collapsing them to `N` before realignment. Defaults to `TRUE`, matching `run_alignment_pipeline()`; see that function for the rationale.
#' @param rooting_pattern Character or `NULL`. Regular expression identifying the terminals the tree will be rooted on. It exempts nothing. Any matching terminal removed by the occupancy filters is named in a warning and flagged in `LOG_SEQ_FILTER_<marker>.csv`. A filter that deletes the rooting outgroup must say so at the moment it does it, not four modules downstream: the five *Portulaca* sequences of `trnL_trnF`, 297 bp against a threshold of about 353, were removed silently and the branch subtending the outgroup collapsed from roughly 600 expected substitutions to 0.03. Defaults to `NULL`.
#' @param protect_pattern Character or `NULL`. Regular expression matched against sequence names; matching terminals that carry at least one non-gap character are retained regardless of the occupancy filters. Intended as a last-resort safeguard for rooting terminals whose sequences are legitimately short (for example `"^(Anacampseros|Grahamia|Talinopsis|Portulaca)_"`). Defaults to `NULL` (no exemption).
#' @param protect_markers Character vector or `NULL`. Names of the markers, as they appear in `input_dir` without the file extension, in which `protect_pattern` is honoured. `NULL`, the default, applies the exemption to every marker. Naming markers restricts it to the loci where a short outgroup sequence is worth its gap cost, instead of retaining every fragment of every terminal across the whole matrix: in the August 2026 dataset only `trnL_trnF` lost rooting terminals, and only there does the exemption buy anything.
#' @return A data frame containing compiled alignment summary statistics across all processed markers.
#' @details
#' This is the step that removes individual taxa from a locus without removing the locus itself.
#' Sequences that entered Stage 4 as short fragments, typically outgroup accessions that were
#' masked separately from the ingroup in `run_alignment_pipeline()`, fail the occupancy filter here
#' and disappear from the supermatrix. Passing `mask_alignment_regions = FALSE` for the outgroup in
#' Stage 2 so that masking happens only once, jointly, at this stage, is preferable to relaxing
#' these thresholds or resorting to `protect_pattern`, because retaining very short sequences
#' inflates the gap fraction of the final supermatrix and can destabilise the affected terminals.
#' @references
#' Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7:
#' Improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772–780.
#' \doi{10.1093/molbev/mst010}
#' @examples
#' \dontrun{
#' run_joint_realignment(
#'   input_dir = "4_Cleaned/cleaned_markers_joint",
#'   output_fasta_dir = "5_MAFFT_Cleaned",
#'   output_aln_dir = "5_MAFFT_Cleaned/aligned_markers"
#' )
#' }
#' @export
run_joint_realignment <- function(input_dir,
                                  output_fasta_dir,
                                  output_aln_dir,
                                  min_non_gap_fraction = 0.30,
                                  max_missing_fraction = 0.30,
                                  preserve_iupac = TRUE,
                                  protect_pattern = NULL,
                                  protect_markers = NULL,
                                  rooting_pattern = NULL) {
  dir.create(output_fasta_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_aln_dir, recursive = TRUE, showWarnings = FALSE)
  
  log_file <- file.path(output_fasta_dir, "LOG_ALN_FINAL_RUN.txt")
  writeLines(paste0("Run started: ", Sys.time()), log_file)
  
  log_write <- function(text) {
    con <- file(log_file, open = "a")
    writeLines(text, con)
    close(con)
  }
  
  files <- sort(list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE))
  log_write(paste("Markers detected:", length(files)))
  log_write(paste("Input dir:", input_dir))
  log_write(paste("min_non_gap_fraction:", min_non_gap_fraction))
  log_write(paste("max_missing_fraction:", max_missing_fraction))
  log_write(paste("preserve_iupac:", preserve_iupac))
  log_write(paste("protect_pattern:", if (is.null(protect_pattern)) "NULL" else protect_pattern))
  log_write(paste("protect_markers:", if (is.null(protect_markers)) "ALL" else paste(protect_markers, collapse = ", ")))
  log_write(paste("rooting_pattern:", if (is.null(rooting_pattern)) "NULL" else rooting_pattern))

  collapsed_markers <- character()
  summary_list <- list()
  
  for (f in files) {
    fn <- basename(f)
    marker_name <- tools::file_path_sans_ext(fn)
    message("Processing realignment: ", fn)
    
    raw <- Biostrings::readDNAStringSet(f)
    raw_clean <- clean_ambiguous(raw, preserve_iupac = preserve_iupac)
    
    temp_clean <- file.path(output_fasta_dir, paste0("TEMP_", fn))
    raw_out <- file.path(output_fasta_dir, paste0("TEMP_ALN_", fn))
    raw_err <- file.path(output_fasta_dir, paste0("TEMP_ALN_", marker_name, ".log"))
    temp_files <- c(temp_clean, raw_out, raw_err)
    
    tryCatch({
      Biostrings::writeXStringSet(raw_clean, temp_clean)
      
      status <- system2(
        command = "mafft",
        args = c("--auto", temp_clean),
        stdout = raw_out,
        stderr = raw_err
      )
      
      if (!identical(status, 0L)) {
        err_msg <- if (file.exists(raw_err)) paste(readLines(raw_err, warn = FALSE), collapse = "\n") else "No stderr captured."
        stop("MAFFT failed for: ", fn, "\n", err_msg)
      }
      
      aln <- Biostrings::readDNAStringSet(raw_out)
      aln <- DECIPHER::RemoveGaps(aln, removeGaps = "common")
      aln <- DECIPHER::MaskAlignment(aln, correction = (length(aln) < 200))
      aln <- methods::as(aln, "DNAStringSet")
      
      aln_char <- as.character(aln)
      aln_len <- nchar(aln_char[1])
      
      non_gaps <- nchar(gsub("-", "", aln_char))
      pctMissing <- vapply(aln_char, function(s) sum(strsplit(s, "")[[1]] == "N") / nchar(s), numeric(1))
      
      keep <- non_gaps >= (min_non_gap_fraction * aln_len) & pctMissing <= max_missing_fraction
      # The exemption is scoped by marker. Applied matrix-wide it retains every short fragment of
      # every protected terminal in every locus, which is the gap inflation the parameter
      # documentation warns about; scoped to the loci that actually lose rooting terminals it
      # buys the separation without that cost.
      protect_here <- !is.null(protect_pattern) &&
        (is.null(protect_markers) || marker_name %in% protect_markers)
      if (protect_here) {
        protected <- grepl(protect_pattern, names(aln_char)) & non_gaps > 0L
        keep <- keep | protected
        if (any(protected & non_gaps < (min_non_gap_fraction * aln_len))) {
          message("Marker '", marker_name, "': protect_pattern retained ",
                  sum(protected & non_gaps < (min_non_gap_fraction * aln_len)),
                  " terminal(s) below the occupancy threshold.")
        }
      }
      # Reporting only. A rooting terminal lost here is not an error, it is a decision, but it
      # has to be a visible one: these are the terminals the root is placed on and the ones the
      # deepest calibrations are addressed by.
      is_rooting <- if (is.null(rooting_pattern)) rep(FALSE, length(aln_char)) else grepl(rooting_pattern, names(aln_char))
      dropped_rooting <- names(aln_char)[is_rooting & !keep]
      if (length(dropped_rooting) > 0L) {
        warning("Marker '", marker_name, "': the occupancy filter removed ", length(dropped_rooting),
                " rooting terminal(s): ", paste(dropped_rooting, collapse = ", "),
                ". Their non-gap lengths are in LOG_SEQ_FILTER_", marker_name,
                ".csv. Consider protect_pattern for this marker, weighing it against the gap ",
                "fraction those terminals add to the supermatrix.", call. = FALSE)
      }

      aln_final <- aln[keep]
      
      if (length(aln_final) == 0) {
        collapsed_markers <- c(collapsed_markers, fn)
      }
      
      if (length(aln_final) > 0) {
        out_fasta <- file.path(output_aln_dir, fn)
        Biostrings::writeXStringSet(aln_final, out_fasta)
      }
      
      # Quality control exports (always export, even if collapsed)
      utils::write.csv(data.frame(
        Seq = names(aln),
        NonGaps = non_gaps,
        pctMissing = pctMissing,
        Retained = keep,
        IsRooting = is_rooting,
        DroppedRooting = is_rooting & !keep,
        ProtectedHere = if (protect_here) grepl(protect_pattern, names(aln_char)) & non_gaps > 0L else rep(FALSE, length(aln_char))
      ), file = file.path(output_fasta_dir, paste0("LOG_SEQ_FILTER_", marker_name, ".csv")),
      row.names = FALSE)
      
      # Calculate stats
      seqs <- as.character(aln_final)
      if (length(seqs) == 0) {
        stats_df <- data.frame(
          marker = marker_name,
          nSeq = 0L,
          minLen = NA_real_,
          maxLen = NA_real_,
          avgLen = NA_real_,
          pctMissing = NA_real_,
          pctGaps = NA_real_
        )
      } else {
        seqlen <- nchar(seqs)
        n_missing <- vapply(strsplit(seqs, "", fixed = TRUE), function(z) sum(z == "N"), numeric(1))
        n_gaps <- vapply(strsplit(seqs, "", fixed = TRUE), function(z) sum(z == "-"), numeric(1))
        stats_df <- data.frame(
          marker = marker_name,
          nSeq = length(seqs),
          minLen = min(seqlen),
          maxLen = max(seqlen),
          avgLen = round(mean(seqlen), 1),
          pctMissing = round(mean(n_missing / seqlen) * 100, 3),
          pctGaps = round(mean(n_gaps / seqlen) * 100, 3)
        )
      }
      
      utils::write.csv(stats_df,
        file = file.path(output_fasta_dir, paste0("TABLE_ALN_SUMMARY_", marker_name, ".csv")),
        row.names = FALSE)
      
      summary_list[[marker_name]] <- stats_df
    }, finally = {
      existing <- temp_files[file.exists(temp_files)]
      if (length(existing) > 0) file.remove(existing)
    })
  }
  
  compiled_summary <- dplyr::bind_rows(summary_list)
  utils::write.csv(
    compiled_summary,
    file = file.path(output_fasta_dir, "TABLE_alignment_summary_all_markers.csv"),
    row.names = FALSE
  )
  
  log_write(paste("Collapsed markers:", paste(collapsed_markers, collapse = ", ")))
  log_write(paste(utils::capture.output(utils::sessionInfo()), collapse = "\n"))
  log_write(paste("Run finished:", Sys.time()))
  
  message("\nJoint realignment completed. \U0001f335")
  return(invisible(compiled_summary))
}



