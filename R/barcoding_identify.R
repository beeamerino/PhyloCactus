# Identification of a user query (11_barcoding/10_identify/).
# Phase 6 bis of PhyloCactus 0.5.0, decisions I1 to I9 of BMM (27-09): IdTaxa at 60, one model per
# locus trained on the whole library, a crop by 20-mer hits for contigs and plastomes, three reasons
# for state 3, one row per query and locus and never a row that combines loci (E8).

#' The IdTaxa model of one locus trained on its whole library
#'
#' The object CN2 has trained since 6B, shared with the identification step so that both answer a
#' query the same way (decision I8). The training seed keeps the key `"cn2"` it had in 6B, so the
#' CN2 table does not change.
#' @noRd
.bc_idtaxa_whole_locus <- function(lib_dna, lib_tab, locus, seed, trained = NULL) {
  species <- stats::setNames(lib_tab$species, lib_tab$sid)
  lib_m <- as.matrix(lib_dna)[lib_tab$sid, , drop = FALSE]
  pool <- .bc_strand_pool(lib_m)
  lib_seqs <- gsub("-", "", toupper(apply(as.character(lib_m), 1, paste, collapse = "")), fixed = TRUE)
  names(lib_seqs) <- rownames(lib_m)
  if (is.null(trained)) {
    trained <- .bc_idtaxa_train(lib_seqs, species[names(lib_seqs)],
                                train_seed = .bc_idtaxa_whole_locus_seed(seed, locus))
  }
  list(species = species[names(lib_seqs)], pool = pool, lib_seqs = lib_seqs, trained = trained,
       core_pool = .bc_core_pool(lib_m), lib_m = lib_m, genus_cores = new.env())
}

#' The 20-mers of the core of a locus, with their distance to each end of the core (amendment J1)
#'
#' The core is the span of alignment columns covered by `min_cover` or more of the library
#' sequences. Each 20-mer of a library sequence inside that span is kept with the number of bases
#' of that sequence before it and after it inside the core, the largest over the sequences that
#' hold it. The query is then cut only where it runs beyond those distances, so the ends of a
#' divergent query that lies inside the core are kept.
#' @noRd
.bc_core_pool <- function(lib_m, min_cover = 0.5, k = 20L) {
  m <- as.character(lib_m)
  covered <- which(colMeans(m != "-") >= min_cover)
  empty <- data.frame(kmer = character(0), before = integer(0), after = integer(0))
  if (length(covered) == 0L) return(empty)
  span <- m[, min(covered):max(covered), drop = FALSE]
  seqs <- gsub("-", "", toupper(apply(span, 1, paste, collapse = "")), fixed = TRUE)
  parts <- lapply(seqs, function(x) {
    n <- nchar(x) - k + 1L
    if (n < 1L) return(NULL)
    w <- substring(x, seq_len(n), seq_len(n) + k - 1L)
    ok <- !grepl("[^ACGT]", w)
    data.frame(kmer = w[ok], before = (seq_len(n) - 1L)[ok], after = (nchar(x) - (seq_len(n) + k - 1L))[ok],
               stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, parts)
  if (is.null(d)) return(empty)
  before <- tapply(d$before, d$kmer, max)
  after <- tapply(d$after, d$kmer, max)
  data.frame(kmer = names(before), before = as.integer(before), after = as.integer(after[names(before)]),
             stringsAsFactors = FALSE)
}

#' Where a query runs beyond a core: first and last base to keep, or NULL with no core 20-mer
#' @noRd
.bc_core_cut <- function(oriented, pool) {
  hits <- .bc_kmer_hits(oriented, pool$kmer)
  if (length(hits) == 0L) return(NULL)
  h1 <- min(hits)
  h2 <- max(hits)
  k1 <- match(substr(oriented, h1, h1 + 19L), pool$kmer)
  k2 <- match(substr(oriented, h2, h2 + 19L), pool$kmer)
  c(h1 - min(h1 - 1L, pool$before[k1]),
    h2 + 19L + min(nchar(oriented) - (h2 + 19L), pool$after[k2]))
}

#' The core pool of one genus of a locus (amendment J2), computed once per call; NULL for a genus
#' with fewer than `min_sequences` sequences in the locus
#' @noRd
.bc_genus_core <- function(wl, genus, min_sequences = 3L) {
  if (exists(genus, envir = wl$genus_cores, inherits = FALSE)) return(get(genus, envir = wl$genus_cores))
  rows <- names(wl$species)[.bc_genus(unname(wl$species)) == genus]
  pool <- if (length(rows) >= min_sequences) .bc_core_pool(wl$lib_m[rows, , drop = FALSE]) else NULL
  assign(genus, pool, envir = wl$genus_cores)
  pool
}

#' @noRd
.bc_idtaxa_whole_locus_seed <- function(seed, locus) .bc_query_seed(seed, locus, "cn2", 0L, "<training>")

#' The model of a locus from the cache of step 10, or trained and cached
#'
#' The cache holds the training with the md5 of its library file and the seed. A different md5 or
#' seed trains again; otherwise the file is read and not rewritten.
#' @noRd
.bc_identify_model <- function(library_dir, lib, locus, seed, models_dir) {
  lib_file <- file.path(library_dir, paste0("LIB_", locus, ".fasta"))
  md5 <- unname(tools::md5sum(lib_file))
  dna <- ape::read.dna(lib_file, format = "fasta")
  rownames(dna) <- .bc_parse_header(labels(dna))$sid
  lib_tab <- lib[lib$locus == locus, , drop = FALSE]
  cache <- file.path(models_dir, paste0("MODEL_idtaxa_", locus, ".rds"))
  cached <- if (file.exists(cache)) readRDS(cache) else NULL
  if (!is.null(cached) && identical(cached$library_md5, md5) && identical(cached$seed, seed)) {
    return(.bc_idtaxa_whole_locus(dna, lib_tab, locus, seed, trained = cached$trained))
  }
  wl <- .bc_idtaxa_whole_locus(dna, lib_tab, locus, seed)
  dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(trained = wl$trained, locus = locus, library_md5 = md5, seed = seed,
               r_version = R.version.string,
               decipher_version = as.character(utils::packageVersion("DECIPHER"))), cache)
  wl
}

#' Start positions of the k-mers of `s` that are in the pool
#' @noRd
.bc_kmer_hits <- function(s, pool, k = 20L) {
  m <- nchar(s) - k + 1L
  if (m < 1L || length(pool) == 0L) return(integer(0))
  w <- substring(s, seq_len(m), seq_len(m) + k - 1L)
  which(w %in% pool)
}

#' The window of width `width` with the most hits, and the other windows with half as many
#'
#' Windows start at a hit. The region runs from the first to the last hit of the best window. The
#' other windows are found the same way on the hits left outside that region, one after the other,
#' while they hold at least half the hits of the best one.
#' @noRd
.bc_best_window <- function(hits, width, k = 20L) {
  if (length(hits) == 0L) return(NULL)
  pick <- function(h) {
    n <- vapply(h, function(p) sum(h >= p & h <= p + width - k), integer(1))
    i <- which.max(n)
    inside <- h[h >= h[i] & h <= h[i] + width - k]
    list(start = min(inside), end = max(inside) + k - 1L, hits = n[i])
  }
  best <- pick(hits)
  rest <- hits[hits < best$start | hits > best$end - k + 1L]
  other <- 0L
  while (length(rest) > 0L) {
    w <- pick(rest)
    if (w$hits < best$hits / 2) break
    other <- other + 1L
    rest <- rest[rest < w$start | rest > w$end - k + 1L]
  }
  best$other <- other
  best
}

#' Answer of one query for one locus (decisions I3 and I4)
#' @noRd
.bc_identify_one <- function(s, id, locus, wl, width, threshold, min_overlap, seed) {
  n <- nchar(s)
  empty <- function(path, reason, start = NA_integer_, end = NA_integer_, other = NA_integer_,
                    orientation = NA_character_, trimmed = NA_integer_) {
    genus_trimmed <- NA_integer_
    cbind(data.frame(path = path, region_start = start, region_end = end,
                     region_length = if (is.na(start)) NA_integer_ else as.integer(end - start + 1L),
                     other_windows = other, orientation = orientation, core_trimmed = trimmed,
                     genus_core_trimmed = genus_trimmed, stringsAsFactors = FALSE),
          .bc_prediction_row(3L, reason = reason),
          data.frame(genus_idtaxa = NA_character_, species_idtaxa = NA_character_,
                     genus_confidence = NA_real_, species_confidence = NA_real_,
                     stringsAsFactors = FALSE))
  }
  path <- "whole"
  start <- 1L
  end <- n
  other <- 0L
  o <- .bc_orient_to_library(.bc_dnabin_row(s, "query"), wl$pool)
  if (o$orientation == "no_match") {
    if (n <= width) return(empty("none", "no_overlap"))
    fwd <- .bc_best_window(.bc_kmer_hits(s, wl$pool), width)
    rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(s)))
    rev <- .bc_best_window(.bc_kmer_hits(rc, wl$pool), width)
    if (is.null(fwd) && is.null(rev)) return(empty("none", "no_overlap"))
    if (is.null(fwd) || (!is.null(rev) && rev$hits > fwd$hits)) {
      start <- n - rev$end + 1L
      end <- n - rev$start + 1L
      other <- rev$other
    } else {
      start <- fwd$start
      end <- fwd$end
      other <- fwd$other
    }
    path <- "crop"
    region <- substr(s, start, end)
    o <- .bc_orient_to_library(.bc_dnabin_row(region, "query"), wl$pool)
    if (o$orientation == "no_match") return(empty(path, "no_overlap", start, end, other))
  }
  # J1: the region is cut where it runs beyond the core of the locus, in the coordinates of the
  # query. Before the first core 20-mer the query keeps at most as many bases as the core holds
  # before that 20-mer, and the same after the last one.
  oriented <- toupper(paste(as.character(as.matrix(o$query)), collapse = ""))
  # Keeps bases t[1] to t[2] of the oriented region and moves start and end in query coordinates
  cut <- function(t) {
    if (o$orientation == "reverse") {
      new_start <- end - t[2] + 1L
      end <<- end - t[1] + 1L
    } else {
      new_start <- start + t[1] - 1L
      end <<- start + t[2] - 1L
    }
    start <<- new_start
    removed <- as.integer((t[1] - 1L) + (nchar(oriented) - t[2]))
    oriented <<- substr(oriented, t[1], t[2])
    removed
  }
  t <- .bc_core_cut(oriented, wl$core_pool)
  if (is.null(t)) return(empty(path, "no_overlap", start, end, other, o$orientation))
  before_j1 <- list(start = start, end = end, oriented = oriented)
  trimmed <- cut(t)
  if (end - start + 1L < min_overlap) {
    return(empty(path, "short_overlap", start, end, other, o$orientation, trimmed))
  }
  classify <- function() {
    .bc_classify_idtaxa(wl$lib_seqs, wl$species, oriented, threshold = threshold,
                        trained = wl$trained,
                        query_seed = .bc_query_seed(seed, locus, "cn2", 0L, id))
  }
  r <- classify()
  # J3: the region as it was before J1, cut to the core of the genus reached; classified again when
  # that is not the region already classified. The core of a genus can be narrower or wider than
  # that of the locus.
  genus_trimmed <- 0L
  # J3b: only after a first pass that assigns the genus (state 1 or 2)
  if (!is.na(r$genus_idtaxa) && r$state %in% 1:2) {
    gpool <- .bc_genus_core(wl, r$genus_idtaxa)
    tg <- if (is.null(gpool)) NULL else .bc_core_cut(before_j1$oriented, gpool)
    if (!is.null(tg) && tg[2] - tg[1] + 1L >= min_overlap) {
      j1 <- list(start = start, end = end, oriented = oriented)
      start <- before_j1$start
      end <- before_j1$end
      oriented <- before_j1$oriented
      genus_trimmed <- cut(tg)
      if (start == j1$start && end == j1$end) {
        genus_trimmed <- 0L
      } else {
        r <- classify()
      }
    }
  }
  cbind(data.frame(path = path, region_start = start, region_end = end,
                   region_length = as.integer(end - start + 1L), other_windows = other,
                   orientation = o$orientation, core_trimmed = trimmed,
                   genus_core_trimmed = genus_trimmed, stringsAsFactors = FALSE), r)
}

#' The wrong-species rates of 6D for each locus, or NA with the reason (decision I6)
#' @noRd
.bc_identify_context <- function(metrics_dir, loci, threshold) {
  na <- data.frame(locus = loci, validation_ws_rate_species_present = NA_real_,
                   validation_ws_rate_species_absent = NA_real_, stringsAsFactors = FALSE)
  if (is.null(metrics_dir)) {
    message("No metrics_dir: the validation columns are NA.")
    return(na)
  }
  f <- file.path(metrics_dir, "TABLE_barcoding_metrics_summary_idtaxa.csv")
  if (!file.exists(f)) {
    message("No 6D metrics table (", f, "): the validation columns are NA.")
    return(na)
  }
  m <- utils::read.csv(f, stringsAsFactors = FALSE)
  m <- m[m$method == "idtaxa" & m$locus %in% loci, , drop = FALSE]
  if (nrow(m) > 0L && any(m$threshold != threshold)) {
    message("The 6D rates were measured at threshold ", paste(unique(m$threshold), collapse = ", "),
            " and this call uses threshold ", threshold, ": the validation columns are NA.")
    return(na)
  }
  get <- function(l, sc) {
    v <- m$wrong_species_rate[m$locus == l & m$scheme == sc]
    if (length(v) == 1L) v else NA_real_
  }
  na$validation_ws_rate_species_present <- vapply(loci, get, numeric(1), sc = "species")
  na$validation_ws_rate_species_absent <- vapply(loci, get, numeric(1), sc = "genus")
  na
}

#' The query as a named character vector, from a FASTA path, a DNAStringSet or a character vector
#' @noRd
.bc_identify_read_query <- function(query) {
  if (is.character(query) && length(query) == 1L && is.null(names(query)) && file.exists(query)) {
    first <- readLines(query, n = 1L, warn = FALSE)
    if (grepl("\\.(fastq|fq)(\\.gz)?$", query, ignore.case = TRUE) ||
        (length(first) == 1L && startsWith(first, "@"))) {
      stop("This is a FASTQ file of raw reads. Raw reads need assembly or mapping first, which is ",
           "Phase 8 of the branch; step 10 takes sequences (FASTA, a DNAStringSet or a character ",
           "vector).", call. = FALSE)
    }
    if (grepl("\\.(gb|gbk|gbff|genbank)$", query, ignore.case = TRUE) ||
        (length(first) == 1L && startsWith(first, "LOCUS"))) {
      pl <- .bc_parse_genbank(readLines(query, warn = FALSE))
      query <- stats::setNames(pl$sequence, pl$accession)
    } else {
      query <- Biostrings::readDNAStringSet(query)
    }
  }
  if (methods::is(query, "DNAStringSet")) query <- stats::setNames(as.character(query), names(query))
  if (!is.character(query) || length(query) == 0L) {
    stop("query must be a FASTA path, a DNAStringSet or a character vector of sequences.", call. = FALSE)
  }
  if (is.null(names(query))) names(query) <- paste0("query_", seq_along(query))
  empty_names <- is.na(names(query)) | !nzchar(names(query))
  names(query)[empty_names] <- paste0("query_", which(empty_names))
  stats::setNames(gsub("[^ACGTRYSWKMBDHVN]", "", toupper(unname(query))), names(query))
}

#' Identify user sequences against the barcoding library (step 10)
#'
#' Classifies each query against the IdTaxa model of each locus of the library, or of the one
#' declared, and answers locus by locus with the three states of the branch. It never combines
#' loci into one answer (amendment E8 of the validation plan).
#'
#' @details
#' For each locus:
#' 1. The whole query is tested against the 20-mers of the locus on both strands, with the rule of
#'    the rest of the branch (share of 0.15, ratio of 3). If it matches, it is classified as it is
#'    (`path = "whole"`). This is the path of a Sanger sequence.
#' 2. If it does not match and it is longer than 1.5 times the longest library sequence of the
#'    locus, the region of the query with the most 20-mer hits in a window of that width is cropped
#'    and tested again (`path = "crop"`). `other_windows` counts the other windows with at least
#'    half as many hits, since an inverted repeat or an rDNA array can hold more than one copy.
#' 3. The region is cut to the core of the locus, the span of alignment columns covered by half or
#'    more of its library sequences: beyond its first and last 20-mer of the core, the query keeps
#'    only as many bases as the core itself holds there (`core_trimmed` gives the bases removed). A
#'    segment held by one or two references only is left out, so it cannot name them.
#' 4. When IdTaxa assigns a genus (state 1 or 2) with 3 or more sequences in the locus, the region as it was before
#'    step 3 is cut to the core of that genus, with the same rule, and classified again if that is
#'    not the region already classified (`genus_core_trimmed` gives the bases removed from it). The
#'    core of a genus can be narrower than that of the locus (two long references of a genus of
#'    short sequences cannot then take the queries of its absent species) or wider (a genus whose
#'    sequences are all longer than the core of the locus keeps them).
#' 5. State 3 is returned with its reason: `no_overlap` (no match), `short_overlap` (a region
#'    shorter than `min_overlap`) or `low_confidence` (IdTaxa's genus confidence under `threshold`).
#'
#' The model of each locus is IdTaxa trained on the whole library of the locus with a fixed seed,
#' the same model as CN2 in step 8. It is cached in `output_dir/models/` with the md5 of its
#' library file and trained again only when the library or the seed change. Each query has its own
#' seed, from `seed`, the locus and its identifier (the accession of a `Genus_species|accession`
#' name, the whole name otherwise). The random state of the session is left as it was.
#'
#' Each row carries, for its locus, the wrong-species rate of IdTaxa measured in step 11 with the
#' species present (`validation_ws_rate_species_present`, scheme E) and absent
#' (`validation_ws_rate_species_absent`, scheme G), when that table exists and was measured at the
#' same threshold.
#'
#' Raw reads (FASTQ) are not accepted: they need assembly or mapping first.
#'
#' @param query A FASTA or GenBank flat-file path (a GenBank record is named by its accession), a
#'   `Biostrings::DNAStringSet` or a named character vector of sequences.
#' @param locus Character or `NULL`. Loci to test; `NULL` tests every locus of the library.
#' @param library_dir Character. Directory of step 4, with the `LIB_<locus>.fasta` files.
#' @param metrics_dir Character. Directory of step 11, with
#'   `TABLE_barcoding_metrics_summary_idtaxa.csv`.
#' @param output_dir Character. Where the table and the model cache are written.
#' @param run_name Character. Suffix of the output table, `TABLE_barcoding_identify_<run_name>.csv`.
#' @param threshold Numeric. IdTaxa confidence threshold. Defaults to 60, the operating threshold of
#'   Phase 6B.
#' @param min_overlap Integer. Minimum length of the region to classify. Defaults to 100.
#' @param seed Integer. Seed of the models and of the queries. Defaults to 1, as in Tutorial 5.
#' @param processors Integer. Processors for `DECIPHER::IdTaxa()`.
#' @return Invisibly, the table of answers: one row per query and locus tested.
#' @seealso [run_barcoding_controls()], [sweep_barcoding_threshold()],
#'   [summarise_barcoding_metrics()].
#' @examples
#' \dontrun{
#' identify_barcoding_query("my_sequences.fasta")
#' identify_barcoding_query("my_plastome.fasta", run_name = "plastome")
#' }
#' @export
identify_barcoding_query <- function(query,
                                     locus = NULL,
                                     library_dir = file.path("11_barcoding", "4_library"),
                                     metrics_dir = file.path("11_barcoding", "11_metrics"),
                                     output_dir = file.path("11_barcoding", "10_identify"),
                                     run_name = "query",
                                     threshold = 60,
                                     min_overlap = 100L,
                                     seed = 1L,
                                     processors = 1L) {
  q <- .bc_identify_read_query(query)
  .bc_assert_output_dir(output_dir)
  lib <- .bc_library_from_dir(library_dir)
  loci_all <- sort(unique(lib$locus), method = "radix")
  if (is.null(locus)) {
    loci <- loci_all
  } else {
    missing_loci <- setdiff(locus, loci_all)
    if (length(missing_loci) > 0L) {
      stop("Locus not in the library: ", paste(missing_loci, collapse = ", "), ".", call. = FALSE)
    }
    loci <- loci_all[loci_all %in% locus]
  }
  restore_rng <- .bc_rng_state()
  on.exit(restore_rng(), add = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  models_dir <- file.path(output_dir, "models")
  ids <- .bc_parse_header(names(q))$sid

  rows <- list()
  for (l in loci) {
    wl <- .bc_identify_model(library_dir, lib, l, seed, models_dir)
    width <- floor(1.5 * max(nchar(wl$lib_seqs)))
    for (i in seq_along(q)) {
      r <- .bc_identify_one(unname(q[i]), ids[i], l, wl, width, threshold, min_overlap, seed)
      rows[[length(rows) + 1L]] <- cbind(
        data.frame(query = names(q)[i], locus = l, query_length = nchar(q[i]), .qi = i,
                   stringsAsFactors = FALSE), r)
    }
  }
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$.qi, tab$locus, method = "radix"), , drop = FALSE]
  tab$threshold <- threshold
  ctx <- .bc_identify_context(metrics_dir, loci, threshold)
  tab <- cbind(tab, ctx[match(tab$locus, ctx$locus), c("validation_ws_rate_species_present",
                                                        "validation_ws_rate_species_absent")])
  cols <- c("query", "locus", "query_length", "path", "region_start", "region_end",
            "region_length", "other_windows", "orientation", "state", "predicted_species",
            "predicted_genus", "candidates", "genus_idtaxa", "species_idtaxa", "genus_confidence",
            "species_confidence", "reason", "core_trimmed", "genus_core_trimmed", "threshold", "validation_ws_rate_species_present",
            "validation_ws_rate_species_absent")
  tab <- tab[, cols]
  tab$query_length <- as.integer(tab$query_length)
  tab$region_start <- as.integer(tab$region_start)
  tab$region_end <- as.integer(tab$region_end)
  tab$region_length <- as.integer(tab$region_length)
  tab$other_windows <- as.integer(tab$other_windows)
  tab$core_trimmed <- as.integer(tab$core_trimmed)
  tab$genus_core_trimmed <- as.integer(tab$genus_core_trimmed)
  rownames(tab) <- NULL

  for (nm in unique(tab$query)) {
    if (all(tab$reason[tab$query == nm] %in% "no_overlap")) {
      message("Not assignable: no overlap with the library (", nm, ").")
    }
  }
  out_file <- file.path(output_dir, paste0("TABLE_barcoding_identify_", run_name, ".csv"))
  utils::write.csv(tab, out_file, row.names = FALSE)
  .bc_identify_banner(tab, out_file)
  invisible(tab)
}

#' @noRd
.bc_identify_banner <- function(tab, out_file) {
  cat("\n  ====================================================\n")
  cat("    Barcoding Identification Complete \U0001F335\n")
  cat("  ====================================================\n")
  cat("    Queries:          ", length(unique(tab$query)), "\n")
  cat("    Loci tested:      ", length(unique(tab$locus)), "\n")
  for (s in 1:3) cat("    State ", s, " rows:     ", sum(tab$state == s), "\n", sep = "")
  cat("    Table:            ", out_file, "\n")
  cat("  ====================================================\n\n")
  invisible(NULL)
}
