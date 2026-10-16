# Cut of long records in step 1 (Phase 10, decision L2 of BMM, 28-09).

#' Records over `max_length` cut to the extent of the shorter records of their locus
#'
#' The mining keeps records of up to 5000 bases (L1). In every locus, a record longer than
#' `max_length` is cut to the span of its 20-mers found in two or more records of `max_length` bases
#' or fewer of the same locus, on either strand; the hits are split where two lie more than
#' `max_gap` bases apart and the group with the most hits is kept, so that a stray 20-mer in a flank
#' does not stretch the cut. The cut is taken on the record's own strand. A long record with no such
#' hit, or with a cut shorter than `min_region`, is left out. A locus keeps one extent: the trnK
#' intron around `matK` does not enter the `matK` library.
#' @return A list with `registry`, `sequences` and `table` (one row per long record: sid, locus,
#'   length, region_start, region_end, region_length, action `cut` or `no_overlap`).
#' @noRd
.bc_cut_long_records <- function(registry, sequences, max_length = 2000L, k = 20L, min_support = 2L,
                                 max_gap = 100L, min_region = 100L) {
  len <- nchar(sequences[registry$sid])
  long <- registry$sid[len > max_length]
  rows <- list()
  drop <- character(0)
  for (l in unique(registry$locus[registry$sid %in% long])) {
    short <- toupper(unname(sequences[registry$sid[registry$locus == l & len <= max_length]]))
    kmers <- unlist(lapply(short, function(s) {
      if (nchar(s) < k) return(character(0))
      rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))
      unique(c(substring(s, 1:(nchar(s) - k + 1L), k:nchar(s)), substring(rc, 1:(nchar(rc) - k + 1L), k:nchar(rc))))
    }), use.names = FALSE)
    tab <- table(kmers)
    pool <- names(tab)[tab >= min_support]
    for (sid in registry$sid[registry$locus == l & registry$sid %in% long]) {
      s <- toupper(sequences[[sid]])
      n <- nchar(s)
      pos <- which(substring(s, 1:(n - k + 1L), k:n) %in% pool)
      reg <- NULL
      if (length(pos)) {
        grp <- cumsum(c(1L, diff(pos) > max_gap))
        best <- as.integer(names(which.max(table(grp))))
        p <- pos[grp == best]
        reg <- c(min(p), max(p) + k - 1L)
        if (reg[2] - reg[1] + 1L < min_region) reg <- NULL
      }
      if (is.null(reg)) {
        drop <- c(drop, sid)
        rows[[length(rows) + 1L]] <- data.frame(sid = sid, locus = l, length = n, region_start = NA_integer_,
                                                region_end = NA_integer_, region_length = NA_integer_,
                                                action = "no_overlap", stringsAsFactors = FALSE)
      } else {
        sequences[[sid]] <- substr(sequences[[sid]], reg[1], reg[2])
        rows[[length(rows) + 1L]] <- data.frame(sid = sid, locus = l, length = n, region_start = reg[1],
                                                region_end = reg[2], region_length = reg[2] - reg[1] + 1L,
                                                action = "cut", stringsAsFactors = FALSE)
      }
    }
  }
  table_out <- if (length(rows)) do.call(rbind, rows) else
    data.frame(sid = character(0), locus = character(0), length = integer(0), region_start = integer(0),
               region_end = integer(0), region_length = integer(0), action = character(0))
  list(registry = registry[!registry$sid %in% drop, , drop = FALSE],
       sequences = sequences[!names(sequences) %in% drop], table = table_out)
}
