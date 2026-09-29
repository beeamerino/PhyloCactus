# Extra records of the library cropped from GenBank plastomes (0_genomic_raw/plastomes/).
# Phase 8 of PhyloCactus 0.5.0, decision F2 and PL1 to PL7 of BMM (28-09): the plastomes fill loci
# missing for species of the library; they enter step 1 beside the phylotaR clusters, with their
# source declared, and nothing is classified here.

#' Default GenBank fetcher: flat files for a query or for a set of accessions
#' @noRd
.bc_fetch_genbank <- function(query = NULL, ids = NULL, batch = 20L) {
  if (!requireNamespace("rentrez", quietly = TRUE)) {
    stop("The default download needs the package rentrez (install.packages(\"rentrez\")), or pass a ",
         "function as `fetch`.", call. = FALSE)
  }
  if (!is.null(ids)) {
    ids <- unique(ids)
    parts <- lapply(split(ids, ceiling(seq_along(ids) / 100)), function(b)
      rentrez::entrez_fetch("nuccore", id = b, rettype = "gb", retmode = "text"))
    return(paste(unlist(parts), collapse = "\n"))
  }
  s <- rentrez::entrez_search("nuccore", query, use_history = TRUE, retmax = 0)
  if (s$count == 0) return("")
  parts <- lapply(seq(0, s$count - 1, by = batch), function(start) {
    Sys.sleep(0.4)
    rentrez::entrez_fetch("nuccore", web_history = s$web_history, rettype = "gb", retmode = "text",
                          retstart = start, retmax = batch)
  })
  paste(unlist(parts), collapse = "\n")
}

#' GenBank flat files split into accession, organism, voucher, sequence and the text of the record
#' @noRd
.bc_parse_genbank <- function(text) {
  lines <- strsplit(paste(text, collapse = "\n"), "\n", fixed = TRUE)[[1]]
  ends <- which(lines == "//")
  starts <- c(1L, utils::head(ends, -1L) + 1L)
  recs <- lapply(seq_along(ends), function(i) lines[starts[i]:ends[i]])
  recs <- recs[vapply(recs, function(r) any(startsWith(r, "LOCUS")), logical(1))]
  do.call(rbind, lapply(recs, function(r) {
    field <- function(q) {
      h <- grep(paste0('/', q, '="'), r, value = TRUE, fixed = TRUE)
      if (length(h)) sub(paste0('.*/', q, '="([^"]*)"?.*'), "\\1", h[1]) else NA_character_
    }
    ver <- grep("^VERSION", r, value = TRUE)
    o <- which(startsWith(r, "ORIGIN"))
    sq <- if (length(o)) gsub("[^A-Za-z]", "", paste(r[(o + 1L):(length(r) - 1L)], collapse = "")) else ""
    data.frame(accession = if (length(ver)) sub("^VERSION +([^ ]+).*$", "\\1", ver[1]) else NA_character_,
               organism = field("organism"), voucher = field("specimen_voucher"),
               sequence = toupper(sq), text = paste(r, collapse = "\n"), stringsAsFactors = FALSE)
  }))
}

#' Same specimen: same collector surname and same number, or the same normalised string
#' @noRd
.bc_same_specimen <- function(a, b) {
  norm <- function(v) {
    v <- tolower(ifelse(is.na(v), "", v))
    key <- ifelse(grepl("[a-z]{3,}", v) & grepl("[0-9]", v),
                  paste0(sub(".*?([a-z]{3,}).*", "\\1", v, perl = TRUE), "#", sub("^[^0-9]*([0-9]+).*$", "\\1", v)),
                  NA_character_)
    list(full = gsub("[^a-z0-9]", "", v), key = key)
  }
  x <- norm(a); y <- norm(b)
  (!is.na(x$key) & !is.na(y$key) & x$key == y$key) | (nzchar(x$full) & x$full == y$full)
}

#' Pools of one locus of the library: strand pool, core pool (J1) and crop window
#' @noRd
.bc_locus_pools <- function(library_dir, locus) {
  dna <- ape::read.dna(file.path(library_dir, paste0("LIB_", locus, ".fasta")), format = "fasta")
  lib_m <- as.matrix(dna)
  seqs <- gsub("-", "", toupper(apply(as.character(lib_m), 1, paste, collapse = "")), fixed = TRUE)
  list(pool = .bc_strand_pool(lib_m), core_pool = .bc_core_pool(lib_m),
       width = floor(1.5 * max(nchar(seqs))), sids = .bc_parse_header(rownames(lib_m))$sid,
       species = .bc_parse_header(rownames(lib_m))$species)
}

#' The region of a locus in a long sequence: the crop of step 10 and the cut of J1
#'
#' Returns NULL when the locus is not found, or a list with the coordinates in `s`, the strand and
#' the region oriented as the library.
#' @noRd
.bc_locate_region <- function(s, lp) {
  n <- nchar(s)
  start <- 1L; end <- n; other <- 0L
  o <- .bc_orient_to_library(.bc_dnabin_row(s, "query"), lp$pool)
  if (o$orientation == "no_match") {
    if (n <= lp$width) return(NULL)
    fwd <- .bc_best_window(.bc_kmer_hits(s, lp$pool), lp$width)
    rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(s)))
    rev <- .bc_best_window(.bc_kmer_hits(rc, lp$pool), lp$width)
    if (is.null(fwd) && is.null(rev)) return(NULL)
    if (is.null(fwd) || (!is.null(rev) && rev$hits > fwd$hits)) {
      start <- n - rev$end + 1L; end <- n - rev$start + 1L; other <- rev$other
    } else {
      start <- fwd$start; end <- fwd$end; other <- fwd$other
    }
    o <- .bc_orient_to_library(.bc_dnabin_row(substr(s, start, end), "query"), lp$pool)
    if (o$orientation == "no_match") return(NULL)
  }
  oriented <- toupper(paste(as.character(as.matrix(o$query)), collapse = ""))
  t <- .bc_core_cut(oriented, lp$core_pool)
  if (is.null(t)) return(NULL)
  if (o$orientation == "reverse") {
    new_start <- end - t[2] + 1L; end <- end - t[1] + 1L
  } else {
    new_start <- start + t[1] - 1L; end <- start + t[2] - 1L
  }
  list(start = as.integer(new_start), end = as.integer(end), strand = o$orientation,
       sequence = substr(oriented, t[1], t[2]), other = as.integer(other))
}

#' Extra library records cropped from GenBank plastomes (step 1, Phase 8)
#'
#' Downloads complete plastomes of Cactaceae from GenBank, crops from each the loci of the library
#' and writes them as extra records that [assemble_barcoding_dataset()] can take with
#' `extra_records_dir`. Nothing is classified: a record takes the name of its plastome.
#'
#' @details
#' - **Names** are matched to the checklist with the rule of step 1 and cut to the binomial, the rank
#'   of the library; unmatched names are left out and listed in `TABLE_unmatched_names.csv`.
#' - **Extraction**, per locus: the crop of [identify_barcoding_query()] (the window with the most
#'   20-mer hits of the locus) and the cut to the core of the locus (columns covered by half or more
#'   of its library sequences), so a record has the extent of the library; regions shorter than
#'   `min_region` are not written. Where the crop finds more than one window (an inverted repeat),
#'   the best one is kept and `other_windows` counts the rest.
#' - **Duplicates:** RefSeq is left out of the query, and a plastome whose sequence (on either strand)
#'   repeats one already downloaded adds nothing; it is listed in `TABLE_duplicate_plastomes.csv`.
#' - **Identifier:** `sid = <accession>__<locus>`, so a sid stays unique across loci.
#' - **Same specimen:** a plastome whose `specimen_voucher` matches that of a library accession of
#'   the same species in the same locus (same collector surname and number, or the same normalised
#'   string) does not add that locus; the pair is listed in `TABLE_skipped_same_specimen.csv`.
#'
#' Writes to `output_dir`: `downloads/<accession>.gb`, `MANIFEST.csv` (accession, organism, length,
#' md5, query, date), `extra_records.fasta` (headers `Genus_species|sid`),
#' `TABLE_genomic_extra_records.csv` (sid, accession, species, organism, locus, region_start,
#' region_end, strand, region_length, other_windows, md5 of the plastome file, library_md5 of the
#' `LIB_<locus>.fasta` that set the extent), `TABLE_skipped_same_specimen.csv` and
#' `TABLE_unmatched_names.csv`.
#'
#' @param query Character. GenBank query of the plastomes.
#' @param library_dir Character. Directory of step 4, with the `LIB_<locus>.fasta` files.
#' @param output_dir Character. Where the downloads and the extra records are written.
#' @param checklist_path Character. The checklist (a path, or a vector of accepted names).
#' @param min_region Integer. Shortest region written. Defaults to 100.
#' @param keep_genera Character. Genera outside the checklist whose plastomes are kept under the
#'   NCBI binomial (for the outgroup: Portulacaceae and Talinaceae, decision L5 of 28-09). A name
#'   without an epithet (`sp.`, `cf.`, `aff.`) is still left out. Defaults to none.
#' @param fetch Function `(query, ids)` returning GenBank flat text; defaults to a fetcher through
#'   `rentrez`.
#' @return Invisibly, a list with `records`, `skipped` and `unmatched`.
#' @seealso [assemble_barcoding_dataset()], [identify_barcoding_query()].
#' @examples
#' \dontrun{
#' extract_barcoding_plastome_loci()
#' }
#' @export
extract_barcoding_plastome_loci <- function(query = paste0("Cactaceae[Organism] AND 100000:170000[SLEN] AND ",
                                                           "(chloroplast OR plastid) AND \"complete genome\" NOT refseq[filter]"),
                                            library_dir = file.path("11_barcoding", "4_library"),
                                            output_dir = file.path("0_genomic_raw", "plastomes"),
                                            checklist_path = system.file("extdata", "CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx",
                                                                         package = "PhyloCactus"),
                                            min_region = 100L,
                                            keep_genera = character(0),
                                            fetch = .bc_fetch_genbank) {
  .bc_assert_output_dir(output_dir)
  dir.create(file.path(output_dir, "downloads"), recursive = TRUE, showWarnings = FALSE)
  pl <- .bc_parse_genbank(fetch(query = query))
  if (is.null(pl) || nrow(pl) == 0L) stop("The query returned no GenBank record.", call. = FALSE)
  pl$md5 <- vapply(seq_len(nrow(pl)), function(i) {
    f <- file.path(output_dir, "downloads", paste0(pl$accession[i], ".gb"))
    writeLines(pl$text[i], f)
    unname(tools::md5sum(f))
  }, "")
  utils::write.csv(data.frame(accession = pl$accession, organism = pl$organism, length = nchar(pl$sequence),
                              md5 = pl$md5, query = query, retrieved = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
                   file.path(output_dir, "MANIFEST.csv"), row.names = FALSE)

  # RefSeq copies (NC_) repeat GenBank records: the same sequence, on either strand, is kept once,
  # the GenBank accession first (28-09: 55 of 161 downloads were such copies)
  pl <- pl[order(grepl("^NC_", pl$accession)), , drop = FALSE]
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(pl$sequence)))
  key <- ifelse(pl$sequence < rc, pl$sequence, rc)
  dup <- duplicated(key)
  duplicates <- data.frame(accession = pl$accession[dup], duplicate_of = pl$accession[match(key[dup], key)],
                           stringsAsFactors = FALSE)
  utils::write.csv(duplicates, file.path(output_dir, "TABLE_duplicate_plastomes.csv"), row.names = FALSE)
  pl <- pl[!dup, , drop = FALSE]

  matched <- .bc_resolve_names(pl$organism, checklist_path)
  # L5 (BMM, 28-09): genera outside the checklist (Portulacaceae, Talinaceae) keep the NCBI binomial,
  # provided it names a species
  w <- strsplit(trimws(pl$organism), "[[:space:]]+")
  epithet <- vapply(w, function(v) if (length(v) >= 2L) v[2] else NA_character_, "")
  keep <- is.na(matched) & vapply(w, `[`, "", 1) %in% keep_genera & !is.na(epithet) & grepl("^[a-z-]+$", epithet) &
    !epithet %in% c("sp", "spp", "cf", "aff", "x")
  matched[keep] <- paste(vapply(w[keep], `[`, "", 1), epithet[keep], sep = "_")
  pl$species <- sub("_(subsp|ssp|var|subvar|f|fo|forma)[.]?_.*$", "", matched)
  unmatched <- pl[is.na(pl$species), c("accession", "organism"), drop = FALSE]
  pl <- pl[!is.na(pl$species), , drop = FALSE]

  loci <- sub("^LIB_(.*)\\.fasta$", "\\1", list.files(library_dir, pattern = "^LIB_.*\\.fasta$"))
  pools <- lapply(stats::setNames(loci, loci), function(l) .bc_locus_pools(library_dir, l))
  # The extent of a record is set by the library it is cut with (N1, 28-09)
  lib_md5 <- vapply(loci, function(l) unname(tools::md5sum(file.path(library_dir, paste0("LIB_", l, ".fasta")))), "")

  # Vouchers of the library accessions of the species that come with a voucher
  lib_v <- data.frame(sid = character(0), voucher = character(0), stringsAsFactors = FALSE)
  sp_v <- unique(pl$species[!is.na(pl$voucher)])
  lib_sids <- unique(unlist(lapply(pools, function(p) p$sids[p$species %in% sp_v])))
  if (length(lib_sids)) {
    lv <- .bc_parse_genbank(fetch(ids = lib_sids))
    if (!is.null(lv)) lib_v <- data.frame(sid = lv$accession, voucher = lv$voucher, stringsAsFactors = FALSE)
  }

  records <- list(); skipped <- list(); seqs <- character(0)
  for (i in seq_len(nrow(pl))) {
    for (l in loci) {
      lp <- pools[[l]]
      reg <- .bc_locate_region(pl$sequence[i], lp)
      if (is.null(reg) || nchar(reg$sequence) < min_region) next
      same <- lp$sids[lp$species == pl$species[i]]
      same <- same[same %in% lib_v$sid]
      hit <- same[.bc_same_specimen(pl$voucher[i], lib_v$voucher[match(same, lib_v$sid)])]
      if (length(hit)) {
        skipped[[length(skipped) + 1L]] <- data.frame(accession = pl$accession[i], species = pl$species[i], locus = l,
                                                      library_sid = hit[1], plastome_voucher = pl$voucher[i],
                                                      library_voucher = lib_v$voucher[match(hit[1], lib_v$sid)],
                                                      stringsAsFactors = FALSE)
        next
      }
      sid <- paste0(pl$accession[i], "__", l)
      records[[length(records) + 1L]] <- data.frame(
        sid = sid, accession = pl$accession[i], species = pl$species[i], organism = pl$organism[i], locus = l,
        region_start = reg$start, region_end = reg$end, strand = reg$strand,
        region_length = nchar(reg$sequence), other_windows = reg$other, md5 = pl$md5[i],
        library_md5 = lib_md5[[l]], stringsAsFactors = FALSE)
      seqs[paste0(pl$species[i], "|", sid)] <- reg$sequence
    }
  }
  empty_rec <- data.frame(sid = character(0), accession = character(0), species = character(0), organism = character(0),
                          locus = character(0), region_start = integer(0), region_end = integer(0), strand = character(0),
                          region_length = integer(0), other_windows = integer(0), md5 = character(0),
                          library_md5 = character(0))
  records <- if (length(records)) do.call(rbind, records) else empty_rec
  skipped <- if (length(skipped)) do.call(rbind, skipped) else
    data.frame(accession = character(0), species = character(0), locus = character(0), library_sid = character(0),
               plastome_voucher = character(0), library_voucher = character(0))
  utils::write.csv(records, file.path(output_dir, "TABLE_genomic_extra_records.csv"), row.names = FALSE)
  utils::write.csv(skipped, file.path(output_dir, "TABLE_skipped_same_specimen.csv"), row.names = FALSE)
  utils::write.csv(unmatched, file.path(output_dir, "TABLE_unmatched_names.csv"), row.names = FALSE)
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(seqs), file.path(output_dir, "extra_records.fasta"))
  message("Plastomes: ", nrow(pl) + nrow(unmatched) + nrow(duplicates), " downloaded, ", nrow(duplicates),
          " duplicates, ", nrow(unmatched), " with a name out of the ",
          "checklist; extra records: ", nrow(records), " in ", length(unique(records$locus)), " loci; skipped as ",
          "the same specimen: ", nrow(skipped), ".")
  invisible(list(records = records, skipped = skipped, unmatched = unmatched, duplicates = duplicates))
}

#' The extra records of `extra_records_dir` appended to the registry and sequences of step 1
#'
#' With `NULL` the registry and the sequences come back unchanged. Otherwise every row gets a
#' `source` column, `phylotaR` or `genbank_plastome`.
#' @noRd
.bc_add_extra_records <- function(registry, sequences, extra_records_dir) {
  if (is.null(extra_records_dir)) return(list(registry = registry, sequences = sequences))
  tab <- utils::read.csv(file.path(extra_records_dir, "TABLE_genomic_extra_records.csv"), stringsAsFactors = FALSE)
  fa <- Biostrings::readDNAStringSet(file.path(extra_records_dir, "extra_records.fasta"))
  ex_seq <- stats::setNames(as.character(fa), .bc_parse_header(names(fa))$sid)
  add <- data.frame(sid = tab$sid, species = tab$species, genus = .bc_genus(tab$species), locus = tab$locus,
                    cluster_id = NA_integer_, genbank_name = gsub(" ", "_", tab$organism),
                    source = "genbank_plastome", stringsAsFactors = FALSE)
  if (is.null(registry$source)) registry$source <- "phylotaR"
  out <- rbind(registry, add)
  if (anyDuplicated(out$sid)) {
    stop("An extra record repeats a sid of the registry: ", paste(utils::head(out$sid[duplicated(out$sid)], 5), collapse = ", "),
         call. = FALSE)
  }
  out <- out[order(out$locus, out$species, out$sid), , drop = FALSE]
  rownames(out) <- NULL
  list(registry = out, sequences = c(sequences, ex_seq[tab$sid]))
}
