# ------------------------------------------------------------------------------
# Molecular diagnostic branch (11_barcoding/), identification of a genome assembly. PhyloCactus 0.5.0.
#
# The crop of step 10 reads every 20-mer of a query and does not scale to an assembly of 1 Gb: the
# loci are located first with minimap2, cut with samtools faidx and handed to step 10.
#
# Governing documents: 01_minutes_and_phases/04_barcoding_0.5.0/ of the audit repository, T3 proposal
# of 2026-09-28 (decisions A1 to A5) and decision T3a of BMM (28-09).
# ------------------------------------------------------------------------------

#' Copy of the assembly as plain FASTA, when it is gzipped
#' @noRd
.bc_assembly_plain <- function(path, tmp) {
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) return(path)
  out <- file.path(tmp, "assembly.fna")
  con_in <- gzfile(path, "rb")
  con_out <- file(out, "wb")
  on.exit({ close(con_in); close(con_out) }, add = TRUE)
  repeat {
    b <- readBin(con_in, "raw", 1e7)
    if (length(b) == 0L) break
    writeBin(b, con_out)
  }
  out
}

#' Best minimap2 hit per locus, and the number of scaffolds hit
#' @noRd
.bc_assembly_hits <- function(paf) {
  if (!file.exists(paf) || file.size(paf) == 0) return(NULL)
  h <- utils::read.delim(paf, header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")[, 1:12]
  names(h) <- c("qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "nmatch",
                "alnlen", "mapq")
  h$locus <- sub("__.*$", "", h$qname)
  n <- vapply(split(h$tname, h$locus), function(x) length(unique(x)), integer(1))
  h <- h[order(h$locus, -h$nmatch, h$tname, h$tstart, method = "radix"), , drop = FALSE]
  top <- h[!duplicated(h$locus), , drop = FALSE]
  top$scaffolds_hit <- unname(n[top$locus])
  rownames(top) <- NULL
  top
}

#' Identify a genome assembly against the library
#'
#' Step 10 for a genome assembly (input type T3). For each locus of the library, its sequences
#' (without gaps) are mapped to the assembly with `minimap2 -x asm20`; the hit with the most matching
#' bases is kept, cut with `samtools faidx` with `flank` bases on each side (within the scaffold) and
#' identified with the rules of [identify_barcoding_query()] for that locus: crop, cut to the core of
#' the locus (J1) and of the genus (J3b), IdTaxa, three states. A locus with no hit is state 3,
#' `no_overlap`.
#'
#' The best hit is kept even when the locus has more than one copy in the assembly (decision T3a):
#' the table gives the scaffold, its length, the mapping quality and the number of scaffolds hit, and
#' a mapping quality under `low_mapq` is flagged and reported, since the region may be a second copy
#' (for instance a plastid insertion in the nucleus) rather than the locus. In the real-file test of
#' 28-09 no locus of two assemblies was named to species; the answers are mostly state 3.
#'
#' Requires `minimap2` and `samtools` on the PATH (for instance from bioconda). The assembly may be
#' gzipped; its folder is not written to (the index goes to a temporary folder).
#'
#' @param assembly Character. Path to the assembly, FASTA, plain or gzipped.
#' @param locus Character or `NULL`. Loci to look for; `NULL`, every locus of the library.
#' @param library_dir,metrics_dir,output_dir,threshold,min_overlap,seed As in
#'   [identify_barcoding_query()].
#' @param run_name Character. The table is `TABLE_barcoding_identify_<run_name>.csv`.
#' @param flank Integer. Bases added on each side of the hit. Defaults to 500.
#' @param low_mapq Integer. Mapping quality under which a hit is flagged. Defaults to 20.
#' @param threads Integer. Threads of `minimap2`.
#' @param minimap2,samtools Character. The executables.
#' @param report Logical. Write the identification report ([report_barcoding_identification()]).
#'   The run record `RUN_<run_name>.csv` is written in any case.
#' @return Invisibly, the table of [identify_barcoding_query()], one row per locus, with the columns
#'   `scaffold`, `scaffold_length`, `hit_start`, `hit_end`, `mapq`, `matching_bases`,
#'   `scaffolds_hit` and `low_mapq`.
#' @seealso [identify_barcoding_query()].
#' @examples
#' \dontrun{
#' identify_barcoding_assembly("GCA_024363205.1_genomic.fna.gz")
#' }
#' @export
identify_barcoding_assembly <- function(assembly,
                                        locus = NULL,
                                        library_dir = file.path("11_barcoding", "4_library"),
                                        metrics_dir = file.path("11_barcoding", "11_metrics"),
                                        output_dir = file.path("11_barcoding", "10_identify"),
                                        run_name = "assembly",
                                        flank = 500L,
                                        low_mapq = 20L,
                                        threads = 1L,
                                        minimap2 = "minimap2",
                                        samtools = "samtools",
                                        threshold = 60,
                                        min_overlap = 100L,
                                        seed = 1L,
                                        report = TRUE) {
  if (!file.exists(assembly)) stop("Assembly not found: ", assembly, ".", call. = FALSE)
  exe <- Sys.which(c(minimap2, samtools))
  if (any(!nzchar(exe))) {
    stop("Not found on the PATH: ", paste(c(minimap2, samtools)[!nzchar(exe)], collapse = ", "),
         ". The identification of an assembly needs minimap2 and samtools (bioconda).", call. = FALSE)
  }
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
  tmp <- tempfile("bc_assembly_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  # 1. Library sequences of the loci, without gaps, named <locus>__<sid>
  qs <- unlist(lapply(loci, function(l) {
    x <- Biostrings::readDNAStringSet(file.path(library_dir, paste0("LIB_", l, ".fasta")))
    s <- gsub("-", "", as.character(x))
    stats::setNames(s, paste0(l, "__", .bc_parse_header(names(x))$sid))
  }))
  qs <- qs[nchar(qs) > 0L]
  qf <- file.path(tmp, "library.fasta")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(qs), qf)

  # 2. minimap2, best hit per locus (A1)
  fna <- .bc_assembly_plain(assembly, tmp)
  paf <- file.path(tmp, "hits.paf")
  st <- system2(exe[[1]], c("-x", "asm20", "-c", "--secondary=no", "-t", as.integer(threads), "-o", shQuote(paf),
                            shQuote(fna), shQuote(qf)), stdout = FALSE, stderr = FALSE)
  if (!identical(as.integer(st), 0L)) stop("minimap2 failed (exit status ", st, ").", call. = FALSE)
  top <- .bc_assembly_hits(paf)
  if (!is.null(top)) top <- top[top$locus %in% loci, , drop = FALSE]

  # 3. Regions with their flanks, cut with samtools faidx; the index stays in the temporary folder
  regions <- character(0)
  if (!is.null(top) && nrow(top) > 0L) {
    top$from <- pmax(1L, top$tstart + 1L - as.integer(flank))
    top$to <- pmin(top$tlen, top$tend + as.integer(flank))
    fai <- file.path(tmp, "assembly.fai")
    st <- system2(exe[[2]], c("faidx", "--fai-idx", shQuote(fai), shQuote(fna)), stdout = FALSE, stderr = FALSE)
    if (!identical(as.integer(st), 0L)) stop("samtools faidx failed (exit status ", st, ").", call. = FALSE)
    for (k in seq_len(nrow(top))) {
      reg <- sprintf("%s:%d-%d", top$tname[k], top$from[k], top$to[k])
      out <- system2(exe[[2]], c("faidx", "--fai-idx", shQuote(fai), shQuote(fna), shQuote(reg)),
                     stdout = TRUE, stderr = FALSE)
      regions[top$locus[k]] <- toupper(paste(out[!startsWith(out, ">")], collapse = ""))
    }
  }

  # 4. Identification of each region for its locus
  rows <- list()
  for (l in loci) {
    if (l %in% names(regions) && nzchar(regions[[l]])) {
      q <- stats::setNames(regions[[l]], paste0(run_name, "__", l))
      r <- .bc_identify_table(q, l, lib, library_dir, metrics_dir, output_dir, threshold, min_overlap, seed)
      h <- top[top$locus == l, ]
      r$scaffold <- h$tname
      r$scaffold_length <- as.integer(h$tlen)
      r$hit_start <- as.integer(h$tstart + 1L)
      r$hit_end <- as.integer(h$tend)
      r$mapq <- as.integer(h$mapq)
      r$matching_bases <- as.integer(h$nmatch)
      r$scaffolds_hit <- as.integer(h$scaffolds_hit)
      r$low_mapq <- h$mapq < low_mapq
      if (r$low_mapq) {
        message("Locus ", l, ": mapping quality ", h$mapq, " (under ", low_mapq, "), ", h$scaffolds_hit,
                " scaffolds hit; the locus has more than one copy in the assembly and the best hit is kept.")
      }
    } else {
      # No hit: the empty query of step 10 gives the row of a locus without overlap
      r <- .bc_identify_table(stats::setNames("", paste0(run_name, "__", l)), l, lib, library_dir, metrics_dir,
                              output_dir, threshold, min_overlap, seed)
      r$scaffold <- NA_character_
      r$scaffold_length <- NA_integer_
      r$hit_start <- NA_integer_
      r$hit_end <- NA_integer_
      r$mapq <- NA_integer_
      r$matching_bases <- NA_integer_
      r$scaffolds_hit <- 0L
      r$low_mapq <- FALSE
    }
    rows[[l]] <- r
  }
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  out_file <- file.path(output_dir, paste0("TABLE_barcoding_identify_", run_name, ".csv"))
  utils::write.csv(tab, out_file, row.names = FALSE)
  .bc_write_run_record(output_dir, run_name, input_path = assembly, input_md5 = unname(tools::md5sum(assembly)),
                       input_type = "assembly",
                       route = c("minimap2 -x asm20", "samtools faidx", "crop", "cut to locus core (J1)", "IdTaxa",
                                 "cut to genus core (J3b)"),
                       library_dir = library_dir, loci = loci, models_dir = file.path(output_dir, "models"),
                       threshold = threshold, seed = seed, min_overlap = min_overlap,
                       tools = c(minimap2 = .bc_tool_version(exe[[1]]), samtools = .bc_tool_version(exe[[2]])),
                       extra = data.frame(kind = "setting", name = c("flank", "low_mapq"), value = c(format(flank), format(low_mapq))),
                       loci_declared = locus)
  .bc_identify_banner(tab, out_file)
  if (isTRUE(report)) {
    report_barcoding_identification(run_name, results_dir = output_dir, library_dir = library_dir, metrics_dir = metrics_dir)
  }
  invisible(tab)
}
