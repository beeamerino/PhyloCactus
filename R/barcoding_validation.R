# ------------------------------------------------------------------------------
# Molecular diagnostic branch (11_barcoding/), validation partitions. PhyloCactus 0.5.0, Phase 3.
#
# Hard rule of the validation plan (sec. 2): no accession may be in the training set and in the
# evaluation set of the same fold. Scheme E leaves one sequence out and asks for the species of the
# accession left out; scheme G leaves one whole species out and asks for its genus. The folds are
# built on the library of 4_library/ (representatives, with the identical sequences already
# collapsed), are deterministic, and are checked before they are written.
# ------------------------------------------------------------------------------

#' Library of the branch, read from the FASTA files of step 4
#' @param library_dir Character. Output directory of [finalize_barcoding_library()].
#' @return Data frame with `sid`, `species`, `genus` and `locus`, one row per sequence.
#' @noRd
.bc_library_from_dir <- function(library_dir) {
  files <- list.files(library_dir, pattern = "^LIB_.*\\.fasta$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No LIB_<locus>.fasta in ", library_dir, ". Run finalize_barcoding_library() first.", call. = FALSE)
  }
  out <- do.call(rbind, lapply(files, function(f) {
    locus <- sub("^LIB_(.*)\\.fasta$", "\\1", basename(f))
    h <- names(Biostrings::readDNAStringSet(f))
    p <- .bc_parse_header(h)
    data.frame(sid = p$sid, species = p$species, genus = .bc_genus(p$species), locus = locus,
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

#' Validation folds of the branch, scheme E (species) and scheme G (genus)
#'
#' Scheme E, one sequence left out: for every species with two or more sequences in the locus, one
#' fold per sequence. Scheme G, one species left out: for every species of a genus with two or more
#' species in the locus, one fold with all its sequences. Species without replica and genera with one
#' species stay in the training set and are never evaluated (validation plan, sec. 3). Deterministic
#' and without a seed (E9): the order is locus, species and `sid`, in the C locale.
#'
#' @param library Data frame with `sid`, `species`, `locus` and, optionally, `genus`.
#' @param scheme `"species"` or `"genus"`.
#' @param loci Character vector or `NULL` for every locus of `library`.
#' @return List of folds. Each fold is a list with `scheme`, `locus`, `stratum`, `target`,
#'   `train_ids` and `test_ids`, the last two as vectors of `sid`.
#' @noRd
.barcoding_folds <- function(library, scheme = c("species", "genus"), loci = NULL) {
  scheme <- match.arg(scheme)
  miss <- setdiff(c("sid", "species", "locus"), names(library))
  if (length(miss) > 0) stop("`library` lacks columns: ", paste(miss, collapse = ", "), call. = FALSE)
  if (is.null(library$genus)) library$genus <- .bc_genus(library$species)
  if (is.null(loci)) loci <- sort(unique(as.character(library$locus)), method = "radix")

  folds <- list()
  for (l in loci) {
    d <- library[library$locus == l, , drop = FALSE]
    d <- d[order(d$species, d$sid, method = "radix"), , drop = FALSE]
    if (nrow(d) == 0) next
    if (scheme == "species") {
      n_by_sp <- table(d$species)
      strata <- sort(names(n_by_sp)[n_by_sp >= 2L], method = "radix")
      for (sp in strata) {
        for (s in sort(d$sid[d$species == sp], method = "radix")) {
          folds[[length(folds) + 1L]] <- list(
            scheme = scheme, locus = l, stratum = sp, target = "species",
            train_ids = setdiff(d$sid, s), test_ids = s
          )
        }
      }
    } else {
      sp_by_gen <- tapply(d$species, d$genus, function(x) length(unique(x)))
      gen <- names(sp_by_gen)[sp_by_gen >= 2L]
      strata <- sort(unique(d$species[d$genus %in% gen]), method = "radix")
      for (sp in strata) {
        test_ids <- sort(d$sid[d$species == sp], method = "radix")
        folds[[length(folds) + 1L]] <- list(
          scheme = scheme, locus = l, stratum = sp, target = "genus",
          train_ids = setdiff(d$sid, test_ids), test_ids = test_ids
        )
      }
    }
  }
  folds
}

#' Check the hard rule of the validation plan on a set of folds
#'
#' Reports, and does not repair: every fold whose training set shares an accession with its
#' evaluation set, and, in scheme G, every fold whose training set holds the species left out.
#'
#' @return Data frame of violations, with zero rows when there is none.
#' @noRd
.barcoding_check_folds <- function(folds, library) {
  rows <- lapply(seq_along(folds), function(i) {
    f <- folds[[i]]
    out <- list()
    shared <- intersect(f$train_ids, f$test_ids)
    if (length(shared) > 0) {
      out[[length(out) + 1L]] <- data.frame(
        fold = i, locus = f$locus, scheme = f$scheme, stratum = f$stratum,
        reason = "accession_in_training_and_test",
        detail = paste(sort(shared, method = "radix"), collapse = ";"), stringsAsFactors = FALSE)
    }
    if (identical(f$scheme, "genus")) {
      d <- library[library$locus == f$locus, , drop = FALSE]
      sp_train <- unique(d$species[d$sid %in% f$train_ids])
      if (f$stratum %in% sp_train) {
        out[[length(out) + 1L]] <- data.frame(
          fold = i, locus = f$locus, scheme = f$scheme, stratum = f$stratum,
          reason = "test_species_in_training",
          detail = paste(sort(d$sid[d$species == f$stratum & d$sid %in% f$train_ids], method = "radix"),
                          collapse = ";"), stringsAsFactors = FALSE)
      }
    }
    if (length(out) == 0) NULL else do.call(rbind, out)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    return(data.frame(fold = integer(0), locus = character(0), scheme = character(0),
                      stratum = character(0), reason = character(0), detail = character(0),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Permute the training labels of one fold (machinery of CN1)
#'
#' Permutes the labels of the training set among themselves, leaves the labels of the evaluation set
#' untouched and restores the session's random state, so that calling it changes nothing outside.
#'
#' @param fold One fold of `.barcoding_folds()`.
#' @param labels Character vector of labels named by `sid`.
#' @param seed Integer. The same seed gives the same permutation.
#' @noRd
.barcoding_permute_labels <- function(fold, labels, seed) {
  ids <- c(fold$train_ids, fold$test_ids)
  out <- labels[ids]
  n <- length(fold$train_ids)
  if (n > 1L) {
    had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had) {
      old <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    } else {
      on.exit(suppressWarnings(rm(list = ".Random.seed", envir = globalenv())), add = TRUE)
    }
    set.seed(seed)
    out[seq_len(n)] <- out[sample.int(n)]
  }
  out
}

#' Nearest neighbour by p distance, reference implementation for the controls
#'
#' Compares only the positions where both sequences carry A, C, G or T, so gaps and ambiguities do
#' not count. Ties are broken by the alphabetically smallest `sid`, so the result is deterministic.
#' This is the reference classifier of the negative controls and of the tests, **not** the classifier
#' of the branch, which is chosen in Phase 6.
#'
#' @param train Character vector of aligned sequences named by `sid`.
#' @param train_labels Character vector of labels named by `sid`.
#' @param query Character. One aligned sequence.
#' @noRd
.bc_nn_predict <- function(train, train_labels, query) {
  q <- strsplit(toupper(query), "", fixed = TRUE)[[1]]
  bases <- c("A", "C", "G", "T")
  d <- vapply(train, function(s) {
    x <- strsplit(toupper(s), "", fixed = TRUE)[[1]]
    n <- min(length(x), length(q))
    if (n == 0L) return(NA_real_)
    xi <- x[seq_len(n)]
    qi <- q[seq_len(n)]
    ok <- xi %in% bases & qi %in% bases
    if (!any(ok)) return(NA_real_)
    mean(xi[ok] != qi[ok])
  }, numeric(1))
  names(d) <- names(train)
  if (all(is.na(d))) return(NA_character_)
  best <- names(d)[!is.na(d) & d == min(d, na.rm = TRUE)]
  unname(train_labels[sort(best, method = "radix")[1]])
}

#' Build the Validation Folds of the Molecular Diagnostic Branch
#'
#' Step 5 of the branch. Reads the library written by [finalize_barcoding_library()], builds the two
#' partition schemes of the validation plan, checks the hard rule on every fold and writes the fold
#' tables and their summary. The training set is not written: it is the rest of the library of the
#' locus, and writing it would multiply the size of the tables without adding information.
#'
#' @param library_dir Character. Output directory of [finalize_barcoding_library()].
#' @param output_dir Character. Destination; refused inside a directory of the phylogeny.
#' @param schemes Character vector. Schemes to build; defaults to both.
#' @return Invisibly, a list with `folds` (by scheme) and `summary`. Writes
#'   `TABLE_barcoding_folds_species.csv`, `TABLE_barcoding_folds_genus.csv` and
#'   `TABLE_barcoding_folds_summary.csv`.
#' @examples
#' \dontrun{
#' build_barcoding_folds(
#'   library_dir = "11_barcoding/4_library",
#'   output_dir = "11_barcoding/5_folds"
#' )
#' }
#' @export
build_barcoding_folds <- function(library_dir = file.path("11_barcoding", "4_library"),
                                  output_dir = file.path("11_barcoding", "5_folds"),
                                  schemes = c("species", "genus")) {
  .bc_assert_output_dir(output_dir)
  schemes <- match.arg(schemes, choices = c("species", "genus"), several.ok = TRUE)
  lib <- .bc_library_from_dir(library_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  folds_by_scheme <- list()
  summaries <- list()
  for (sc in schemes) {
    folds <- .barcoding_folds(lib, scheme = sc)
    if (length(folds) == 0) {
      stop("Scheme '", sc, "': no evaluable fold in ", library_dir, ".", call. = FALSE)
    }
    violations <- .barcoding_check_folds(folds, lib)
    if (nrow(violations) > 0) {
      stop("Scheme '", sc, "': ", nrow(violations), " fold(s) break the rule that no accession is in ",
           "training and evaluation at once. First: ", violations$locus[1], ", ", violations$stratum[1],
           ", ", violations$reason[1], ".", call. = FALSE)
    }
    folds_by_scheme[[sc]] <- folds

    loci <- vapply(folds, function(f) f$locus, character(1))
    tab <- do.call(rbind, lapply(sort(unique(loci), method = "radix"), function(l) {
      ff <- folds[loci == l]
      do.call(rbind, lapply(seq_along(ff), function(i) {
        data.frame(locus = l, fold = i, stratum = ff[[i]]$stratum, sid = ff[[i]]$test_ids,
                   stringsAsFactors = FALSE)
      }))
    }))
    utils::write.csv(tab, file.path(output_dir, paste0("TABLE_barcoding_folds_", sc, ".csv")), row.names = FALSE)

    summaries[[sc]] <- do.call(rbind, lapply(sort(unique(loci), method = "radix"), function(l) {
      ff <- folds[loci == l]
      d <- lib[lib$locus == l, , drop = FALSE]
      strata <- unique(vapply(ff, function(f) f$stratum, character(1)))
      gen_eval <- unique(d$genus[d$species %in% strata])
      ctx <- .bc_locus_summary(d)
      data.frame(
        locus = l, scheme = sc, folds = length(ff),
        accessions_test = sum(vapply(ff, function(f) length(f$test_ids), integer(1))),
        species_testable = length(strata),
        species_training_only = length(unique(d$species)) - length(strata),
        genera_testable = length(gen_eval),
        genera_training_only = length(unique(d$genus)) - length(gen_eval),
        ctx[, setdiff(names(ctx), "locus"), drop = FALSE],
        stringsAsFactors = FALSE
      )
    }))
    message("Scheme '", sc, "': ", length(folds), " folds in ", length(unique(loci)), " loci; ",
            "no fold breaks the hard rule.")
  }

  summary_tab <- do.call(rbind, summaries)
  rownames(summary_tab) <- NULL
  utils::write.csv(summary_tab, file.path(output_dir, "TABLE_barcoding_folds_summary.csv"), row.names = FALSE)

  invisible(list(folds = folds_by_scheme, summary = summary_tab))
}
