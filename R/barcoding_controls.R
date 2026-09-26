# Negative controls of the molecular diagnostic branch (11_barcoding/8_controls/).
# Phase 5B of PhyloCactus 0.5.0, validation plan sec. 4.
#
# The metric functions live here because CN1 needs them. In this phase they are applied only to
# permuted labels, to outgroup queries and to resubstitution: never to the true labels, which is
# Phase 6. The rule is procedural and it is written down rather than assumed.

#' Accuracy of a prediction table, with its base rate beside it
#'
#' What counts as a hit is fixed by sec. 1 of the validation plan and is not reinvented here. For
#' the species target, only state 1 with the right species; a candidate set never counts as a
#' species hit. For the genus target, the right genus, whether it comes from state 1 or from a
#' state 2 whose genus is unique. State 3 is never a hit, and it is not an error either: it is the
#' third state, and it is counted as itself.
#'
#' The base rate is the frequency of the majority class **among the queries evaluated**, not in the
#' library. Publishing an accuracy without it says nothing with this class imbalance.
#' @noRd
.bc_metric_accuracy <- function(pred) {
  if (nrow(pred) == 0L) {
    return(data.frame(locus = character(0), scheme = character(0), method = character(0),
                      queries = integer(0), hits = integer(0), accuracy = numeric(0),
                      base_rate = numeric(0), stringsAsFactors = FALSE))
  }
  key <- paste(pred$locus, pred$scheme, pred$method, sep = "\r")
  out <- do.call(rbind, lapply(sort(unique(key), method = "radix"), function(k) {
    d <- pred[key == k, , drop = FALSE]
    species_scheme <- d$scheme[1] == "species"
    truth <- if (species_scheme) d$true_species else d$true_genus
    is_hit <- if (species_scheme) {
      d$state == 1L & !is.na(d$predicted_species) & d$predicted_species == d$true_species
    } else {
      d$state %in% c(1L, 2L) & !is.na(d$predicted_genus) & d$predicted_genus == d$true_genus
    }
    data.frame(locus = d$locus[1], scheme = d$scheme[1], method = d$method[1],
               queries = nrow(d), hits = as.integer(sum(is_hit)),
               accuracy = sum(is_hit) / nrow(d),
               base_rate = max(table(truth)) / length(truth),
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

#' Verdict of CN1, with the rule of the plan written out
#'
#' The accuracy with permuted labels must fall to the base rate. The threshold is the base rate plus
#' three standard deviations of the permutation distribution, and what is published is the whole
#' distribution, never a single value. Above the threshold there is leakage, the real result is not
#' interpreted and the branch goes back to Phase 3.
#' @noRd
.bc_cn1_verdict <- function(accuracies, base_rate) {
  accuracies <- accuracies[!is.na(accuracies)]
  sd_perm <- if (length(accuracies) > 1L) stats::sd(accuracies) else 0
  threshold <- base_rate + 3 * sd_perm
  data.frame(permutations = length(accuracies), mean = mean(accuracies), sd = sd_perm,
             base_rate = base_rate, threshold = threshold, leakage = mean(accuracies) > threshold,
             stringsAsFactors = FALSE)
}

#' Negative Controls of the Molecular Diagnostic Branch
#'
#' Step 8 of the branch, and the gate in front of every real figure (validation plan, sec. 4). Three
#' controls, none of which is a result of the method: the three are the condition for a result of
#' the method to be interpretable.
#'
#' **CN1, permuted labels.** The whole procedure is repeated with the species labels permuted inside
#' the training set. The accuracy must fall to the base rate; if it does not, there is leakage.
#' Written as one row per permutation plus its verdict.
#'
#' **CN2, a query alien to the clade.** Sequences of the outgroup, which are not in the library, are
#' passed as queries. The table publishes how many outgroup **species** each locus had, and how many
#' sequences, zeros included. The size of the control is the species count: GenBank holds population
#' level studies, so a locus can carry dozens of sequences of a single species and control almost
#' nothing. A control over one species is an anecdote and has to be read as one.
#'
#' **CN3, resubstitution.** The same queries as the folds, in both schemes, each one classified
#' against its own training set with the query put back into it. Nothing else changes between this
#' figure and the honest one, which is what makes the two comparable: until 2026-09-25 this control
#' ran over the whole library instead, a different set of queries, and the two figures could not be
#' subtracted. Written with the label `irreproducible_resubstitution` in the table itself, so the
#' figure cannot travel without it. **Its value is published in Phase 6 and only next to the honest
#' figure**, which does not exist yet.
#'
#' @param library_dir Character. Output directory of [finalize_barcoding_library()].
#' @param folds_dir Character. Output directory of [build_barcoding_folds()].
#' @param output_dir Character. Destination; refused inside a directory of the phylogeny.
#' @param outgroup_dir Character or `NULL`. Directory with `<locus>.fasta` of the outgroup, as the
#'   assembly writes them. `NULL` skips CN2 and says so.
#' @param method Character. `"nn"` or `"idtaxa"`. With `"idtaxa"` only CN2 runs (decision K6 of
#'   BMM, 2026-09-26): CN1 would retrain IdTaxa once per fold and permutation, and CN3 measures a
#'   bias of the nearest neighbour. The outgroup queries are oriented against the library with the
#'   rule of the nearest neighbour, classified by IdTaxa trained on the whole locus at threshold 0
#'   with a seed per query, and written to `TABLE_barcoding_cn2_queries_idtaxa.csv`; the tables of
#'   the nearest neighbour are not touched. MAFFT is not needed.
#' @param threshold Numeric. Confidence threshold applied to IdTaxa's output, for `"idtaxa"`, as in
#'   [classify_barcoding_folds()].
#' @param permutations Integer. Permutations of CN1; at least 10 in a real run.
#' @param seed Integer. Seed of the permutations, declared so they are reproducible.
#' @param model,min_comparable Distance model and minimum of comparable positions, as in step 7.
#' @param schemes,loci Schemes and loci; `NULL` for every locus of the library.
#' @param max_folds Integer or `NULL`. Declared subset of folds, for the expensive method.
#' @param controls Character vector. Which of `"CN1"`, `"CN2"` and `"CN3"` to run.
#' @param mafft_exec Character. Command or path of the `MAFFT` binary. CN2 adds every outgroup
#'   query to the alignment of the library with `--add --keeplength`, one call per query, so the
#'   columns of the library stay fixed and the distance of one query is comparable with the
#'   distance of another. Without `outgroup_dir` MAFFT is never called.
#' @param mafft_opts Character. Further options for `MAFFT`. `--add` and `--keeplength` are added
#'   by the function and are not optional.
#' @return Invisibly, a list with the tables. Writes `TABLE_barcoding_cn1_permutations.csv`,
#'   `TABLE_barcoding_cn1_verdict.csv`, `TABLE_barcoding_cn2_outgroup.csv` and
#'   `TABLE_barcoding_cn3_resubstitution.csv`.
#' @examples
#' \dontrun{
#' run_barcoding_controls(
#'   library_dir = "11_barcoding/4_library",
#'   folds_dir = "11_barcoding/5_folds",
#'   output_dir = "11_barcoding/8_controls",
#'   outgroup_dir = "11_barcoding/8_controls/outgroup/1_assembly"
#' )
#' }
#' @export
run_barcoding_controls <- function(library_dir = file.path("11_barcoding", "4_library"),
                                   folds_dir = file.path("11_barcoding", "5_folds"),
                                   output_dir = file.path("11_barcoding", "8_controls"),
                                   outgroup_dir = NULL,
                                   method = c("nn", "idtaxa"),
                                   permutations = 10L,
                                   seed = 1L,
                                   model = "raw",
                                   min_comparable = 100L,
                                   schemes = c("species", "genus"),
                                   loci = NULL,
                                   max_folds = NULL,
                                   controls = c("CN1", "CN2", "CN3"),
                                   mafft_exec = "mafft",
                                   mafft_opts = "--auto",
                                   threshold = 60) {
  method <- match.arg(method)
  if (method == "idtaxa" && any(c("CN1", "CN3") %in% controls)) {
    stop("CN1 and CN3 run with the nearest neighbour only (decision K6); with method = \"idtaxa\" ",
         "ask for controls = \"CN2\".", call. = FALSE)
  }
  .bc_assert_output_dir(output_dir)
  lib <- .bc_library_from_dir(library_dir)
  for (sc in schemes) {
    f <- file.path(folds_dir, paste0("TABLE_barcoding_folds_", sc, ".csv"))
    if (!file.exists(f)) {
      stop("No ", basename(f), " in ", folds_dir, ". Run build_barcoding_folds() first.", call. = FALSE)
    }
  }
  if ("CN2" %in% controls && !is.null(outgroup_dir) && method == "nn") {
    # CN2 aligns every query against the library, so the controls need MAFFT at run time. The
    # check belongs here and not inside the CN2 loop: an absent binary has to cost seconds, not
    # the hour that CN1 takes before reaching it.
    .bc_assert_mafft(mafft_exec)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  loci_todos <- sort(unique(lib$locus), method = "radix")
  loci <- if (is.null(loci)) loci_todos else intersect(loci_todos, loci)

  dist_matrices <- list()
  sequences <- list()
  for (l in loci) {
    dna <- ape::read.dna(file.path(library_dir, paste0("LIB_", l, ".fasta")), format = "fasta")
    rownames(dna) <- .bc_parse_header(labels(dna))$sid
    if (method == "nn") dist_matrices[[l]] <- .bc_classifier_matrix(dna, model, min_comparable)
    sequences[[l]] <- dna
  }

  out <- list()

  if ("CN1" %in% controls) {
    rows <- list()
    for (sc in schemes) {
      tab <- utils::read.csv(file.path(folds_dir, paste0("TABLE_barcoding_folds_", sc, ".csv")),
                             stringsAsFactors = FALSE)
      for (l in loci) {
        d <- lib[lib$locus == l, , drop = FALSE]
        species <- stats::setNames(d$species, d$sid)
        folds <- .bc_sample_folds(.bc_folds_from_table(tab[tab$locus == l, , drop = FALSE], d$sid),
                                     max_folds, seed)
        if (length(folds) == 0L) next
        for (k in seq_len(permutations)) {
          pred <- do.call(rbind, lapply(folds, function(p) {
            label_set <- .barcoding_permute_labels(
              list(train_ids = p$train_ids, test_ids = p$test_ids), species, seed = seed + k)
            do.call(rbind, lapply(p$test_ids, function(q) {
              r <- .bc_classify_nn(dist_matrices[[l]], p$train_ids, q, label_set)
              cbind(data.frame(locus = l, scheme = sc, method = method,
                               true_species = unname(species[q]),
                               true_genus = .bc_genus(unname(species[q])),
                               stringsAsFactors = FALSE), r)
            }))
          }))
          m <- .bc_metric_accuracy(pred)
          rows[[length(rows) + 1L]] <- cbind(data.frame(permutation = k, stringsAsFactors = FALSE), m)
        }
        message(sprintf("CN1, locus '%s', scheme %s: %d permutations over %d folds.",
                        l, sc, permutations, length(folds)))
      }
    }
    perm <- do.call(rbind, rows)
    rownames(perm) <- NULL
    key <- paste(perm$locus, perm$scheme, perm$method, sep = "\r")
    verdict <- do.call(rbind, lapply(sort(unique(key), method = "radix"), function(k) {
      d <- perm[key == k, , drop = FALSE]
      cbind(data.frame(locus = d$locus[1], scheme = d$scheme[1], method = d$method[1],
                       stringsAsFactors = FALSE),
            .bc_cn1_verdict(d$accuracy, d$base_rate[1]))
    }))
    rownames(verdict) <- NULL
    out$cn1 <- perm
    out$cn1_verdict <- verdict
    utils::write.csv(perm, file.path(output_dir, "TABLE_barcoding_cn1_permutations.csv"), row.names = FALSE)
    utils::write.csv(verdict, file.path(output_dir, "TABLE_barcoding_cn1_verdict.csv"), row.names = FALSE)
  }

  if ("CN2" %in% controls && method == "idtaxa") {
    queries <- list()
    for (l in loci) {
      f <- if (is.null(outgroup_dir)) NA_character_ else file.path(outgroup_dir, paste0(l, ".fasta"))
      if (is.na(f) || !file.exists(f)) next
      queries[[length(queries) + 1L]] <- .bc_cn2_queries_idtaxa(
        ape::read.dna(f, format = "fasta"), sequences[[l]], lib[lib$locus == l, , drop = FALSE], l,
        seed = seed, threshold = threshold)
    }
    cn2_q <- if (length(queries) > 0L) do.call(rbind, queries) else NULL
    if (!is.null(cn2_q)) rownames(cn2_q) <- NULL
    out$cn2_queries_idtaxa <- cn2_q
    utils::write.csv(if (is.null(cn2_q)) data.frame() else cn2_q,
                     file.path(output_dir, "TABLE_barcoding_cn2_queries_idtaxa.csv"), row.names = FALSE)
  }

  if ("CN2" %in% controls && method == "nn") {
    queries <- list()
    summary_df <- list()
    for (l in loci) {
      f <- if (is.null(outgroup_dir)) NA_character_ else file.path(outgroup_dir, paste0(l, ".fasta"))
      if (is.na(f) || !file.exists(f)) {
        # A locus with no outgroup sequence still appears, and says so. A control over zero
        # sequences is not a control, and hiding the row would hide exactly that.
        summary_df[[length(summary_df) + 1L]] <- data.frame(
          locus = l, queries = 0L, comparable = 0L, reversed = 0L, no_match = 0L,
          n_species = 0L, state_1 = 0L, state_2 = 0L, state_3 = 0L,
          prop_unassigned = NA_real_, median_nn_distance = NA_real_,
          reason = "no outgroup sequences", stringsAsFactors = FALSE)
        next
      }
      q <- .bc_cn2_queries(ape::read.dna(f, format = "fasta"), sequences[[l]],
                           lib[lib$locus == l, , drop = FALSE], model, min_comparable, l,
                           mafft_exec = mafft_exec, mafft_opts = mafft_opts)
      queries[[length(queries) + 1L]] <- q
      # A sequence that matches the locus in neither direction is published and counted, and stays
      # out of everything else: it is not a query of this control. Measured on 2026-09-25, all 83
      # outgroup sequences of trnS-trnG are in that state, so that locus tests nothing.
      comp <- q[q$orientation != "no_match", , drop = FALSE]
      summary_df[[length(summary_df) + 1L]] <- data.frame(
        locus = l, queries = nrow(q), comparable = nrow(comp),
        reversed = as.integer(sum(q$orientation == "reverse")),
        no_match = as.integer(sum(q$orientation == "no_match")),
        # The size of this control is the number of species, not the number of sequences. GenBank
        # holds population level studies: in the run of 2026-09-25 the 83 outgroup sequences of
        # trnS-trnG were all Talinopsis frutescens, one species queried 83 times. Publishing only
        # the sequences would make that look like the best controlled locus instead of the worst.
        n_species = length(unique(comp$query_species)),
        state_1 = as.integer(sum(comp$state == 1L)), state_2 = as.integer(sum(comp$state == 2L)),
        state_3 = as.integer(sum(comp$state == 3L)),
        prop_unassigned = if (nrow(comp) == 0L) NA_real_ else mean(comp$state == 3L),
        median_nn_distance = if (nrow(comp) == 0L) NA_real_ else
          stats::median(comp$nn_distance, na.rm = TRUE),
        reason = NA_character_, stringsAsFactors = FALSE)
    }
    cn2 <- do.call(rbind, summary_df)
    rownames(cn2) <- NULL
    cn2_q <- if (length(queries) > 0L) do.call(rbind, queries) else
      .bc_cn2_queries_empty()
    rownames(cn2_q) <- NULL
    out$cn2 <- cn2
    out$cn2_queries <- cn2_q
    utils::write.csv(cn2, file.path(output_dir, "TABLE_barcoding_cn2_outgroup.csv"), row.names = FALSE)
    utils::write.csv(cn2_q, file.path(output_dir, "TABLE_barcoding_cn2_queries.csv"), row.names = FALSE)
  }

  if ("CN3" %in% controls) {
    rows <- list()
    for (sc in schemes) {
      tab <- utils::read.csv(file.path(folds_dir, paste0("TABLE_barcoding_folds_", sc, ".csv")),
                             stringsAsFactors = FALSE)
      for (l in loci) {
        d <- lib[lib$locus == l, , drop = FALSE]
        species <- stats::setNames(d$species, d$sid)
        folds <- .bc_sample_folds(.bc_folds_from_table(tab[tab$locus == l, , drop = FALSE], d$sid),
                                     max_folds, seed)
        if (length(folds) == 0L) next
        pred <- do.call(rbind, lapply(folds, function(p) {
          do.call(rbind, lapply(p$test_ids, function(q) {
            # The training set of the honest run with the query put back into it. One difference,
            # and only one, between this figure and the honest one.
            r <- .bc_classify_nn(dist_matrices[[l]], c(p$train_ids, q), q, species,
                                 allow_resubstitution = TRUE)
            cbind(data.frame(locus = l, scheme = sc, method = method,
                             true_species = unname(species[q]),
                             true_genus = .bc_genus(unname(species[q])),
                             stringsAsFactors = FALSE), r)
          }))
        }))
        rows[[length(rows) + 1L]] <- .bc_metric_accuracy(pred)
      }
    }
    cn3 <- do.call(rbind, rows)
    cn3$label <- "irreproducible_resubstitution"
    rownames(cn3) <- NULL
    out$cn3 <- cn3
    utils::write.csv(cn3, file.path(output_dir, "TABLE_barcoding_cn3_resubstitution.csv"), row.names = FALSE)
    message("CN3 written with its label. Its figure is published in Phase 6, and never without the honest one beside it.")
  }

  .bc_controls_banner(out, output_dir, method)
  invisible(out)
}

#' One row per outgroup query, with the distance the threshold of Phase 6 will be swept over
#'
#' The outgroup sequences come from the branch's own assembly, so they share the locus but not the
#' alignment. Positions that one of the two does not have are not comparable, which is the same
#' rule the rest of the branch uses (E13).
#'
#' The state alone does not measure this control. `.bc_classify_nn()` has no rule of remoteness by
#' design, so an alien query still gets its nearest neighbour and can come back as state 1 with a
#' species that means nothing. What the control measures is **how far** these queries fall: if they
#' do not fall beyond the library, no threshold will ever separate them, and that is what sec. 4 of
#' the plan is asking. Hence the distance, published per query.
#' @noRd
.bc_assert_mafft <- function(mafft_exec) {
  ok <- tryCatch(length(check_mafft_available(mafft_exec)) > 0L,
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!isTRUE(ok)) {
    stop("CN2 aligns every outgroup query against the library and needs MAFFT. ",
         "It is not available as '", mafft_exec, "'.", call. = FALSE)
  }
  invisible(TRUE)
}

# Orientation of a query against the library, with the rule of .normalise_strand().
#
# GenBank stores a record on whichever strand the submitter deposited, and MAFFT --add compares only
# the orientation it is given. In the run of 2026-09-25 four of the thirteen rbcL sequences of the
# outgroup, MH767707.1 to MH767710.1, were deposited reversed: they came out at distance 0.48 while
# their conspecifics sat at 0.011, and nothing in the table said why.
#
# The reference pool is built from the library alone and handed in already built, which is what
# keeps one query from voting on the direction of another and what keeps the cost to one pass per
# locus. A query matching the locus in neither direction is left exactly as it arrived and reported,
# the same thing .normalise_strand() does: that is a homology problem, not an orientation problem.
# A query too short to hold five k-mers cannot be oriented either, and falls in the same bucket.
.bc_strand_pool <- function(lib, k = 20L, min_share = 0.15) {
  m <- as.character(as.matrix(lib))
  seqs <- gsub("-", "", toupper(apply(m, 1, paste, collapse = "")), fixed = TRUE)
  sets <- .kmer_sets(seqs, k)
  usable <- .strand_usable(sets)
  if (!any(usable)) return(character(0))
  .strand_reference_pool(sets, nchar(seqs), usable, min_share)
}

.bc_orient_to_library <- function(query, pool, k = 20L, min_share = 0.15, ratio = 3) {
  m <- as.matrix(query)
  seq_name <- rownames(m)
  s <- gsub("-", "", toupper(paste(as.character(m), collapse = "")), fixed = TRUE)
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(s)))
  sets <- .kmer_sets(c(s, rc), k)
  usable <- .strand_usable(sets)[1]

  fwd <- .kmer_share(sets[[1]], pool)
  rev <- .kmer_share(sets[[2]], pool)
  action <- .strand_action(fwd, rev, usable, min_share, ratio)

  orientation <- switch(action, kept = "forward", reverse_complemented = "reverse",
                        "no_match")
  sequence <- if (action == "reverse_complemented") rc else s
  list(query = .bc_dnabin_row(sequence, seq_name), orientation = orientation,
       share_forward = fwd, share_reverse = rev)
}

.bc_dnabin_row <- function(s, seq_name) {
  ape::as.DNAbin(matrix(strsplit(tolower(s), "")[[1]], nrow = 1,
                        dimnames = list(seq_name, NULL)))
}

# Adds one query to the alignment of the library without moving a single column of it.
#
# Until 2026-09-25 CN2 truncated both alignments to their smallest width and compared them position
# by position. That assumes column i of one alignment is homologous to column i of another, which is
# false: in the real run the library of matK had 2464 columns and the outgroup sequences ran from
# 326 to 1533 bases, each starting at a different place. The defect showed up as a crash only
# because the sequences of step 1 are not aligned, and the crash was the lucky part.
#
# MAFFT --add places the query against the profile of the library, and --keeplength forbids it from
# opening new columns: the insertions of the query are discarded. That is what keeps the distance of
# one query comparable with the distance of another, and it is the same operation the identification
# step will need when a user brings a sequence of their own.
#
# run_mafft() splits its options by whitespace before pasting the input file, so a path carrying a
# space would break the call. The work happens inside a temporary directory with relative names, so
# no path of the caller ever reaches the command line.
.bc_align_to_library <- function(lib, query, mafft_exec = "mafft", mafft_opts = "--auto") {
  lib_m <- as.matrix(lib)
  qn <- labels(query)
  if (length(qn) != 1L) stop("The query has to be a single sequence, not ", length(qn), ".", call. = FALSE)

  dir <- tempfile("bc_align_")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  previous_wd <- setwd(dir)
  on.exit({ setwd(previous_wd); unlink(dir, recursive = TRUE) }, add = TRUE)

  ape::write.FASTA(lib_m, "library.fasta")
  ape::write.FASTA(query, "query.fasta")
  run_mafft("library.fasta", "aligned.fasta", mafft_exec = mafft_exec,
            mafft_opts = paste("--add query.fasta --keeplength", mafft_opts))
  al <- ape::read.dna("aligned.fasta", format = "fasta", as.matrix = TRUE)

  if (ncol(al) != ncol(lib_m)) {
    stop("MAFFT --keeplength returned ", ncol(al), " columns for an alignment of ", ncol(lib_m),
         ". The library would stop being fixed between queries.", call. = FALSE)
  }
  missing <- setdiff(c(rownames(lib_m), qn), rownames(al))
  if (length(missing) > 0L) {
    stop("MAFFT did not return ", length(missing), " of the sequences it was given, the first being '",
         missing[1], "'.", call. = FALSE)
  }
  al <- al[c(rownames(lib_m), qn), , drop = FALSE]
  rownames(al) <- c(rownames(lib_m), "query")
  al
}

.bc_cn2_queries <- function(og, lib_dna, lib_tab, model, min_comparable, locus,
                            mafft_exec = "mafft", mafft_opts = "--auto") {
  species <- stats::setNames(lib_tab$species, lib_tab$sid)
  h <- .bc_parse_header(labels(og))
  lib_m <- as.matrix(lib_dna)
  # The pool of the locus is built once; the alignment is one call per query, decision of BMM of
  # 2026-09-25, so that no query is ever aligned in the company of another, which is how a real
  # query will arrive.
  #
  # Since 2026-09-26 the orientation, the alignment and the classification of one query are
  # .bc_classify_by_add(), shared with the add path of step 7, and a query that matches the locus
  # in neither direction comes back as state 3 with no species and no distance (decision D4).
  lib_m <- lib_m[lib_tab$sid, , drop = FALSE]
  pool <- .bc_strand_pool(lib_m)
  queries <- as.character(og)
  if (!is.list(queries)) queries <- lapply(seq_len(nrow(queries)), function(i) queries[i, ])
  rows <- lapply(seq_along(queries), function(i) {
    r <- .bc_classify_by_add(lib_m, species, paste(queries[[i]], collapse = ""), model,
                             min_comparable, pool = pool, mafft_exec = mafft_exec,
                             mafft_opts = mafft_opts)
    cbind(data.frame(locus = locus, sid = h$sid[i], query_species = h$species[i],
                     stringsAsFactors = FALSE), r)
  })
  do.call(rbind, rows)
}

#' CN2 with IdTaxa: the outgroup queries of one locus, oriented and classified
#'
#' The orientation is the one of the nearest neighbour, against the strand pool of the library of
#' the locus; a query that matches the locus in neither direction is state 3 with the reason
#' `no_match` and no confidence, and IdTaxa is not called for it (decision D4). The others are
#' classified by IdTaxa trained once on the whole locus, at threshold 0, each with the seed of
#' `.bc_query_seed(seed, locus, "cn2", 0, sid)`.
#' @noRd
.bc_cn2_queries_idtaxa <- function(og, lib_dna, lib_tab, locus, seed = 1L, threshold = 60) {
  restore_rng <- .bc_rng_state()
  on.exit(restore_rng(), add = TRUE)
  species <- stats::setNames(lib_tab$species, lib_tab$sid)
  h <- .bc_parse_header(labels(og))
  lib_m <- as.matrix(lib_dna)[lib_tab$sid, , drop = FALSE]
  pool <- .bc_strand_pool(lib_m)
  lib_seqs <- gsub("-", "", toupper(apply(as.character(lib_m), 1, paste, collapse = "")), fixed = TRUE)
  names(lib_seqs) <- rownames(lib_m)
  trained <- .bc_idtaxa_train(lib_seqs, species[names(lib_seqs)],
                              train_seed = .bc_query_seed(seed, locus, "cn2", 0L, "<training>"))
  queries <- as.character(og)
  if (!is.list(queries)) queries <- lapply(seq_len(nrow(queries)), function(i) queries[i, ])
  rows <- lapply(seq_along(queries), function(i) {
    s <- gsub("-", "", paste(queries[[i]], collapse = ""), fixed = TRUE)
    o <- .bc_orient_to_library(.bc_dnabin_row(s, "query"), pool)
    r <- if (o$orientation == "no_match") {
      cbind(.bc_prediction_row(3L, reason = "no_match"),
            data.frame(genus_idtaxa = NA_character_, species_idtaxa = NA_character_,
                       genus_confidence = NA_real_, species_confidence = NA_real_,
                       stringsAsFactors = FALSE))
    } else {
      oriented <- toupper(paste(as.character(as.matrix(o$query)), collapse = ""))
      .bc_classify_idtaxa(lib_seqs, species[names(lib_seqs)], oriented, threshold = threshold,
                          trained = trained,
                          query_seed = .bc_query_seed(seed, locus, "cn2", 0L, h$sid[i]))
    }
    cbind(data.frame(locus = locus, sid = h$sid[i], query_species = h$species[i],
                     orientation = o$orientation, stringsAsFactors = FALSE), r)
  })
  do.call(rbind, rows)
}

#' The empty shape of that table, so a run with no outgroup at all still writes its columns
#' @noRd
.bc_cn2_queries_empty <- function() {
  cbind(data.frame(locus = character(0), sid = character(0), query_species = character(0),
                   orientation = character(0), stringsAsFactors = FALSE),
        .bc_prediction_row(1L)[0, , drop = FALSE])
}

#' Closing banner of step 8
#' @noRd
.bc_controls_banner <- function(out, output_dir, method) {
  cat("\n====================================================\n")
  cat("  Barcoding Negative Controls Complete \U0001f335\n")
  cat("====================================================\n")
  cat("  Method:                ", method, "\n")
  if (!is.null(out$cn1_verdict)) {
    cat("  CN1: ", nrow(out$cn1_verdict), " locus-scheme pairs, leakage flagged in ",
        sum(out$cn1_verdict$leakage), "\n", sep = "")
  }
  if (!is.null(out$cn2)) {
    cat("  CN2: ", sum(out$cn2$n_species), " outgroup species (", sum(out$cn2$comparable),
        " comparable sequences of ", sum(out$cn2$queries), ") over ",
        sum(out$cn2$comparable > 0), " loci\n", sep = "")
    cat("       ", sum(out$cn2$reversed), " reversed before aligning; ",
        sum(out$cn2$no_match), " match the locus in neither direction\n", sep = "")
  }
  if (!is.null(out$cn2_queries_idtaxa)) {
    q <- out$cn2_queries_idtaxa
    cat("  CN2 with IdTaxa: ", sum(q$orientation != "no_match"), " comparable outgroup queries of ",
        nrow(q), "; ", sum(q$orientation == "no_match"), " match the locus in neither direction\n", sep = "")
  }
  if (!is.null(out$cn3)) cat("  CN3: written and labelled; its figure belongs to Phase 6\n")
  cat("  Output directory:      ", output_dir, "\n")
  cat("====================================================\n\n")
  invisible(TRUE)
}
