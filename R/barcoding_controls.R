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
    return(data.frame(locus = character(0), esquema = character(0), metodo = character(0),
                      consultas = integer(0), aciertos = integer(0), exactitud = numeric(0),
                      tasa_base = numeric(0), stringsAsFactors = FALSE))
  }
  key <- paste(pred$locus, pred$esquema, pred$metodo, sep = "\r")
  out <- do.call(rbind, lapply(sort(unique(key), method = "radix"), function(k) {
    d <- pred[key == k, , drop = FALSE]
    especie <- d$esquema[1] == "species"
    verdad <- if (especie) d$especie_verdadera else d$genero_verdadero
    acierto <- if (especie) {
      d$estado == 1L & !is.na(d$especie_predicha) & d$especie_predicha == d$especie_verdadera
    } else {
      d$estado %in% c(1L, 2L) & !is.na(d$genero_predicho) & d$genero_predicho == d$genero_verdadero
    }
    data.frame(locus = d$locus[1], esquema = d$esquema[1], metodo = d$metodo[1],
               consultas = nrow(d), aciertos = as.integer(sum(acierto)),
               exactitud = sum(acierto) / nrow(d),
               tasa_base = max(table(verdad)) / length(verdad),
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
.bc_cn1_verdict <- function(exactitudes, tasa_base) {
  exactitudes <- exactitudes[!is.na(exactitudes)]
  desv <- if (length(exactitudes) > 1L) stats::sd(exactitudes) else 0
  umbral <- tasa_base + 3 * desv
  data.frame(permutaciones = length(exactitudes), media = mean(exactitudes), sd = desv,
             tasa_base = tasa_base, umbral = umbral, hay_fuga = mean(exactitudes) > umbral,
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
#' subtracted. Written with the label `irreproducible_resustitucion` in the table itself, so the
#' figure cannot travel without it. **Its value is published in Phase 6 and only next to the honest
#' figure**, which does not exist yet.
#'
#' @param library_dir Character. Output directory of [finalize_barcoding_library()].
#' @param folds_dir Character. Output directory of [build_barcoding_folds()].
#' @param output_dir Character. Destination; refused inside a directory of the phylogeny.
#' @param outgroup_dir Character or `NULL`. Directory with `<locus>.fasta` of the outgroup, as the
#'   assembly writes them. `NULL` skips CN2 and says so.
#' @param method Character. `"nn"` or `"idtaxa"`.
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
                                   mafft_opts = "--auto") {
  method <- match.arg(method)
  .bc_assert_output_dir(output_dir)
  lib <- .bc_library_from_dir(library_dir)
  for (sc in schemes) {
    f <- file.path(folds_dir, paste0("TABLE_barcoding_folds_", sc, ".csv"))
    if (!file.exists(f)) {
      stop("No ", basename(f), " in ", folds_dir, ". Run build_barcoding_folds() first.", call. = FALSE)
    }
  }
  if ("CN2" %in% controls && !is.null(outgroup_dir)) {
    # CN2 aligns every query against the library, so the controls need MAFFT at run time. The
    # check belongs here and not inside the CN2 loop: an absent binary has to cost seconds, not
    # the hour that CN1 takes before reaching it.
    .bc_assert_mafft(mafft_exec)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  loci_todos <- sort(unique(lib$locus), method = "radix")
  loci <- if (is.null(loci)) loci_todos else intersect(loci_todos, loci)

  matrices <- list()
  secuencias <- list()
  for (l in loci) {
    dna <- ape::read.dna(file.path(library_dir, paste0("LIB_", l, ".fasta")), format = "fasta")
    rownames(dna) <- .bc_parse_header(labels(dna))$sid
    matrices[[l]] <- .bc_classifier_matrix(dna, model, min_comparable)
    secuencias[[l]] <- dna
  }

  out <- list()

  if ("CN1" %in% controls) {
    filas <- list()
    for (sc in schemes) {
      tab <- utils::read.csv(file.path(folds_dir, paste0("TABLE_barcoding_folds_", sc, ".csv")),
                             stringsAsFactors = FALSE)
      for (l in loci) {
        d <- lib[lib$locus == l, , drop = FALSE]
        species <- stats::setNames(d$species, d$sid)
        pliegues <- .bc_sample_folds(.bc_folds_from_table(tab[tab$locus == l, , drop = FALSE], d$sid),
                                     max_folds, seed)
        if (length(pliegues) == 0L) next
        for (k in seq_len(permutations)) {
          pred <- do.call(rbind, lapply(pliegues, function(p) {
            etiquetas <- .barcoding_permute_labels(
              list(train_ids = p$train_ids, test_ids = p$test_ids), species, seed = seed + k)
            do.call(rbind, lapply(p$test_ids, function(q) {
              r <- .bc_classify_nn(matrices[[l]], p$train_ids, q, etiquetas)
              cbind(data.frame(locus = l, esquema = sc, metodo = method,
                               especie_verdadera = unname(species[q]),
                               genero_verdadero = .bc_genus(unname(species[q])),
                               stringsAsFactors = FALSE), r)
            }))
          }))
          m <- .bc_metric_accuracy(pred)
          filas[[length(filas) + 1L]] <- cbind(data.frame(permutacion = k, stringsAsFactors = FALSE), m)
        }
        message(sprintf("CN1, locus '%s', scheme %s: %d permutations over %d folds.",
                        l, sc, permutations, length(pliegues)))
      }
    }
    perm <- do.call(rbind, filas)
    rownames(perm) <- NULL
    key <- paste(perm$locus, perm$esquema, perm$metodo, sep = "\r")
    veredicto <- do.call(rbind, lapply(sort(unique(key), method = "radix"), function(k) {
      d <- perm[key == k, , drop = FALSE]
      cbind(data.frame(locus = d$locus[1], esquema = d$esquema[1], metodo = d$metodo[1],
                       stringsAsFactors = FALSE),
            .bc_cn1_verdict(d$exactitud, d$tasa_base[1]))
    }))
    rownames(veredicto) <- NULL
    out$cn1 <- perm
    out$cn1_verdict <- veredicto
    utils::write.csv(perm, file.path(output_dir, "TABLE_barcoding_cn1_permutations.csv"), row.names = FALSE)
    utils::write.csv(veredicto, file.path(output_dir, "TABLE_barcoding_cn1_verdict.csv"), row.names = FALSE)
  }

  if ("CN2" %in% controls) {
    consultas <- list()
    resumen <- list()
    for (l in loci) {
      f <- if (is.null(outgroup_dir)) NA_character_ else file.path(outgroup_dir, paste0(l, ".fasta"))
      if (is.na(f) || !file.exists(f)) {
        # A locus with no outgroup sequence still appears, and says so. A control over zero
        # sequences is not a control, and hiding the row would hide exactly that.
        resumen[[length(resumen) + 1L]] <- data.frame(
          locus = l, consultas = 0L, comparables = 0L, invertidas = 0L, sin_coincidencia = 0L,
          especies = 0L, estado_1 = 0L, estado_2 = 0L, estado_3 = 0L,
          prop_no_asignable = NA_real_, mediana_distancia_vecino = NA_real_,
          motivo = "sin secuencias de grupo externo", stringsAsFactors = FALSE)
        next
      }
      q <- .bc_cn2_queries(ape::read.dna(f, format = "fasta"), secuencias[[l]],
                           lib[lib$locus == l, , drop = FALSE], model, min_comparable, l,
                           mafft_exec = mafft_exec, mafft_opts = mafft_opts)
      consultas[[length(consultas) + 1L]] <- q
      # A sequence that matches the locus in neither direction is published and counted, and stays
      # out of everything else: it is not a query of this control. Measured on 2026-09-25, all 83
      # outgroup sequences of trnS-trnG are in that state, so that locus tests nothing.
      comp <- q[q$orientacion != "sin_coincidencia", , drop = FALSE]
      resumen[[length(resumen) + 1L]] <- data.frame(
        locus = l, consultas = nrow(q), comparables = nrow(comp),
        invertidas = as.integer(sum(q$orientacion == "reversa")),
        sin_coincidencia = as.integer(sum(q$orientacion == "sin_coincidencia")),
        # The size of this control is the number of species, not the number of sequences. GenBank
        # holds population level studies: in the run of 2026-09-25 the 83 outgroup sequences of
        # trnS-trnG were all Talinopsis frutescens, one species queried 83 times. Publishing only
        # the sequences would make that look like the best controlled locus instead of the worst.
        especies = length(unique(comp$especie_consulta)),
        estado_1 = as.integer(sum(comp$estado == 1L)), estado_2 = as.integer(sum(comp$estado == 2L)),
        estado_3 = as.integer(sum(comp$estado == 3L)),
        prop_no_asignable = if (nrow(comp) == 0L) NA_real_ else mean(comp$estado == 3L),
        mediana_distancia_vecino = if (nrow(comp) == 0L) NA_real_ else
          stats::median(comp$distancia_vecino, na.rm = TRUE),
        motivo = NA_character_, stringsAsFactors = FALSE)
    }
    cn2 <- do.call(rbind, resumen)
    rownames(cn2) <- NULL
    cn2_q <- if (length(consultas) > 0L) do.call(rbind, consultas) else
      .bc_cn2_queries_empty()
    rownames(cn2_q) <- NULL
    out$cn2 <- cn2
    out$cn2_queries <- cn2_q
    utils::write.csv(cn2, file.path(output_dir, "TABLE_barcoding_cn2_outgroup.csv"), row.names = FALSE)
    utils::write.csv(cn2_q, file.path(output_dir, "TABLE_barcoding_cn2_queries.csv"), row.names = FALSE)
  }

  if ("CN3" %in% controls) {
    filas <- list()
    for (sc in schemes) {
      tab <- utils::read.csv(file.path(folds_dir, paste0("TABLE_barcoding_folds_", sc, ".csv")),
                             stringsAsFactors = FALSE)
      for (l in loci) {
        d <- lib[lib$locus == l, , drop = FALSE]
        species <- stats::setNames(d$species, d$sid)
        pliegues <- .bc_sample_folds(.bc_folds_from_table(tab[tab$locus == l, , drop = FALSE], d$sid),
                                     max_folds, seed)
        if (length(pliegues) == 0L) next
        pred <- do.call(rbind, lapply(pliegues, function(p) {
          do.call(rbind, lapply(p$test_ids, function(q) {
            # The training set of the honest run with the query put back into it. One difference,
            # and only one, between this figure and the honest one.
            r <- .bc_classify_nn(matrices[[l]], c(p$train_ids, q), q, species,
                                 allow_resubstitution = TRUE)
            cbind(data.frame(locus = l, esquema = sc, metodo = method,
                             especie_verdadera = unname(species[q]),
                             genero_verdadero = .bc_genus(unname(species[q])),
                             stringsAsFactors = FALSE), r)
          }))
        }))
        filas[[length(filas) + 1L]] <- .bc_metric_accuracy(pred)
      }
    }
    cn3 <- do.call(rbind, filas)
    cn3$etiqueta <- "irreproducible_resustitucion"
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
  nombre <- rownames(m)
  s <- gsub("-", "", toupper(paste(as.character(m), collapse = "")), fixed = TRUE)
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(s)))
  sets <- .kmer_sets(c(s, rc), k)
  usable <- .strand_usable(sets)[1]

  fwd <- .kmer_share(sets[[1]], pool)
  rev <- .kmer_share(sets[[2]], pool)
  accion <- .strand_action(fwd, rev, usable, min_share, ratio)

  orientacion <- switch(accion, kept = "directa", reverse_complemented = "reversa",
                        "sin_coincidencia")
  secuencia <- if (accion == "reverse_complemented") rc else s
  list(query = .bc_dnabin_row(secuencia, nombre), orientacion = orientacion,
       share_directa = fwd, share_reversa = rev)
}

.bc_dnabin_row <- function(s, nombre) {
  ape::as.DNAbin(matrix(strsplit(tolower(s), "")[[1]], nrow = 1,
                        dimnames = list(nombre, NULL)))
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
  anterior <- setwd(dir)
  on.exit({ setwd(anterior); unlink(dir, recursive = TRUE) }, add = TRUE)

  ape::write.FASTA(lib_m, "biblioteca.fasta")
  ape::write.FASTA(query, "consulta.fasta")
  run_mafft("biblioteca.fasta", "alineado.fasta", mafft_exec = mafft_exec,
            mafft_opts = paste("--add consulta.fasta --keeplength", mafft_opts))
  al <- ape::read.dna("alineado.fasta", format = "fasta", as.matrix = TRUE)

  if (ncol(al) != ncol(lib_m)) {
    stop("MAFFT --keeplength returned ", ncol(al), " columns for an alignment of ", ncol(lib_m),
         ". The library would stop being fixed between queries.", call. = FALSE)
  }
  faltan <- setdiff(c(rownames(lib_m), qn), rownames(al))
  if (length(faltan) > 0L) {
    stop("MAFFT did not return ", length(faltan), " of the sequences it was given, the first being '",
         faltan[1], "'.", call. = FALSE)
  }
  al <- al[c(rownames(lib_m), qn), , drop = FALSE]
  rownames(al) <- c(rownames(lib_m), "consulta")
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
  pool <- .bc_strand_pool(lib_m)
  consultas <- as.list(og)
  filas <- lapply(seq_along(consultas), function(i) {
    q <- consultas[i]
    names(q) <- "consulta"
    o <- .bc_orient_to_library(q, pool)
    junto <- .bc_align_to_library(lib_m, o$query, mafft_exec = mafft_exec, mafft_opts = mafft_opts)
    dm <- .bc_classifier_matrix(junto, model, min_comparable)
    sp <- c(species, consulta = "consulta_externa")
    r <- .bc_classify_nn(dm, lib_tab$sid, "consulta", sp)
    cbind(data.frame(locus = locus, sid = h$sid[i], especie_consulta = h$species[i],
                     orientacion = o$orientacion, stringsAsFactors = FALSE), r)
  })
  do.call(rbind, filas)
}

#' The empty shape of that table, so a run with no outgroup at all still writes its columns
#' @noRd
.bc_cn2_queries_empty <- function() {
  cbind(data.frame(locus = character(0), sid = character(0), especie_consulta = character(0),
                   orientacion = character(0), stringsAsFactors = FALSE),
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
        sum(out$cn1_verdict$hay_fuga), "\n", sep = "")
  }
  if (!is.null(out$cn2)) {
    cat("  CN2: ", sum(out$cn2$especies), " outgroup species (", sum(out$cn2$comparables),
        " comparable sequences of ", sum(out$cn2$consultas), ") over ",
        sum(out$cn2$comparables > 0), " loci\n", sep = "")
    cat("       ", sum(out$cn2$invertidas), " reversed before aligning; ",
        sum(out$cn2$sin_coincidencia), " match the locus in neither direction\n", sep = "")
  }
  if (!is.null(out$cn3)) cat("  CN3: written and labelled; its figure belongs to Phase 6\n")
  cat("  Output directory:      ", output_dir, "\n")
  cat("====================================================\n\n")
  invisible(TRUE)
}
