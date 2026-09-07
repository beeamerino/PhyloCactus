# Regression tests for the marker-set defects audited on 2026-09-01.
#
# Three filters decided the composition of the supermatrix using the ingroup alone, and the
# outgroup was joined afterwards. The August 2026 matrix reached RAxML-NG with the branch
# subtending its outgroup carrying 0.03 expected substitutions, against roughly 600 in the two
# preceding runs, and every calibrated node of the resulting chronogram returned its own bound
# instead of an estimate. These tests hold the three guards that make that failure visible.

write_fasta <- function(seqs, path) {
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(seqs), path)
  invisible(path)
}

# A locus with enough signal to survive every criterion except the sequence count, so the count
# is the only thing under test.
informative_alignment <- function(labels, len = 400L, seed = 1L) {
  set.seed(seed)
  vapply(seq_along(labels), function(i) {
    paste(sample(c("A", "C", "G", "T"), len, replace = TRUE), collapse = "")
  }, character(1)) |> stats::setNames(labels)
}

test_that("the outgroup counts are reported but never decisive", {
  skip_if_not_installed("Biostrings")

  tmp <- withr::local_tempdir()
  ing <- file.path(tmp, "ingroup");  dir.create(ing)
  out <- file.path(tmp, "outgroup"); dir.create(out)

  # 8 ingroup and 6 outgroup sequences against a threshold of 12. The joint total would clear it;
  # the ingroup count does not, and the ingroup count is what decides. This module asks whether a
  # locus resolves the ingroup radiation, and that cannot depend on how many outgroup accessions
  # exist. Bringing such a locus back is integrate_and_clean_markers(readmit_markers = ).
  write_fasta(informative_alignment(paste0("Ingroup_sp", 1:8)),
              file.path(ing, "ALN_masked_final_shared.fasta"))
  write_fasta(informative_alignment(paste0("Portulaca_sp", 1:6), seed = 2L),
              file.path(out, "ALN_masked_final_shared.fasta"))

  # Random fixture sequences carry no phylogenetic structure, so phangorn returns NaN distances
  # and ggplot drops the corresponding points. That noise belongs to the fixture.
  res <- suppressWarnings(run_marker_screening(
    fasta_folder = ing, out_base = file.path(tmp, "res"),
    min_nseq_to_retain = 12L, min_aln_len_to_retain = 100L,
    outgroup_folder = out
  ))

  expect_equal(res$n_ingroup[1], 8L)
  expect_equal(res$n_outgroup[1], 6L)
  expect_equal(res$n_total[1], 14L)
  expect_equal(res$decision[1], "NO")
  expect_true(grepl("nseq_lt_12", res$decision_reason[1]))
})

test_that("report_marker_group_coverage names the loci with no outgroup coverage", {
  skip_if_not_installed("Biostrings")

  tmp <- withr::local_tempdir()

  # `nuclear_like` reproduces the phyC pattern: dense in one clade, absent in the outgroup.
  gapped <- paste(rep("-", 200), collapse = "")
  real   <- paste(rep("A", 200), collapse = "")
  write_fasta(stats::setNames(c(real, real, gapped, gapped),
                              c("Opuntia_a", "Opuntia_b", "Portulaca_a", "Portulaca_b")),
              file.path(tmp, "nuclear_like.fasta"))
  write_fasta(stats::setNames(rep(real, 4),
                              c("Opuntia_a", "Opuntia_b", "Portulaca_a", "Portulaca_b")),
              file.path(tmp, "shared.fasta"))

  expect_warning(
    cov <- report_marker_group_coverage(tmp, outgroup_pattern = "^Portulaca_"),
    "nuclear_like"
  )
  expect_equal(cov$n_outgroup_covered[cov$marker == "nuclear_like"], 0L)
  expect_equal(cov$n_outgroup_covered[cov$marker == "shared"], 2L)
  expect_equal(cov$n_ingroup_covered[cov$marker == "nuclear_like"], 2L)
})

test_that("rooting_pattern reports without exempting, and protect_markers scopes the exemption", {
  skip_if_not_installed("Biostrings")
  skip_if(Sys.which("mafft") == "", "requires MAFFT")

  tmp <- withr::local_tempdir()
  ind <- file.path(tmp, "in"); dir.create(ind)

  # The fixture has to make the filter fire, and that depends on the alignment width DECIPHER
  # leaves behind. With few sequences, a single short terminal puts a gap in a quarter of every
  # column it does not cover, MaskAlignment removes those columns, and coercing back to a
  # DNAStringSet drops them: aln_len collapses to the stretch the fragment covers and its
  # occupancy becomes 1. Twenty near-identical ingroup sequences keep that share below 5%, so
  # the columns survive and the 40 bp terminal is genuinely under 0.80 of the width.
  set.seed(11L)
  base_seq <- paste(sample(c("A", "C", "G", "T"), 800L, replace = TRUE), collapse = "")
  mutate <- function(x, n) {
    ch <- strsplit(x, "", fixed = TRUE)[[1]]
    pos <- sample(seq_along(ch), n)
    ch[pos] <- sample(c("A", "C", "G", "T"), n, replace = TRUE)
    paste(ch, collapse = "")
  }
  ingroup_seqs <- vapply(seq_len(20L), function(i) mutate(base_seq, 30L), character(1))
  short_seq    <- substr(base_seq, 1L, 40L)

  for (m in c("trnL_trnF", "matK")) {
    write_fasta(stats::setNames(c(ingroup_seqs, short_seq),
                                c(paste0("Opuntia_sp", seq_len(20L)), "Portulaca_a")),
                file.path(ind, paste0(m, ".fasta")))
  }

  run_one <- function(outdir, ...) {
    suppressMessages(run_joint_realignment(
      input_dir = ind,
      output_fasta_dir = outdir,
      output_aln_dir = file.path(outdir, "aligned_markers"),
      min_non_gap_fraction = 0.80,
      rooting_pattern = "^Portulaca_",
      ...
    ))
    lapply(stats::setNames(c("trnL_trnF", "matK"), c("trnL_trnF", "matK")), function(m) {
      utils::read.csv(file.path(outdir, paste0("LOG_SEQ_FILTER_", m, ".csv")))
    })
  }

  # Without protect_pattern the short rooting terminal is dropped in both loci, and each drop is
  # named in a warning. The ingroup terminals must survive, or the fixture is testing nothing.
  expect_warning(logs_none <- run_one(file.path(tmp, "none")), "rooting terminal")
  for (m in c("trnL_trnF", "matK")) {
    lg <- logs_none[[m]]
    # Precondition. If masking collapsed the alignment the threshold would fall below 40 bp and
    # the test would silently stop exercising the filter, which is how the first version of this
    # fixture passed nothing while appearing to.
    expect_gt(max(lg$NonGaps), 400L)
    expect_true(all(lg$Retained[startsWith(lg$Seq, "Opuntia_")]))
    expect_false(lg$Retained[lg$Seq == "Portulaca_a"])
    expect_true(lg$IsRooting[lg$Seq == "Portulaca_a"])
    expect_true(lg$DroppedRooting[lg$Seq == "Portulaca_a"])
  }

  # Scoped to trnL_trnF, the exemption applies there and nowhere else.
  logs_scoped <- suppressWarnings(run_one(
    file.path(tmp, "scoped"),
    protect_pattern = "^Portulaca_",
    protect_markers = "trnL_trnF"
  ))
  kept <- vapply(logs_scoped, function(lg) isTRUE(lg$Retained[lg$Seq == "Portulaca_a"]), logical(1))
  expect_true(kept[["trnL_trnF"]])
  expect_false(kept[["matK"]])
  expect_true(logs_scoped$trnL_trnF$ProtectedHere[logs_scoped$trnL_trnF$Seq == "Portulaca_a"])
  expect_false(logs_scoped$matK$ProtectedHere[logs_scoped$matK$Seq == "Portulaca_a"])
})

# Added after the 2026-09-01 audit of target_genes.txt and genes_map.csv. The ingroup pepC
# accessions are "phosphoenolpyruvate carboxylase-like" (Hernandez-Hernandez, ~489 bp) and the
# outgroup ones are "phosphoenolpyruvate carboxylase (ppc-2)" (Ogburn & Edwards, ~760 bp): two
# paralogues of the PEPC gene family. The map sent both to `pepC`, once by an explicit rule and
# once because the token "phosphoenolpyruvate carboxylase" matches as a substring of
# "...carboxylase-like" and the qualifier disappeared. MAFFT then aligned them into one block, so
# 506 columns of the supermatrix carried outgroup characters with no positional homology, on the
# seven terminals that define the root. The two sets share 0.018 of their 10-mers; phyC, also
# nuclear, shares 0.715.

read_genes_map <- function() {
  path <- system.file("extdata", "genes_map.csv", package = "PhyloCactus")
  skip_if(!nzchar(path), "genes_map.csv not installed")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

test_that("the two PEPC paralogues are kept apart", {
  gm <- read_genes_map()
  lookup <- stats::setNames(trimws(gm$replace), trimws(gm$search))

  expect_equal(unname(lookup[["phosphoenolpyruvate carboxylase"]]), "pepC_like")
  expect_equal(unname(lookup[["phosphoenolpyruvate carboxylase-like"]]), "pepC_like")
  expect_equal(unname(lookup[["ppc-2"]]), "ppc2")
  expect_equal(unname(lookup[["phosphoenolpyruvate carboxylase, ppc-2"]]), "ppc2")

  # No rule may route a ppc-2 record to the ingroup paralogue again.
  ppc2_rules <- lookup[grepl("ppc-2", names(lookup), fixed = TRUE)]
  expect_true(all(ppc2_rules == "ppc2"))
  expect_false(any(ppc2_rules == "pepC_like"))
})

test_that("the qualifier that distinguishes the paralogues is extractable", {
  path <- system.file("extdata", "target_genes.txt", package = "PhyloCactus")
  skip_if(!nzchar(path), "target_genes.txt not installed")
  markers <- unique(trimws(gsub('"', "", readLines(path))))
  markers <- markers[nzchar(markers) & markers != "markers"]

  # Without the long form in the target list, the pattern matches the short one inside it and
  # "-like" is never seen. Both forms have to be declared for the alternation to tell them apart.
  expect_true("phosphoenolpyruvate carboxylase" %in% markers)
  expect_true("phosphoenolpyruvate carboxylase-like" %in% markers)
  expect_true("ppc-2" %in% markers)
})

test_that("marker_aliases maps outgroup names onto their ingroup counterpart", {
  gm <- read_genes_map()
  # Every destination that is itself a search key must map to itself, or a two-step rename would
  # depend on row order.
  lookup <- stats::setNames(trimws(gm$replace), trimws(gm$search))
  destinations <- unique(unname(lookup))
  circular <- destinations[destinations %in% names(lookup) & lookup[destinations] != destinations]
  expect_equal(length(circular), 0L,
               info = paste("destinations that are re-mapped:", paste(circular, collapse = ", ")))
})

test_that("the rbcL gene and the atpB-rbcL spacer stay separate loci", {
  gm <- read_genes_map()
  lookup <- stats::setNames(trimws(gm$replace), trimws(gm$search))

  # The two were merged on 2026-09-01 on the argument that 135 of the 138 spacer records already
  # read into the gene, and the merge was reverted the same day. Merged, the locus reached 460
  # sequences and 459 species but failed the saturation screen: slope 0.216 against a keep cutoff
  # of 0.5. The flag is not an artefact of the ragged coverage the merge introduces. Restricted to
  # the 489 columns every sequence covers completely, the slope is 0.205, against 0.931 for matK
  # and 0.646 for ITS on the same statistic. Whether that is real saturation or spurious homology
  # from aligning spacer-bearing against gene-only records, a block with that slope does not
  # belong in a matrix whose failure mode is collapsing deep branches. The merge bought 70 species
  # that no other marker had; the other 389 were already in the supermatrix.
  #
  # The distinction the map encodes is what the record names, not where the amplicon sits: a
  # record naming rbcL alone is the gene, one naming both genes spans the intergenic spacer.
  expect_equal(unname(lookup[["rbcL"]]), "rbcL")
  expect_equal(unname(lookup[["ribulose-1,5-bisphosphate carboxylase"]]), "rbcL")
  expect_equal(unname(lookup[["atpB, rbcL"]]), "atpB-rbcL")
  expect_equal(unname(lookup[["atpB, rbcL, ribulose-1,5-bisphosphate carboxylase"]]), "atpB-rbcL")

  # No rule may send a gene-only record back to the spacer locus.
  gene_only <- lookup[names(lookup) %in% c("rbcL", "ribulose-1,5-bisphosphate carboxylase",
                                           "ribulose-1,5-bisphosphate carboxylase, rbcL")]
  expect_false(any(gene_only == "atpB-rbcL"))

  # Ogburn & Edwards deposited the whole PEPC family. Each paralogue is its own locus, and none
  # of them shares a 20-mer with the ingroup locus.
  for (p in c("ppc-1E1a", "ppc-1E1d", "ppc-1E2")) {
    expect_equal(unname(lookup[[p]]), sub("-", "", p, fixed = TRUE))
  }
  expect_false(any(lookup[grepl("ppc-1E", names(lookup), fixed = TRUE)] == "pepC_like"))
})

test_that("readmission validates its arguments before touching disk", {
  skip_if_not_installed("Biostrings")
  tmp <- withr::local_tempdir()

  # The directories are empty on purpose: an incoherent argument pair has to fail on the argument,
  # not on whichever directory scan happens to run first, or the message sends you to the wrong
  # place.
  expect_error(
    integrate_and_clean_markers(
      ingroup_dir = tmp, outgroup_dir = tmp, output_dir = file.path(tmp, "out"),
      accepted_list_file = tmp, metadata_in_file = tmp, metadata_out_file = tmp,
      readmit_markers = "trnT-psbD"
    ),
    "readmit_dir"
  )
})

# Added after the third defect of the 2026-09-01 audit. Portulaca oleracea and P. pilosa deposit
# their rbcL and matK on the reverse strand. MAFFT aligned them forward, and the resulting rows sat
# at 0.51 and 0.40 observed divergence from Cactaceae where P. grandiflora, the same genus, sits at
# 0.030 and 0.066. P. oleracea was at that moment the terminal the acceptance table nominated for
# rooting the tree. The ingroup was clean: 0 of 3207 sequences reversed.

test_that(".normalise_strand flips the reversed sequences and leaves the rest alone", {
  skip_if_not_installed("Biostrings")

  set.seed(31L)
  locus <- paste(sample(c("A", "C", "G", "T"), 1200L, replace = TRUE), collapse = "")
  mutate <- function(x, n) {
    ch <- strsplit(x, "", fixed = TRUE)[[1]]
    pos <- sample(seq_along(ch), n)
    ch[pos] <- sample(c("A", "C", "G", "T"), n, replace = TRUE)
    paste(ch, collapse = "")
  }
  revcomp <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))

  forward <- vapply(1:8, function(i) mutate(locus, 60L), character(1))
  reversed <- vapply(1:2, function(i) revcomp(mutate(locus, 60L)), character(1))
  # A sequence from somewhere else entirely: it matches in neither direction and must be reported
  # as a homology problem rather than quietly flipped.
  alien <- paste(sample(c("A", "C", "G", "T"), 1200L, replace = TRUE), collapse = "")

  dna <- Biostrings::DNAStringSet(stats::setNames(
    c(forward, reversed, alien),
    c(paste0("Opuntia_sp", 1:8), "Portulaca_oleracea", "Portulaca_pilosa", "Alien_sp")))

  res <- PhyloCactus:::.normalise_strand(dna)
  act <- stats::setNames(res$log$Action, res$log$Seq)

  expect_true(all(act[paste0("Opuntia_sp", 1:8)] == "kept"))
  expect_equal(unname(act["Portulaca_oleracea"]), "reverse_complemented")
  expect_equal(unname(act["Portulaca_pilosa"]), "reverse_complemented")
  expect_equal(unname(act["Alien_sp"]), "no_match_either_direction")

  # The flip has to be the actual sequence, not just a label. After correction the two reversed
  # terminals must share k-mers with the ingroup, which is what they failed to do before.
  kset <- function(s, k = 20L) {
    m <- nchar(s) - k + 1L
    unique(substring(s, seq_len(m), seq_len(m) + k - 1L))
  }
  pool <- unique(unlist(lapply(forward, kset)))
  fixed <- as.character(res$dna)
  for (t in c("Portulaca_oleracea", "Portulaca_pilosa")) {
    after <- kset(fixed[[t]])
    expect_gt(length(intersect(after, pool)) / length(after), 0.20)
  }

  # A set already on one strand must come back untouched. A normaliser that flips on a tie would
  # corrupt every clean marker in the pipeline, which is the failure worth guarding against.
  clean <- Biostrings::DNAStringSet(stats::setNames(forward, paste0("Opuntia_sp", 1:8)))
  res2 <- PhyloCactus:::.normalise_strand(clean)
  expect_true(all(res2$log$Action == "kept"))
  expect_identical(as.character(res2$dna), as.character(clean))
})

# Added after trnT-psbD, the second false pair of the 2026-09-01 audit and the one that survived
# every count-based check. It was readmitted precisely because its counts were the best in the
# dataset: 49 ingroup terminals against 49 outgroup at 0.998 median occupancy. The sequences share
# 0.039 of their 20-mers, against 0.28-0.60 for every genuine counterpart, and sit at 0.440
# observed divergence from Cactaceae where matK gives 0.077. Counts cannot see this; only the
# sequences can.

test_that(".marker_homology_share separates a shared region from a shared name", {
  skip_if_not_installed("Biostrings")

  tmp <- withr::local_tempdir()
  set.seed(21L)
  draw <- function(n) paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
  mutate <- function(x, n) {
    ch <- strsplit(x, "", fixed = TRUE)[[1]]
    pos <- sample(seq_along(ch), n)
    ch[pos] <- sample(c("A", "C", "G", "T"), n, replace = TRUE)
    paste(ch, collapse = "")
  }
  revcomp <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))

  locus <- draw(900L)
  ing <- file.path(tmp, "in"); dir.create(ing)
  out <- file.path(tmp, "out"); dir.create(out)

  # A real counterpart: the same region, diverged. 8% of sites differ, which is more than the
  # 3.3% rbcL shows across this root.
  write_fasta(stats::setNames(vapply(1:12, function(i) mutate(locus, 72L), character(1)),
                              paste0("Opuntia_sp", 1:12)),
              file.path(ing, "true_pair.fasta"))
  write_fasta(stats::setNames(vapply(1:6, function(i) mutate(locus, 72L), character(1)),
                              paste0("Portulaca_sp", 1:6)),
              file.path(out, "true_pair.fasta"))

  # The trnT-psbD case: same marker name, unrelated sequence, and the outgroup twice as long.
  write_fasta(stats::setNames(vapply(1:12, function(i) mutate(draw(600L), 30L), character(1)),
                              paste0("Opuntia_sp", 1:12)),
              file.path(ing, "false_pair.fasta"))
  write_fasta(stats::setNames(vapply(1:6, function(i) mutate(draw(1340L), 30L), character(1)),
                              paste0("Portulaca_sp", 1:6)),
              file.path(out, "false_pair.fasta"))

  # Same region, deposited in the opposite orientation. Two independent phylotaR runs do not
  # agree on strand, and a check that missed this would condemn good loci.
  write_fasta(stats::setNames(vapply(1:6, function(i) revcomp(mutate(locus, 72L)), character(1)),
                              paste0("Portulaca_sp", 1:6)),
              file.path(out, "revcomp_pair.fasta"))

  true_share  <- PhyloCactus:::.marker_homology_share(
    file.path(ing, "true_pair.fasta"), file.path(out, "true_pair.fasta"), k = 20L)$share
  false_share <- PhyloCactus:::.marker_homology_share(
    file.path(ing, "false_pair.fasta"), file.path(out, "false_pair.fasta"), k = 20L)$share
  rc_share    <- PhyloCactus:::.marker_homology_share(
    file.path(ing, "true_pair.fasta"), file.path(out, "revcomp_pair.fasta"), k = 20L)$share

  expect_gt(true_share, 0.20)
  expect_lt(false_share, 0.05)
  # Orientation must not matter. Only the outgroup side is reverse-complemented, so this is the
  # asymmetry most likely to break.
  expect_gt(rc_share, 0.20)

  # Unrelated sequences share nothing at k = 20 by construction; the pair has to fail at k = 10
  # too, or a divergent real homologue would be condemned alongside it.
  expect_lt(PhyloCactus:::.marker_homology_share(
    file.path(ing, "false_pair.fasta"), file.path(out, "false_pair.fasta"), k = 10L)$share, 0.15)
})

test_that("exclude_markers leaves the locus out of the matrix and the file on disk", {
  skip_if_not_installed("Biostrings")

  tmp <- withr::local_tempdir()
  ind <- file.path(tmp, "aligned_markers"); dir.create(ind)

  # `nuclear_like` is the phyC case: present in the ingroup, absent from the outgroup. It has to
  # be excludable from the dating matrix without being deleted, because the same alignment is the
  # barcoding locus.
  taxa <- c("Opuntia_a", "Opuntia_b", "Portulaca_a", "Portulaca_b")
  set.seed(3L)
  aln <- function(n) vapply(seq_len(n), function(i)
    paste(sample(c("A", "C", "G", "T"), 300L, replace = TRUE), collapse = ""), character(1))
  for (m in c("shared", "nuclear_like")) {
    write_fasta(stats::setNames(aln(4L), taxa), file.path(ind, paste0(m, ".fasta")))
  }

  out <- file.path(tmp, "concat")
  res <- suppressMessages(suppressWarnings(run_concatenation_pipeline(
    input_dir = ind, output_dir = out, exclude_markers = "nuclear_like"
  )))

  ranges <- utils::read.delim(file.path(out, "final_tables", "TABLE_marker_ranges.tsv"))
  expect_equal(sort(ranges$marker), "shared")
  expect_false("nuclear_like" %in% ranges$marker)

  # The exclusion is recorded, not just performed.
  dropped <- utils::read.csv(file.path(out, "logs_and_qc", "TABLE_markers_excluded.csv"))
  expect_equal(dropped$marker, "nuclear_like")

  # And the alignment survives for the analyses that do want it.
  expect_true(file.exists(file.path(ind, "nuclear_like.fasta")))
})

test_that("a marker named in exclude_markers but absent is a warning, not a silent no-op", {
  skip_if_not_installed("Biostrings")

  tmp <- withr::local_tempdir()
  ind <- file.path(tmp, "aligned_markers"); dir.create(ind)
  set.seed(4L)
  write_fasta(stats::setNames(
    vapply(seq_len(4L), function(i)
      paste(sample(c("A", "C", "G", "T"), 300L, replace = TRUE), collapse = ""), character(1)),
    c("Opuntia_a", "Opuntia_b", "Portulaca_a", "Portulaca_b")),
    file.path(ind, "shared.fasta"))

  # A typo in the call site would otherwise leave the locus in the matrix while the script reads
  # as though it had been removed. The pipeline emits other warnings on a fixture this small, so
  # collect them all rather than matching whichever arrives first.
  seen <- character(0)
  withCallingHandlers(
    suppressMessages(run_concatenation_pipeline(
      input_dir = ind, output_dir = file.path(tmp, "concat"), exclude_markers = "phyc"
    )),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("absent from", seen, fixed = TRUE)),
              info = paste("warnings raised:", paste(seen, collapse = " | ")))

  # The real marker survived the typo.
  ranges <- utils::read.delim(file.path(tmp, "concat", "final_tables", "TABLE_marker_ranges.tsv"))
  expect_equal(ranges$marker, "shared")
})
