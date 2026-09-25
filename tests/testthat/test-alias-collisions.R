# Regression tests for defect C3 of the 2026-09-10 audit.
#
# `marker_aliases = c(trnL = "trnL_trnF")` pointed two outgroup files at one marker_key. Both were
# read, both were kept, and the exported marker carried Portulaca_oleracea and P. quadrifida twice.
# Matrix assembly then resolved the duplication by assigning rows by name, which takes the first
# match: the supermatrix kept the 274-site copies and discarded the 463 and 314-site ones, in the
# one locus carrying the exemption put there so as not to lose those rooting terminals. Nothing in
# the run said so. These tests hold both halves of the fix: the rule that resolves the collision,
# and the guard that refuses to arbitrate one silently further downstream.

make_set <- function(seqs) Biostrings::DNAStringSet(seqs)

test_that(".resolve_alias_collisions() keeps the longest record and reports every competitor", {
  skip_if_not_installed("Biostrings")

  aln <- make_set(stats::setNames(
    c(paste(rep("A", 274), collapse = ""),
      paste(rep("A", 463), collapse = ""),
      paste(rep("C", 300), collapse = "")),
    c("Portulaca_oleracea", "Portulaca_oleracea", "Opuntia_ficus_indica")
  ))
  src <- c("trnL-trnF.fasta", "trnL.fasta", "trnL-trnF.fasta")

  expect_warning(
    res <- .resolve_alias_collisions(aln, src, "trnL_trnF"),
    "Portulaca_oleracea"
  )

  # No header may appear twice in the exported marker.
  expect_equal(anyDuplicated(names(res$aln)), 0L)
  expect_equal(length(res$aln), 2L)

  # The retained copy is the informative one, not the one whose file sorts first.
  kept <- res$aln[names(res$aln) == "Portulaca_oleracea"]
  expect_equal(unname(Biostrings::width(kept)), 463L)

  # Both competitors are reported, with their provenance, and exactly one is marked retained.
  expect_equal(nrow(res$collisions), 2L)
  expect_equal(sort(res$collisions$source_file), c("trnL-trnF.fasta", "trnL.fasta"))
  expect_equal(sum(res$collisions$retained), 1L)
  expect_equal(res$collisions$source_file[res$collisions$retained], "trnL.fasta")
  expect_equal(res$collisions$marker_key, rep("trnL_trnF", 2L))
})

test_that(".resolve_alias_collisions() counts sites and not columns", {
  skip_if_not_installed("Biostrings")

  # The wider record is mostly gaps and carries fewer sites. Length in columns would keep it; the
  # declared rule is non-gap sites, so it does not.
  aln <- make_set(stats::setNames(
    c(paste0(paste(rep("-", 400), collapse = ""), paste(rep("A", 50), collapse = "")),
      paste(rep("A", 200), collapse = "")),
    c("Portulaca_quadrifida", "Portulaca_quadrifida")
  ))
  src <- c("a_wide.fasta", "b_narrow.fasta")

  suppressWarnings(res <- .resolve_alias_collisions(aln, src, "trnL_trnF"))

  expect_equal(length(res$aln), 1L)
  expect_equal(unname(Biostrings::width(res$aln)), 200L)
  expect_equal(res$collisions$n_non_gap[res$collisions$retained], 200L)
})

test_that(".resolve_alias_collisions() is silent and inert when the aliases collide with nothing", {
  skip_if_not_installed("Biostrings")

  aln <- make_set(stats::setNames(
    c(paste(rep("A", 100), collapse = ""), paste(rep("C", 100), collapse = "")),
    c("Opuntia_ficus_indica", "Portulaca_oleracea")
  ))

  expect_silent(res <- .resolve_alias_collisions(aln, c("x.fasta", "y.fasta"), "matK"))
  expect_equal(names(res$aln), names(aln))
  expect_equal(nrow(res$collisions), 0L)
  # The empty table still carries its columns, so an empty export is a statement and not a shape.
  expect_equal(names(res$collisions),
               c("marker_key", "species", "fasta_name", "source_file", "n_non_gap", "retained"))
})

test_that("the resolution does not depend on the order the files were listed in", {
  skip_if_not_installed("Biostrings")

  long <- paste(rep("A", 463), collapse = "")
  short <- paste(rep("A", 274), collapse = "")

  a <- make_set(stats::setNames(c(short, long), rep("Portulaca_oleracea", 2)))
  b <- make_set(stats::setNames(c(long, short), rep("Portulaca_oleracea", 2)))

  suppressWarnings(res_a <- .resolve_alias_collisions(a, c("f1.fasta", "f2.fasta"), "trnL_trnF"))
  suppressWarnings(res_b <- .resolve_alias_collisions(b, c("f2.fasta", "f1.fasta"), "trnL_trnF"))

  expect_equal(as.character(res_a$aln), as.character(res_b$aln))
  expect_equal(res_a$collisions$source_file[res_a$collisions$retained],
               res_b$collisions$source_file[res_b$collisions$retained])
})


# Added on 2026-09-25, closing the carry-over declared on 2026-09-21. extract_species_binomial()
# turned every hyphen into an underscore and then kept the first two fields, so
# "Opuntia_ficus-indica" came back as "Opuntia_ficus". Counted over the files on 2026-09-25: 13 of
# the 1065 species of the library and 11 of the 1024 tips of the tree carry a hyphenated epithet.
#
# The note of 2026-09-21 said those species would be fused with others. That was wrong and is
# corrected here: applying the old rule to the 1065 species gives 1065 distinct names, no collision.
# The damage is a truncated and false name, not a merge.
#
# The function is not exported and reaches the pipeline only as the default of species_fn in
# .resolve_alias_collisions(), whose single real call passes extract_species(). That is why nothing
# downstream changes, and why this is a latent defect rather than a wrong result already published.

test_that("extract_species_binomial keeps a hyphenated epithet whole", {
  expect_equal(extract_species_binomial("Opuntia_ficus-indica|KJ773783.1"), "Opuntia_ficus-indica")
  expect_equal(extract_species_binomial("Astrophytum_caput-medusae|X1.1"), "Astrophytum_caput-medusae")
  expect_equal(extract_species_binomial("Pachycereus pecten-aboriginum"), "Pachycereus_pecten-aboriginum")
})

test_that("extract_species_binomial still does everything else it did", {
  # Two fields, with a space or with an underscore
  expect_equal(extract_species_binomial("Opuntia robusta|AY1.1"), "Opuntia_robusta")
  expect_equal(extract_species_binomial("Opuntia_robusta|AY1.1"), "Opuntia_robusta")
  # Anything past the epithet is dropped, which is what makes it a binomial
  expect_equal(extract_species_binomial("Opuntia_robusta_var_algo|AY1.1"), "Opuntia_robusta")
  # A single field comes back as it is
  expect_equal(extract_species_binomial("Opuntia|AY1.1"), "Opuntia")
  # Empty and missing
  expect_true(is.na(extract_species_binomial(NA_character_)))
  expect_true(is.na(extract_species_binomial("")))
  expect_true(is.na(extract_species_binomial("|AY1.1")))
  # Vectorised, without names
  expect_equal(extract_species_binomial(c("Opuntia_ficus-indica|A", "Opuntia robusta|B")),
               c("Opuntia_ficus-indica", "Opuntia_robusta"))
})

test_that("the rule separates the thirteen hyphenated species instead of truncating them", {
  trece <- c("Astrophytum_caput-medusae", "Cephalocereus_columna-trajani", "Cereus_pierre-braunianus",
             "Coryphantha_maiz-tablasensis", "Eulychnia_saint-pieana", "Opuntia_ficus-indica",
             "Opuntia_santa-rita", "Pachycereus_pecten-aboriginum", "Peniocereus_lazaro-cardenasii",
             "Pereskia_diaz-romeroana", "Rhipsalis_campos-portoana", "Rhipsalis_neves-armondii",
             "Rhipsalis_pacheco-leonis")
  salida <- extract_species_binomial(paste0(trece, "|sid.1"))
  expect_equal(salida, trece)
  expect_equal(length(unique(salida)), 13L)
})
