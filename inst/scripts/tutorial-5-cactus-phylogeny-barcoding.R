# -------------------------------------------------------------
# PhyloCactus: Tutorial 5 - DNA Barcoding (section under construction)
# -------------------------------------------------------------
# This section ships no implementation. Its analytical strategy has not been
# defined, so no runnable script is provided rather than one whose output could
# be mistaken for a validated result.
#
# Scope of the section, decisions already taken, open questions and the criteria
# for reopening it are documented in the manifest:
#
#   vignette("tutorial-5-cactus-phylogeny-barcoding", package = "PhyloCactus")
#
# Calling evaluate_dna_barcoding() signals an error pointing to that manifest.
#
# Single-locus gene tree inference is NOT part of the withdrawal. It was
# relocated to the phylogenetic inference module because its products feed the
# gene tree / species tree discordance analyses of Module 13 (ASTRAL-III), which
# do not depend on the barcoding strategy. It remains available:
#
#   library(PhyloCactus)
#
#   gene_trees_manifest <- infer_gene_trees(
#     fasta_dir        = "4_Cleaned/cleaned_markers_ingroup",
#     output_dir       = "4b_Gene_Trees",
#     include_outgroup = FALSE,   # TRUE when preparing unrooted trees for ASTRAL-III
#     method           = "auto",
#     model            = "GTR+G",
#     threads          = 2L
#   )
# -------------------------------------------------------------
