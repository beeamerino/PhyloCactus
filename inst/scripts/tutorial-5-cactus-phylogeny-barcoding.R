# -------------------------------------------------------------
# PhyloCactus: Tutorial 5. DNA Barcoding (section under construction)
# -------------------------------------------------------------
# This section contains no implementation, so no runnable script is provided.
# Scope, decisions taken, open questions and conditions for implementation:
#
#   vignette("tutorial-5-cactus-phylogeny-barcoding", package = "PhyloCactus")
#
# evaluate_dna_barcoding() signals an error that points to that page.
#
# Single-locus gene trees can still be estimated with infer_gene_trees(). Its
# gene trees can be used as input to ASTRAL-III (Module 13):
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
