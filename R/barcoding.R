#' DNA Barcoding Evaluation (Section Under Construction)
#'
#' Placeholder for the molecular diagnostic section of `PhyloCactus`. The analytical
#' strategy for DNA barcoding in Cactaceae has not been settled, so no implementation is
#' shipped: the previous routines were withdrawn in v0.4.2 rather than left in place as
#' working code whose results could be mistaken for validated output.
#'
#' The scope of the section, the decisions already taken, and the questions that remain
#' open are set out in the manifest, `vignette("tutorial-5-cactus-phylogeny-barcoding")`.
#' Single-locus gene tree inference is unaffected and remains available through
#' [infer_gene_trees()], which is documented with the phylogenetic inference module.
#'
#' @name barcoding
#' @seealso [evaluate_dna_barcoding()], [infer_gene_trees()]
NULL

#' Evaluate DNA Barcoding Resolution and Taxonomic Classification Power
#'
#' Not implemented. This function is a documented placeholder that signals an error when
#' called. It is retained so that the section keeps a stable entry point in the reference
#' index and in the tutorial series while its analytical strategy is being defined.
#'
#' The implementation withdrawn in v0.4.2 computed intra- and inter-specific pairwise
#' distance distributions, derived the DNA barcode gap per locus, and trained a
#' `DECIPHER::LearnTaxa` classifier against the botanical backbone of Korotkova
#' *et al*. (2021). It was withdrawn for two reasons. First, its accuracy estimates were
#' produced by evaluating the trained classifier on the same sequence set used for
#' training, which measures memorisation rather than diagnostic performance and therefore
#' reported values that could not be interpreted. Second, the barcode gap thresholds it
#' applied were fixed rather than derived from the empirical distributions of the family,
#' a decision that cannot be made before the phylogenetic backbone and its visualisation
#' are complete.
#'
#' The full rationale, the design decisions already taken, and the criteria that must be
#' met before the section is reopened are documented in
#' `vignette("tutorial-5-cactus-phylogeny-barcoding")`. The withdrawn implementation
#' remains recoverable from the repository history.
#'
#' @param fasta_dir Character. Directory containing curated ingroup FASTA alignments.
#' @param checklist_path Character. Path to the accepted botanical checklist (Korotkova *et al*. 2021).
#' @param output_dir Character. Destination directory for barcoding summary tables.
#' @param distance_models Character vector. Substitution models for pairwise distance estimation.
#' @param train_decipher Logical. Train a taxonomic classifier?
#' @param confidence_threshold Numeric. Bootstrap confidence threshold for classification.
#' @param max_subsample Integer. Maximum number of sequences per species used in distance calculations.
#' @return None. The function signals an error.
#' @references
#' Korotkova, N., Aquino, D., Arias, S., Eggli, U., Franck, A., Gómez-Hinostrosa, C., Guerrero, P. C.,
#' Hernández, H. M., Kohlbecker, A., Köhler, M., Luna, R., Machado, M., Merclinger, M., Nyffeler, R.,
#' Salvador-Montiel, S., Sánchez, D., Schlumpberger, B. O., & Berendsohn, W. G. (2021). Cactaceae at
#' Caryophyllales.org - a dynamic online species-level taxonomic backbone for the family.
#' *Willdenowia*, 51(2), 251-270. \doi{10.3372/wi.51.51208}
#' @seealso [infer_gene_trees()] for single-locus gene tree inference, which remains available.
#' @examples
#' \dontrun{
#' # Signals an error and points to the manifest:
#' evaluate_dna_barcoding("4_Cleaned/cleaned_markers_ingroup", "checklist.csv", "out")
#' }
#' @export
evaluate_dna_barcoding <- function(
  fasta_dir,
  checklist_path,
  output_dir,
  distance_models = c("raw", "K80"),
  train_decipher = TRUE,
  confidence_threshold = 50,
  max_subsample = 20L
) {
  stop(
    "evaluate_dna_barcoding() is not implemented in this version of PhyloCactus.\n",
    "The DNA barcoding section is under construction: its analytical strategy has not ",
    "been defined, so no implementation is shipped rather than one whose output could be ",
    "mistaken for a validated result.\n",
    "Scope, decisions taken and open questions: ",
    "vignette(\"tutorial-5-cactus-phylogeny-barcoding\").\n",
    "Single-locus gene tree inference is unaffected: see infer_gene_trees().",
    call. = FALSE
  )
}
