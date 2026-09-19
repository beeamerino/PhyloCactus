#' DNA Barcoding Evaluation (Section Under Construction)
#'
#' Placeholder for the molecular diagnostic section of `PhyloCactus`. The analytical
#' strategy for DNA barcoding in Cactaceae has not been defined, so the section contains no
#' implementation.
#'
#' The scope of the section, the decisions already taken, and the questions that remain
#' open are set out in the manifest, `vignette("tutorial-5-cactus-phylogeny-barcoding")`.
#'
#' @name barcoding
#' @seealso [evaluate_dna_barcoding()]
NULL

#' Evaluate DNA Barcoding Resolution and Taxonomic Classification Power
#'
#' Not implemented. This function is a documented placeholder that signals an error when
#' called. It keeps a stable entry point for the section in the reference index and in the
#' tutorial series while the analytical strategy is being defined. The scope of the section
#' and its open questions are set out in `vignette("tutorial-5-cactus-phylogeny-barcoding")`.
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
#' Hernández, H. M., Kohlbecker, A., Köhler, M., Luther, K., Majure, L. C., Müller, A., Metzing, D.,
#' Nyffeler, R., Sánchez, D., Schlumpberger, B. O., & Berendsohn, W. G. (2021). Cactaceae at
#' Caryophyllales.org - a dynamic online species-level taxonomic backbone for the family.
#' *Willdenowia*, 51(2), 251-270. \doi{10.3372/wi.51.51208}
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
    "The DNA barcoding section is under construction and its analytical strategy has not ",
    "been defined.\n",
    "Scope and open questions: vignette(\"tutorial-5-cactus-phylogeny-barcoding\").",
    call. = FALSE
  )
}
