# Tutorial 5: DNA Barcoding, Section Manifest

## Abstract

This section of `PhyloCactus` is under construction and contains no
implementation. This document sets out what the molecular diagnostic
section is intended to do, the design decisions already taken, the
questions that remain open, and the conditions for implementing it.
[`evaluate_dna_barcoding()`](https://beeamerino.github.io/PhyloCactus/reference/evaluate_dna_barcoding.md)
signals an error that points to this page. Single-locus gene trees can
still be estimated with
[`infer_gene_trees()`](https://beeamerino.github.io/PhyloCactus/reference/infer_gene_trees.md).

## Scope

The section is intended to answer one question: given a curated
multilocus reference library for Cactaceae, how reliably can a single
locus, or a combination of loci, assign an unknown sample to genus and
to species, and in which groups does that assignment fail?

The planned components are:

1.  **Empirical barcode gap per locus.** Intraspecific and interspecific
    pairwise distance distributions for each marker, with the gap
    reported as an observed quantity and not tested against a fixed
    cutoff.
2.  **Probabilistic classification with cross-validation.** A taxonomic
    classifier trained on sequences of accepted species and evaluated by
    cross-validation stratified by species. Species represented by a
    single accession cannot be evaluated by leaving one sequence out,
    and are reported separately.
3.  **Monophyly in single-locus gene trees.** Recovery of species and
    genus monophyly for each marker. Non-monophyly is interpreted in
    terms of incomplete lineage sorting, hybridization and plastid
    capture, and paraphyletic assemblages are reported as grades.

## Decisions Taken

- **Taxonomic backbone.** Names are reconciled against the Cactaceae
  checklist of Caryophyllales.org (Korotkova *et al*. 2021), the same
  backbone used in the assembly modules.
- **Input.** The section uses the curated ingroup alignments in
  `4_Cleaned/cleaned_markers_ingroup` (Module 4) and the phylogenetic
  products of Modules 7 to 10. It does not query GenBank again.
- **Gene trees.**
  [`infer_gene_trees()`](https://beeamerino.github.io/PhyloCactus/reference/infer_gene_trees.md)
  belongs to the inference module. Its gene trees can be used as input
  to `ASTRAL-III` (Module 13) and do not depend on the design of this
  section.
- **No fixed distance thresholds.** Any decision rule must be derived
  from the observed distance distributions of the family, and its
  derivation must be reported with the result. In a family with recent
  radiations and low plastid divergence, a fixed threshold assumes a
  level of divergence that is not uniform across genera.
- **Validation.** Classification accuracy is reported only from
  cross-validation stratified by species. Accuracy measured on the
  training sequences reflects how well the model fits those sequences,
  not its performance on new samples.

## Open Questions

1.  **Loci and combinations.** Whether resolution is assessed per locus,
    on concatenated subsets, or both, and whether the standard plant
    barcode combination is evaluated as a unit for comparison with
    published surveys.
2.  **Unit of assignment.** Whether the target is assignment to species,
    assignment to genus with an explicit statement of ambiguity at
    species level, or a graded output. Many plant species are not
    monophyletic in single-locus trees, so the output must allow for
    justified non-assignment.
3.  **Reporting ambiguity.** Whether a sample that cannot be separated
    from a set of candidates is reported as a candidate set, as an
    assignment to genus, or as a failure, and how this affects the
    conservation applications of the package.
4.  **Relationship to genome skimming.** Whether the section remains an
    evaluation of Sanger markers or becomes the reference-library
    component of a plastome-based (super-barcoding) workflow. The two
    designs require different outputs and different validation.
5.  **Paralogy and organelle capture.** Whether loci with cytonuclear
    conflict are excluded from diagnostic use, flagged, or reported as
    evidence of reticulation.

## Conditions for Implementation

Two conditions are met in version 0.4.5: a dated and supported phylogeny
is available (Modules 7 to 10), so monophyly per locus can be assessed
against a reference topology, and the metadata and visualization modules
are available (Modules 11 and 12), so diagnostic results can be related
to conservation status and distribution. Two conditions remain:

- the five open questions above are answered and recorded;
- a validation plan is written before any code, specifying the
  cross-validation scheme, the treatment of species with a single
  accession, and the negative controls.

[Continue to Tutorial 6: Core Pipeline Functions and Methodological
Dictionary](https://beeamerino.github.io/PhyloCactus/articles/tutorial-6-cactus-phylogeny-functions.html)

## References

- Korotkova *et al*. 2021. Cactaceae at Caryophyllales.org - A dynamic
  online species-level taxonomic backbone for the family. *Willdenowia*,
  51(2), 251–270. <https://doi.org/10.3372/wi.51.51208>
