# Tutorial 5: DNA Barcoding, Section Manifest

## Abstract

This section of `PhyloCactus` is under construction and ships no
implementation. The document that follows is a manifest, not a tutorial:
it states what the molecular diagnostic section is intended to do, which
design decisions have already been taken, which remain open, and what
must be true before the section is reopened for implementation. Calling
[`evaluate_dna_barcoding()`](https://beeamerino.github.io/PhyloCactus/reference/evaluate_dna_barcoding.md)
signals an error and points here. Single-locus gene tree inference is
unaffected and remains available through
[`infer_gene_trees()`](https://beeamerino.github.io/PhyloCactus/reference/infer_gene_trees.md),
documented with the phylogenetic inference module.

## Why the Section Is Empty Rather Than Provisional

An earlier implementation existed and was withdrawn in v0.4.2. It
computed intra- and inter-specific pairwise distance distributions,
derived a DNA barcode gap per locus, and trained a
[`DECIPHER::LearnTaxa`](https://rdrr.io/pkg/DECIPHER/man/LearnTaxa.html)
classifier against the botanical backbone of Korotkova *et al*. (2021).
Two defects made it unsafe to keep in place.

The classifier was evaluated on the same sequence set used to train it.
That is resubstitution, not validation: it measures how well the model
memorised its training data, and the resulting per-species and per-genus
accuracies were inflated by construction. A user reading those figures
would have had no way to know they did not describe diagnostic
performance.

The barcode gap thresholds were fixed rather than derived from the
empirical distance distributions of the family. In a clade with recent
rapid radiations and low plastid divergence, a fixed threshold encodes
an assumption about divergence that Cactaceae does not satisfy uniformly
across genera.

Neither defect is difficult to correct in isolation. Both are downstream
of a strategic question that is still open, so correcting them now would
produce a working function resting on an undecided design. Shipping
nothing is preferable to shipping output that could be mistaken for a
validated result.

## Scope of the Section

The section is intended to answer one question: given a curated
multilocus reference library for Cactaceae, how reliably can an
individual locus, or a combination of loci, assign an unknown sample to
genus and to species, and where does that resolution break down.

The intended components are the following.

1.  **Empirical barcode gap per locus.** Intra-specific and
    inter-specific pairwise distance distributions computed per marker,
    with the gap reported as an observed quantity rather than tested
    against a fixed cutoff.
2.  **Probabilistic classification with honest validation.** A taxonomic
    classifier trained against the accepted checklist and evaluated
    under stratified cross-validation, with species represented by a
    single accession declared separately because they admit no
    leave-one-out estimate.
3.  **Monophyly in single-locus gene trees.** Recovery of specific and
    generic monophyly marker by marker, interpreted through biological
    processes (incomplete lineage sorting, hybridisation, plastid
    capture) rather than as analytical noise, and reporting paraphyletic
    assemblages as evolutionary grades.

## Decisions Already Taken

These are settled and constrain any future implementation.

- **Taxonomic backbone.** Nomenclatural reconciliation is anchored in
  the dynamic Cactaceae checklist of Caryophyllales.org (Korotkova *et
  al*. 2021), the same backbone used throughout the assembly modules. No
  parallel taxonomy will be introduced.
- **Input.** The section consumes the curated ingroup alignments in
  `4_Cleaned/cleaned_markers_ingroup`, produced by Module 4, and the
  phylogenetic products of Module 8 onwards in `7_Phylogenetics/`. It
  does not re-mine GenBank.
- **Gene tree inference stays out.**
  [`infer_gene_trees()`](https://beeamerino.github.io/PhyloCactus/reference/infer_gene_trees.md)
  was relocated to the phylogenetic inference module because its
  products feed the gene tree versus species tree discordance analyses
  of Module 13, which do not depend on the barcoding strategy.
- **No fixed distance thresholds.** Any decision rule must be derived
  from the observed distributions of the family, and its derivation must
  be reported alongside the result.
- **Validation is not optional.** No classifier accuracy will be
  reported without cross-validation stratified by species, and the
  treatment of singleton species will be declared explicitly.

## Open Questions

These must be resolved before implementation begins.

1.  **Which loci, and combined how.** Whether resolution is assessed per
    locus, on concatenated subsets, or on both, and whether the standard
    barcode combination is evaluated as a unit for comparability with
    published surveys.
2.  **The unit of assignment.** Whether the target is species-level
    assignment, generic assignment with an explicit statement of
    specific ambiguity, or a graded output. The biological ceiling on
    monophyly in plants makes universal species assignment unattainable,
    so the reporting format has to accommodate justified non-assignment.
3.  **How ambiguity is reported.** Whether a sample that cannot be
    separated from a set of candidates is reported as a candidate set,
    as a genus-level assignment, or as a failure, and how that interacts
    with the conservation applications of the package.
4.  **Relationship to the genome skimming extension.** Whether this
    section remains a Sanger-era diagnostic evaluation, or is reframed
    as the reference-library component of a super-barcoding workflow.
    The two designs imply different outputs and different validation
    regimes.
5.  **Handling of paralogy and organelle capture.** Whether loci showing
    cyto-nuclear conflict are excluded from diagnostic use, flagged, or
    reported as evidence of reticulation in their own right.

## Reopening Criteria

The section will be reopened for implementation when all of the
following hold.

- Modules 7 to 10 are complete: a dated, supported phylogenetic backbone
  exists, so that monophyly assessments per locus can be read against a
  reference topology rather than in isolation.
- Modules 11 and 12 are complete: the visualisation and metadata
  integration layer exists, so that diagnostic results can be presented
  against the conservation and distribution attributes they are meant to
  inform.
- The five open questions above are answered and recorded.
- A validation plan is written before any code, specifying the
  cross-validation scheme, the treatment of singleton species, and the
  negative controls.

## Conclusion

Until those conditions are met, the correct state of this section is
empty. The entry point, the documentation and the tutorial slot are
retained so that the structure of the package reflects the intended
workflow, and so that a user who reaches this stage is told what is
missing rather than handed an unvalidated result.

[Continue to Tutorial 6: Function Reference &
Dictionary](https://beeamerino.github.io/PhyloCactus/articles/tutorial-6-cactus-phylogeny-functions.html)

## References

Korotkova, N., Aquino, D., Arias, S., Eggli, U., Franck, A.,
Gómez-Hinostrosa, C., Guerrero, P. C., Hernández, H. M., Kohlbecker, A.,
Köhler, M., Luther, K., Majure, L. C., Müller, A., Metzing, D.,
Nyffeler, R., Sánchez, D., Schlumpberger, B. O., & Berendsohn, W. G.
(2021). Cactaceae at Caryophyllales.org, a dynamic online species-level
taxonomic backbone for the family. *Willdenowia*, 51(2), 251-270.
<https://doi.org/10.3372/wi.51.51208>

Wright, E. S. (2016). Using DECIPHER v2.0 to analyze big biological
sequence data in R. *The R Journal*, 8(1), 352-359.
<https://doi.org/10.32614/RJ-2016-025>
