# Changelog

## PhyloCactus 0.4.4

- Updated outgroup sampling in data preparation and tree rooting to
  include Talinaceae (*Talinum*, *Talinella*) alongside Portulacaceae
  and Anacampserotaceae.
- Standardized primary clade support metric to Felsenstein Bootstrap
  Proportions (FBP) across documentation, functions, and collapsing
  workflows
  ([`integrate_publication_tree()`](https://beeamerino.github.io/PhyloCactus/reference/integrate_publication_tree.md)),
  retaining Transfer Bootstrap Expectation (TBE) as a complementary
  diagnostic metric.
- Added support for identifying and distinguishing topologically
  constrained nodes in tree visualization
  ([`classify_constrained_nodes()`](https://beeamerino.github.io/PhyloCactus/reference/classify_constrained_nodes.md)).
- Replaced shell wrapper execution of treePL with direct system calls
  ([`run_treePL_direct()`](https://beeamerino.github.io/PhyloCactus/reference/run_treePL_direct.md))
  and native R cross-validation optimization
  ([`run_treePL_cv()`](https://beeamerino.github.io/PhyloCactus/reference/run_treePL_cv.md)).
- Retired legacy `run_treePL()` wrapper function.
- Improved documentation, metadata manifests, and error handling across
  all analytical modules.
