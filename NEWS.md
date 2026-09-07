# PhyloCactus 0.4.4

* Updated outgroup sampling in data preparation and tree rooting to include Talinaceae (*Talinum*, *Talinella*) alongside Portulacaceae and Anacampserotaceae.
* Standardized primary clade support metric to Felsenstein Bootstrap Proportions (FBP) across documentation, functions, and collapsing workflows (`integrate_publication_tree()`), retaining Transfer Bootstrap Expectation (TBE) as a complementary diagnostic metric.
* Added support for identifying and distinguishing topologically constrained nodes in tree visualization (`classify_constrained_nodes()`).
* Replaced shell wrapper execution of treePL with direct system calls (`run_treePL_direct()`) and native R cross-validation optimization (`run_treePL_cv()`).
* Retired legacy `run_treePL()` wrapper function.
* Improved documentation, metadata manifests, and error handling across all analytical modules.
