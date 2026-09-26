# Monophyly per locus of the molecular diagnostic branch (11_barcoding/12_monophyly/).
# Phase 6E of PhyloCactus 0.5.0, decisions N1 to N10 of BMM (26-09); section 8, point 5, and
# section 6.5 of the validation plan. One gene tree per locus, inferred on the library alignment of
# step 4 with one tip per library sequence; monophyly read on its splits, grades on a root taken
# from the reference ML tree.

#' SLURM header shared by the job scripts of the branch
#' @noRd
.bc_slurm_header <- function(name, log, mail_type, cluster_partition, cluster_mem, cluster_time,
                             cluster_queue, cluster_mail_user, load_module, cpus = 1L, array = NULL) {
  h <- c("#!/bin/bash", paste0("#SBATCH -J ", name), paste0("#SBATCH -p ", cluster_partition),
         "#SBATCH --nodes=1", "#SBATCH --ntasks-per-node=1", paste0("#SBATCH --cpus-per-task=", cpus),
         paste0("#SBATCH --mem=", cluster_mem), paste0("#SBATCH -t ", cluster_time))
  if (!is.null(array)) h <- c(h, paste0("#SBATCH --array=", array))
  h <- c(h, paste0("#SBATCH -o ", log, ".out"), paste0("#SBATCH -e ", log, ".err"))
  if (!is.null(cluster_queue) && nzchar(cluster_queue)) h <- c(h, paste0("#SBATCH -q ", cluster_queue))
  if (!is.null(cluster_mail_user) && nzchar(trimws(cluster_mail_user))) {
    h <- c(h, paste0("#SBATCH --mail-type=", mail_type), paste0("#SBATCH --mail-user=", trimws(cluster_mail_user)))
  }
  h <- c(h, "#", "")
  if (!is.null(load_module) && length(load_module) > 0L && nzchar(trimws(paste(load_module, collapse = " ")))) {
    h <- c(h, paste0("module load ", paste(load_module, collapse = " ")), "")
  }
  h
}

#' Species of a tip labelled `species|sid`, in the comparison key of the branch
#' @noRd
.bc_tip_species <- function(tips) .bc_checklist_key(sub("\\|.*$", "", tips))

#' The backbone of the constraint of Module 9, as the nesting of its clades
#'
#' The same clades and the same nesting as `build_constraint_scaffold()` for Cactaceae, which the
#' library of the branch is limited to. A clade absent from a locus is dropped, and a node left with
#' one child is replaced by it.
#' @noRd
.bc_module9_backbone <- function() {
  list("Leuenbergeria",
       list("Pereskia",
            list(list("Cylindropuntieae", list("Tephrocacteae", "Opuntieae")),
                 list("Maihuenia",
                      list("Blossfeldia",
                           list("Cacteae",
                                list("Frailea", "Calymmanthium", "Copiapoa",
                                     list(list("Rhipsalideae", list("Notocacteae", "BCT")), "Core I"))))))))
}

#' Constraint tree of one locus, one tip per library sequence, from the clades of Module 9
#'
#' Each tip joins the clade of its species in `constraints_csv`; tips of species absent from the
#' table stay out of the constraint and are free in the search, which RAxML-NG allows.
#' @noRd
.bc_gene_tree_constraint <- function(tips, constraints_csv) {
  tab <- utils::read.csv(constraints_csv, stringsAsFactors = FALSE)
  tab$key <- .bc_checklist_key(tab$Specie_name)
  tab <- unique(tab[, c("key", "Clade")])
  dup <- unique(tab$key[duplicated(tab$key)])
  if (length(dup) > 0L) {
    stop("Species with more than one clade in ", constraints_csv, ": ", paste(utils::head(dup, 5), collapse = ", "),
         call. = FALSE)
  }
  clade_of <- stats::setNames(tab$Clade, tab$key)[.bc_tip_species(tips)]
  build <- function(node) {
    if (is.character(node)) {
      t <- sort(tips[!is.na(clade_of) & clade_of == node], method = "radix")
      if (length(t) == 0L) return("")
      if (length(t) == 1L) return(t)
      return(paste0("(", paste(t, collapse = ","), ")"))
    }
    parts <- vapply(node, build, character(1))
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0L) return("")
    if (length(parts) == 1L) return(parts)
    paste0("(", paste(parts, collapse = ","), ")")
  }
  s <- build(.bc_module9_backbone())
  if (!nzchar(s)) return(NULL)
  if (!startsWith(s, "(")) s <- paste0("(", s, ")")
  paste0(s, ";")
}

#' Generate the SLURM Job Scripts of the Gene Trees of the Molecular Diagnostic Branch
#'
#' Writes one job per locus that infers the gene tree of its library alignment (step 4) with
#' RAxML-NG: an unconstrained or constrained search from 20 starting trees, bootstrap replicates
#' until the autoMRE test converges (at most 1000), and the bootstrap support (FBP) mapped on the
#' best tree, all in one `--all` call after the alignment is parsed to RBA. Phase 6E, decisions N1,
#' N2 and N10 of BMM. The job follows the conventions of [generate_ml_search_script()]: the MPI
#' build and its modules, workers over the starting trees, `--force perf_threads` and
#' `--extra thread-nopin`, thread binding off.
#'
#' @details
#' Each locus gets a fixed seed derived from its name and `seed`, the same on any machine. With
#' `constraint = "module9"` each tip joins the clade of its species in the constraint table of the
#' phylogeny, with the nesting of `build_constraint_scaffold()`; tips of species absent from the
#' table stay free. The paths are written as given and must be those of the cluster.
#'
#' @param library_dir Character. Directory of step 4, with `LIB_<locus>.fasta`.
#' @param trees_dir Character. Where RAxML-NG writes, as `<locus>.raxml.*`.
#' @param loci Character vector of loci, or `NULL` for all.
#' @param constraint Character. `"none"` or `"module9"`.
#' @param constraints_csv Character. Constraint table of the phylogeny. Defaults to the one of the
#'   package.
#' @param model,starting_trees,bs_trees Character. RAxML-NG options. Default `"GTR+G4"`,
#'   `"pars{10},rand{10}"` and `"autoMRE{1000}"`.
#' @param bs_cutoff Numeric. Cutoff of the autoMRE convergence test (`--bs-cutoff`), written in the
#'   job rather than left to the RAxML-NG default. Defaults to 0.03, the value of the bootstrap
#'   convergence test of the phylogeny.
#' @param seed Integer. Seed of the run, from which each locus derives its own.
#' @param threads Integer. Cores per job (`--cpus-per-task` and `--threads`). Defaults to 32.
#' @param workers,min_threads_per_worker Integer. RAxML-NG workers over the starting trees; `NULL`
#'   plans them as [generate_ml_search_script()] does.
#' @param job_dir Character. Where the scripts are written. Defaults to `job/` inside `trees_dir`.
#' @param cluster_job_name,cluster_partition,cluster_mem,cluster_time,cluster_queue,cluster_mail_user,load_module
#'   As in [generate_barcoding_job_scripts()]; with an address every SLURM mail is sent
#'   (`--mail-type=ALL`), as in the jobs of the phylogeny.
#' @param raxml_exec Character. The RAxML-NG command. Defaults to the MPI build used by the
#'   phylogeny on Leftraru, `raxml-ng-mpi`, with its modules in `load_module`.
#' @return Invisibly, the paths of the job scripts.
#' @seealso [assess_barcoding_monophyly()], [generate_barcoding_job_scripts()].
#' @examples
#' \dontrun{
#' generate_barcoding_gene_tree_scripts(
#'   library_dir = "/home/user/PhyloCactus_Tutorial/11_barcoding/4_library",
#'   trees_dir = "/home/user/PhyloCactus_Tutorial/11_barcoding/12_monophyly/trees",
#'   loci = "trnL-trnF"
#' )
#' }
#' @export
generate_barcoding_gene_tree_scripts <- function(library_dir, trees_dir, loci = NULL,
                                                 constraint = c("none", "module9"),
                                                 constraints_csv = system.file("extdata", "cactus_constraints.csv",
                                                                               package = "PhyloCactus"),
                                                 model = "GTR+G4", starting_trees = "pars{10},rand{10}",
                                                 bs_trees = "autoMRE{1000}", bs_cutoff = 0.03,
                                                 seed = 1L, threads = 32L,
                                                 workers = NULL, min_threads_per_worker = 4L,
                                                 job_dir = file.path(trees_dir, "job"),
                                                 cluster_job_name = "cactus_genetree",
                                                 cluster_partition = "main", cluster_mem = "16G",
                                                 cluster_time = "48:00:00", cluster_queue = NULL,
                                                 cluster_mail_user = Sys.getenv("MY_EMAIL", ""),
                                                 load_module = c("gcc/14.2.0-nlhpc", "openmpi/5.0.3-o",
                                                                 "raxml-ng/1.1.0-mpi-zen4-n"),
                                                 raxml_exec = "raxml-ng-mpi") {
  constraint <- match.arg(constraint)
  files <- list.files(library_dir, pattern = "^LIB_.*\\.fasta$", full.names = TRUE)
  if (length(files) == 0L) {
    stop("No LIB_<locus>.fasta in ", library_dir, ". Run finalize_barcoding_library() first.", call. = FALSE)
  }
  all_loci <- sub("^LIB_(.*)\\.fasta$", "\\1", basename(files))
  loci <- if (is.null(loci)) sort(all_loci, method = "radix") else intersect(loci, all_loci)
  if (length(loci) == 0L) stop("None of the loci asked for is in ", library_dir, ".", call. = FALSE)
  if (constraint == "module9" && !file.exists(constraints_csv)) {
    stop("No constraint table at ", constraints_csv, ".", call. = FALSE)
  }
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
  threads <- as.integer(threads)
  n_trees <- .parse_n_start_trees(starting_trees)
  plan <- if (is.null(workers)) .plan_ml_workers(n_trees, threads, min_threads_per_worker) else
    list(workers = max(1L, as.integer(workers)))
  jobs <- character(0)
  for (l in loci) {
    msa <- file.path(library_dir, paste0("LIB_", l, ".fasta"))
    l_seed <- .bc_query_seed(seed, l, "gene_tree", 0L, "<raxml>")
    cons_line <- character(0)
    if (constraint == "module9") {
      tips <- sub("^>", "", grep("^>", readLines(msa), value = TRUE))
      nwk <- .bc_gene_tree_constraint(tips, constraints_csv)
      if (!is.null(nwk)) {
        f_cons <- file.path(job_dir, paste0("constraint_", l, ".tree"))
        writeLines(nwk, f_cons)
        cons_line <- paste0(" --tree-constraint ", shQuote(f_cons))
      }
    }
    job <- c(.bc_slurm_header(paste0(cluster_job_name, "_", l), paste0("slurm_", l, "_%j"), "ALL",
                              cluster_partition, cluster_mem, cluster_time, cluster_queue,
                              cluster_mail_user, load_module, cpus = threads),
             paste0("# Gene tree of ", l, ", Phase 6E. Written by generate_barcoding_gene_tree_scripts() on ",
                    format(Sys.time(), "%Y-%m-%d %H:%M"), "; constraint: ", constraint, "; seed ", l_seed, "."),
             "echo \"=== RAxML-NG version used by this job ===\"",
             paste0(raxml_exec, " --version | head -n 2"),
             "echo \"=========================================\"",
             "",
             "# Thread binding off, as in the ML search of the phylogeny on Leftraru",
             "export OMP_PROC_BIND=false",
             "export OMPI_MCA_hwloc_base_binding_policy=none",
             "",
             paste0("mkdir -p ", shQuote(trees_dir)),
             paste0("RBA_PREFIX=", shQuote(file.path(trees_dir, paste0(l, "_parsed")))),
             "RBA_FILE=\"${RBA_PREFIX}.raxml.rba\"",
             "if [ ! -f \"$RBA_FILE\" ]; then",
             paste0("    ", raxml_exec, " --parse --msa ", shQuote(msa), " --model ", model, " --prefix \"$RBA_PREFIX\""),
             "fi",
             "",
             paste0("echo \"Gene tree of ", l, ": ", n_trees, " starting trees over ", plan$workers,
                    " workers, bootstraps ", bs_trees, ", autoMRE cutoff ", format(bs_cutoff), "\""),
             paste0(raxml_exec, " --all --msa \"$RBA_FILE\" --tree ", shQuote(starting_trees, type = "sh"),
                    " --bs-trees ", shQuote(bs_trees, type = "sh"), " --bs-cutoff ", format(bs_cutoff),
                    " --bs-metric fbp --seed ", l_seed,
                    " --threads ${SLURM_CPUS_PER_TASK:-", threads, "} --workers ", plan$workers,
                    " --force perf_threads --extra thread-nopin",
                    " --prefix ", shQuote(file.path(trees_dir, l)), cons_line))
    f_job <- file.path(job_dir, paste0("job_", l, ".sh"))
    writeLines(job, f_job)
    Sys.chmod(f_job, mode = "0755")
    jobs <- c(jobs, f_job)
  }
  submit <- c("#!/bin/bash", "# Submits one gene-tree job per locus. Written by generate_barcoding_gene_tree_scripts().",
              "set -euo pipefail", paste0("cd ", shQuote(job_dir)),
              paste0("sbatch ", basename(jobs)))
  writeLines(submit, file.path(job_dir, "submit.sh"))
  Sys.chmod(file.path(job_dir, "submit.sh"), mode = "0755")
  message("Gene-tree jobs for ", length(loci), " loci written to ", job_dir,
          ". Submit on the cluster with: bash ", file.path(job_dir, "submit.sh"))
  invisible(jobs)
}

#' The non-trivial splits of a tree, with the support of each
#'
#' One row per internal branch: which tips are on its side, taken as the descendants of its lower
#' node in the tree as written, and its support from the node label (NA if absent). A split and its
#' complement are the same bipartition, so the root of the written tree does not matter.
#' @noRd
.bc_tree_splits <- function(tree) {
  n <- length(tree$tip.label)
  parts <- ape::prop.part(tree)
  sup <- if (is.null(tree$node.label)) rep(NA_real_, length(parts)) else
    suppressWarnings(as.numeric(tree$node.label))
  m <- matrix(FALSE, nrow = length(parts), ncol = n)
  for (k in seq_along(parts)) m[k, parts[[k]]] <- TRUE
  keep <- rowSums(m) >= 2L & rowSums(!m) >= 2L
  list(members = m[keep, , drop = FALSE], support = sup[keep])
}

#' Status of one group in an unrooted tree: monophyletic, rejected or not rejected
#' @noRd
.bc_group_status <- function(g, splits, cutoff) {
  if (sum(g) < 2L) return(list(status = "not_evaluable", support = NA_real_))
  if (sum(!g) <= 1L) return(list(status = "monophyletic", support = NA_real_))
  # The four intersections of each split with the group and its complement, counted at once
  m <- splits$members
  n <- length(g); ng <- sum(g); ns <- rowSums(m)
  sg <- as.vector(m %*% as.numeric(g))
  s_not_g <- ns - sg; not_s_g <- ng - sg; not_s_not_g <- n - ns - not_s_g
  same <- (sg == ng & s_not_g == 0) | (sg == 0 & not_s_g == ng & not_s_not_g == 0)
  if (any(same)) return(list(status = "monophyletic", support = splits$support[which(same)[1]]))
  clash <- sg > 0 & s_not_g > 0 & not_s_g > 0 & not_s_not_g > 0
  strong <- clash & !is.na(splits$support) & splits$support >= cutoff
  list(status = if (any(strong)) "rejected" else "not_rejected", support = NA_real_)
}

#' Grade of a rejected group in the rooted tree: paraphyletic when the other tips under its most
#' recent common ancestor form one clade (one nested lineage), polyphyletic otherwise
#' @noRd
.bc_group_grade <- function(rooted, members) {
  mrca <- ape::getMRCA(rooted, members)
  desc <- ape::extract.clade(rooted, mrca)$tip.label
  others <- setdiff(desc, members)
  if (length(others) == 0L) return(NA_character_)
  if (length(others) == 1L || ape::is.monophyletic(rooted, others, reroot = FALSE)) "paraphyletic" else "polyphyletic"
}

#' The earliest lineage of the rooted reference present among `present` species
#'
#' From the root of the ingroup, the smaller side of each node is the lineage that diverged first;
#' if none of its species is in the locus, the walk goes on along the other side. Returns the
#' species of the lineage, or NULL.
#' @noRd
.bc_root_lineage <- function(ref, ingroup, present) {
  node <- ape::getMRCA(ref, ingroup)
  n <- length(ref$tip.label)
  tips_of <- function(x) if (x <= n) ref$tip.label[x] else ape::extract.clade(ref, x)$tip.label
  queue <- ref$edge[ref$edge[, 1] == node, 2]
  while (length(queue) > 0L) {
    sizes <- vapply(queue, function(x) length(tips_of(x)), integer(1))
    first <- vapply(queue, function(x) sort(tips_of(x), method = "radix")[1], character(1))
    k <- order(sizes, first, method = "radix")[1]
    lin <- tips_of(queue[k])
    if (any(lin %in% present)) return(lin)
    queue <- queue[-k]
    if (length(queue) == 1L) {
      x <- queue
      queue <- if (x <= n) integer(0) else ref$edge[ref$edge[, 1] == x, 2]
      if (x <= n && ref$tip.label[x] %in% present) return(ref$tip.label[x])
    }
  }
  NULL
}

#' Assess the Monophyly of Species and Genera per Locus
#'
#' Step 12 of the branch (Phase 6E). Reads the gene tree of each locus with its bootstrap support
#' (`<locus>.raxml.support` from [generate_barcoding_gene_tree_scripts()]) and writes, per locus,
#' whether each species and each genus is monophyletic, and sets each genus beside the reference ML
#' tree of the package and beside the other loci of its compartment.
#'
#' @details
#' A group with two or more tips is **monophyletic** when its tips are one side of a split of the
#' unrooted gene tree, whatever the root. It is **rejected** when it conflicts with a split of
#' support `support_cutoff` or more, and **not rejected** when it conflicts only with weaker
#' splits. A species with one tip, and a genus with one species in the locus, are **not
#' evaluable**.
#'
#' A rejected group is classified as a **paraphyletic** grade when the other tips under its most
#' recent common ancestor form one clade, and as **polyphyletic** otherwise. Those two need a root:
#' the gene tree is rooted on the tips of the earliest-diverging lineage of the reference present in
#' the locus. If those tips are not a split of the gene tree, the tree is left unrooted, flagged in
#' the summary, and no grade is classified.
#'
#' @param trees_dir Character. Directory with `<locus>.raxml.support`.
#' @param output_dir Character. Where the tables are written.
#' @param reference_tree Character. The reference ML tree, one tip per species. Defaults to the one
#'   of the package.
#' @param reference_outgroup Character vector. Tips of the reference that root it. `NULL` takes the
#'   outgroup species of the species table of the package.
#' @param support_cutoff Numeric. Bootstrap support (FBP, 0 to 100) from which a conflicting split
#'   rejects monophyly. Defaults to 70.
#' @param nuclear_loci Character vector. Loci of the nuclear compartment; the rest are plastid.
#' @return Invisibly, a list with the species, genus, summary and reference tables.
#' @seealso [generate_barcoding_gene_tree_scripts()], [summarise_barcoding_metrics()].
#' @examples
#' \dontrun{
#' assess_barcoding_monophyly()
#' }
#' @export
assess_barcoding_monophyly <- function(trees_dir = file.path("11_barcoding", "12_monophyly", "trees"),
                                       output_dir = file.path("11_barcoding", "12_monophyly"),
                                       reference_tree = system.file("extdata", "phylocactus_ml_tree.tree",
                                                                    package = "PhyloCactus"),
                                       reference_outgroup = NULL, support_cutoff = 70,
                                       nuclear_loci = c("ITS", "phyC", "pepC_like")) {
  .bc_assert_output_dir(output_dir)
  files <- list.files(trees_dir, pattern = "\\.raxml\\.support$", full.names = TRUE)
  if (length(files) == 0L) {
    stop("No <locus>.raxml.support in ", trees_dir, ". Run the jobs of generate_barcoding_gene_tree_scripts() first.",
         call. = FALSE)
  }
  ref <- ape::read.tree(reference_tree)
  if (is.null(reference_outgroup)) {
    sp <- utils::read.csv(system.file("extdata", "phylocactus_table_species.csv", package = "PhyloCactus"),
                          stringsAsFactors = FALSE)
    reference_outgroup <- sp$species[sp$is_outgroup %in% TRUE]
  }
  ref <- root_on_clade(ref, intersect(reference_outgroup, ref$tip.label))
  og_key <- .bc_checklist_key(reference_outgroup)
  ref$tip.label <- .bc_checklist_key(ref$tip.label)
  ingroup_ref <- setdiff(ref$tip.label, og_key)
  ref_genus <- .bc_genus(ref$tip.label)

  sp_rows <- list(); ge_rows <- list(); sm_rows <- list(); rf_rows <- list()
  for (f in sort(files, method = "radix")) {
    l <- sub("\\.raxml\\.support$", "", basename(f))
    tree <- ape::read.tree(f)
    tips <- tree$tip.label
    species <- .bc_tip_species(tips)
    genus <- .bc_genus(species)
    splits <- .bc_tree_splits(tree)

    # Root for the grades
    lin <- .bc_root_lineage(ref, ingroup_ref, unique(species))
    rooted <- NULL; lineage_label <- NA_character_
    if (!is.null(lin)) {
      lt <- tips[species %in% lin]
      st <- .bc_group_status(tips %in% lt, splits, Inf)
      if (length(lt) == 1L || st$status == "monophyletic") {
        rooted <- tryCatch(root_on_clade(tree, lt), error = function(e) NULL)
        if (!is.null(rooted)) lineage_label <- paste(sort(unique(.bc_genus(lin)), method = "radix"), collapse = ";")
      }
    }

    one <- function(groups, key, kind) {
      out <- lapply(groups, function(k) {
        g <- key == k
        n_sp <- length(unique(species[g]))
        st <- if (kind == "genus" && n_sp < 2L) list(status = "not_evaluable", support = NA_real_) else
          .bc_group_status(g, splits, support_cutoff)
        grade <- if (st$status == "rejected" && !is.null(rooted)) .bc_group_grade(rooted, tips[g]) else NA_character_
        data.frame(locus = l, group = k, tips = sum(g), species_in_locus = n_sp, status = st$status,
                   support = st$support, grade = grade, stringsAsFactors = FALSE)
      })
      do.call(rbind, out)
    }
    s_tab <- one(sort(unique(species), method = "radix"), species, "species")
    s_tab$species_in_locus <- NULL
    s_tab$in_reference <- s_tab$group %in% ref$tip.label
    g_tab <- one(sort(unique(genus), method = "radix"), genus, "genus")
    g_tab$compartment <- if (l %in% nuclear_loci) "nuclear" else "plastid"
    sp_rows[[l]] <- s_tab; ge_rows[[l]] <- g_tab

    ref_mono <- vapply(g_tab$group, function(gname) {
      m <- ref$tip.label[ref_genus == gname & ref$tip.label %in% ingroup_ref]
      if (length(m) < 2L) NA else ape::is.monophyletic(ref, m, reroot = FALSE)
    }, logical(1))
    rf_rows[[l]] <- data.frame(genus = g_tab$group, locus = l, compartment = g_tab$compartment,
                               gene_tree_status = g_tab$status, reference_monophyletic = unname(ref_mono),
                               stringsAsFactors = FALSE)

    count <- function(t, what) c(evaluable = sum(t$status != "not_evaluable"),
                                 monophyletic = sum(t$status == "monophyletic"),
                                 rejected = sum(t$status == "rejected"),
                                 not_rejected = sum(t$status == "not_rejected"),
                                 not_evaluable = sum(t$status == "not_evaluable"),
                                 paraphyletic = sum(t$grade %in% "paraphyletic"),
                                 polyphyletic = sum(t$grade %in% "polyphyletic"))
    cs <- count(s_tab); cg <- count(g_tab)
    sm <- data.frame(locus = l, tips = length(tips), rooted = !is.null(rooted), root_lineage = lineage_label,
                     support_cutoff = support_cutoff, stringsAsFactors = FALSE)
    for (nm in names(cs)) sm[[paste0("species_", nm)]] <- as.integer(cs[[nm]])
    sm$species_prop_monophyletic <- if (cs[["evaluable"]] > 0) cs[["monophyletic"]] / cs[["evaluable"]] else NA_real_
    for (nm in names(cg)) sm[[paste0("genus_", nm)]] <- as.integer(cg[[nm]])
    sm$genus_prop_monophyletic <- if (cg[["evaluable"]] > 0) cg[["monophyletic"]] / cg[["evaluable"]] else NA_real_
    sm_rows[[l]] <- sm
  }
  bind <- function(x) { d <- do.call(rbind, x); rownames(d) <- NULL; d }
  out <- list(species = bind(sp_rows), genus = bind(ge_rows), summary = bind(sm_rows), reference = bind(rf_rows))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(out)) {
    utils::write.csv(out[[nm]], file.path(output_dir, paste0("TABLE_barcoding_monophyly_", nm, ".csv")), row.names = FALSE)
  }
  cat("\n====================================================\n")
  cat("  Barcoding Monophyly Complete \U0001f335\n")
  cat("====================================================\n")
  cat("  Loci:                  ", nrow(out$summary), "(rooted for grades:", sum(out$summary$rooted), ")\n")
  cat("  Support cutoff:        ", support_cutoff, "\n")
  cat("  Output directory:      ", output_dir, "\n")
  cat("====================================================\n\n")
  invisible(out)
}
