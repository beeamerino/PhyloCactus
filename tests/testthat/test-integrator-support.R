test_that(".collapse_weak_support_nodes normalizes percentage support to the 0-1 scale (real function call)", {
  skip_if_not_installed("ape")

  # Node labels on a percentage scale. With collapse_cutoff = 0.70 the 40 percent node
  # must be collapsed and the 90, 95 and 100 percent nodes must survive, which is only
  # true if the function divides by 100 first. The weak node is deliberately an internal
  # node rather than the root: collapsing works by zeroing the incident edge and calling
  # ape::di2multi(), and the root has no incident edge, so a weak root cannot be
  # collapsed by construction (see the dedicated test below).
  tree <- ape::read.tree(text = "(((A:1,B:1)40:1,(C:1,D:1)95:1)90:1,E:2)100;")

  collapsed <- .collapse_weak_support_nodes(tree, collapse_cutoff = 0.70)

  expect_s3_class(collapsed, "phylo")
  expect_lt(collapsed$Nnode, tree$Nnode)
  expect_setequal(collapsed$tip.label, tree$tip.label)
})

test_that(".collapse_weak_support_nodes does not double-normalize 0-1 support (real function call)", {
  skip_if_not_installed("ape")

  # Already on the 0-1 scale. Dividing again by 100 would push every node below the
  # cutoff and collapse the whole topology, so the node count must be preserved here.
  tree <- ape::read.tree(text = "((A:1,B:1)0.95:1,(C:1,D:1)0.90:1)0.85;")

  collapsed <- .collapse_weak_support_nodes(tree, collapse_cutoff = 0.70)

  expect_equal(collapsed$Nnode, tree$Nnode)
  expect_setequal(collapsed$tip.label, tree$tip.label)
})

test_that(".collapse_weak_support_nodes warns on non-numeric support instead of treating it as strong (real function call)", {
  skip_if_not_installed("ape")

  # which(NA < cutoff) silently drops NA, so a node with an unreadable label would be
  # kept as if it were well supported. The function must say so out loud.
  tree <- ape::read.tree(text = "((A:1,B:1)NA_LABEL:1,(C:1,D:1)0.95:1)0.90;")

  expect_warning(
    .collapse_weak_support_nodes(tree, collapse_cutoff = 0.70),
    "missing or non-numeric support"
  )
})

test_that(".collapse_weak_support_nodes collapses every node below the cutoff (real function call)", {
  skip_if_not_installed("ape")

  weak <- ape::read.tree(text = "((A:1,B:1)0.10:1,(C:1,D:1)0.20:1)0.30;")
  strong <- ape::read.tree(text = "((A:1,B:1)0.99:1,(C:1,D:1)0.98:1)0.97;")

  expect_lt(.collapse_weak_support_nodes(weak, collapse_cutoff = 0.70)$Nnode, weak$Nnode)
  expect_equal(.collapse_weak_support_nodes(strong, collapse_cutoff = 0.70)$Nnode, strong$Nnode)
})

test_that("a weakly supported root is not collapsed, because it has no incident edge (real function call)", {
  skip_if_not_installed("ape")

  # Collapsing works by setting the incident edge length to zero and letting
  # ape::di2multi() merge the node. The root has no incident edge, so a root below the
  # cutoff survives silently. This test pins that limitation rather than asserting a
  # behaviour the function does not have: a root-level polytomy has to be handled by
  # rooting, not by this function.
  tree <- ape::read.tree(text = "((A:1,B:1)0.95:1,(C:1,D:1)0.90:1)0.30;")

  collapsed <- .collapse_weak_support_nodes(tree, collapse_cutoff = 0.70)

  expect_equal(collapsed$Nnode, tree$Nnode)
})
