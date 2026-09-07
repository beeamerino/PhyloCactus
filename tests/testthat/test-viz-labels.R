test_that(".collapse_weak_support_nodes() normalizes percentage-scale (0-100) support and collapses weak nodes accordingly (real function call)", {
  skip_if_not_installed("ape")

  # Node labels: root = 90, (A,B) = 50, (C,D) = 100 -- percentage scale (max > 1.0).
  tree <- ape::read.tree(text = "((A:1,B:1)50:1,(C:1,D:1)100:1)90;")
  expect_true(ape::is.binary(tree))

  collapsed <- .collapse_weak_support_nodes(tree, collapse_cutoff = 0.70)

  # (A,B) support normalizes to 0.50, below the 0.70 cutoff: this internal node must be
  # collapsed into the root, turning the root into a soft polytomy (A, B, and the (C,D) clade).
  expect_false(ape::is.binary(collapsed))
  expect_equal(collapsed$Nnode, 2L)

  # (C,D) support normalizes to 1.00, at/above cutoff: this clade must remain intact.
  expect_true(ape::is.monophyletic(collapsed, c("C", "D")))
})

test_that(".collapse_weak_support_nodes() does not double-normalize support values already on a 0-1 scale (real function call)", {
  skip_if_not_installed("ape")

  # Same topology and same underlying support (50% / 100% / 90%), but already expressed on a
  # 0-1 scale. If the normalization step incorrectly divided these by 100 again, every node
  # would fall below any reasonable cutoff and (C,D) would also be collapsed.
  tree <- ape::read.tree(text = "((A:1,B:1)0.50:1,(C:1,D:1)1.00:1)0.90;")

  collapsed <- .collapse_weak_support_nodes(tree, collapse_cutoff = 0.70)

  expect_false(ape::is.binary(collapsed))
  expect_equal(collapsed$Nnode, 2L)
  expect_true(ape::is.monophyletic(collapsed, c("C", "D")))
})

test_that(".collapse_weak_support_nodes() warns and preserves nodes with missing or non-numeric support labels rather than silently treating them as supported (real function call)", {
  skip_if_not_installed("ape")

  # (A,B) has a non-numeric support label ("XYZ"); the corrected package behavior must warn
  # explicitly and leave this node uncollapsed (never silently pass it through as "supported").
  tree <- ape::read.tree(text = "((A:1,B:1)XYZ:1,(C:1,D:1)90:1)80;")

  expect_warning(
    collapsed <- .collapse_weak_support_nodes(tree, collapse_cutoff = 0.70),
    "missing or non-numeric support"
  )

  # The (A,B) clade must survive uncollapsed despite its unparseable support label.
  expect_true(ape::is.monophyletic(collapsed, c("A", "B")))
})
