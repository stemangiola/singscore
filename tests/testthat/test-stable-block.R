test_that(".rank_block_by_stable_genes computes expected ranks on a small block", {
  # small toy block: 3 genes (rows) x 2 samples (cols)
  block <- matrix(c(1, 3, 2,
                    5, 4, 6), nrow = 3, byrow = FALSE)
  rownames(block) <- c("g1", "g2", "g3")
  colnames(block) <- c("s1", "s2")
  # stable genes are g1 and g3 -> indices 1 and 3
  st_idx <- c(1L, 3L)
  # call internal helper directly
  res <- singscore:::`.rank_block_by_stable_genes`(block, st_idx)
  # expected raw ranks per column (counts of stable s < x + 1)
  # col1 (s = [1,2]): x=[1,3,2] -> [1,3,2]
  # col2 (s = [5,6]): x=[5,4,6] -> [1,1,2]
  exp <- matrix(c(1, 3, 2,
                  1, 1, 2), nrow = 3, byrow = FALSE,
                dimnames = list(rownames(block), colnames(block)))
  expect_equal(unname(res), unname(exp))
})

test_that("rankGenes stable-genes normalization matches block helper", {
  mat <- matrix(c(1, 3, 2,
                  5, 4, 6), nrow = 3, byrow = FALSE)
  rownames(mat) <- c("g1", "g2", "g3")
  colnames(mat) <- c("s1", "s2")
  # stable genes g1 and g3
  st <- c("g1", "g3")
  # expected normalized = raw_ranks / (|stable|+1) = / 3
  raw <- singscore:::`.rank_block_by_stable_genes`(mat, c(1L, 3L))
  exp_norm <- raw / 3
  ranked <- rankGenes(mat, stableGenes = st)
  expect_true(isTRUE(attr(ranked, "stable")))
  expect_true(isTRUE(all.equal(unname(ranked), unname(exp_norm), check.attributes = FALSE)))
})


