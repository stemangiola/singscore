test_that("rankGenes works on HDF5-backed DelayedMatrix", {
  testthat::skip_if_not_installed("HDF5Array")
  testthat::skip_if_not_installed("DelayedArray")

  # create small matrix and persist as HDF5-backed DelayedMatrix
  set.seed(1)
  mat <- matrix(runif(5 * 4), nrow = 5)
  rownames(mat) <- paste0("g", seq_len(nrow(mat)))
  colnames(mat) <- paste0("s", seq_len(ncol(mat)))

  h5_file <- tempfile(fileext = ".h5")
  on.exit(unlink(h5_file), add = TRUE)
  dm <- HDF5Array::writeHDF5Array(mat, filepath = h5_file, name = "expr")

  # default ranking (no stable genes)
  ranked <- rankGenes(dm)
  expect_equal(dim(ranked), dim(mat))
  expect_equal(rownames(ranked), rownames(mat))
  expect_equal(colnames(ranked), colnames(mat))
  expect_true(all(ranked >= 1 & ranked <= nrow(mat)))

  # stable-genes ranking on DelayedMatrix (uses on-disk friendly backend)
  stable <- rownames(mat)[1:2]
  ranked_stable <- rankGenes(dm, stableGenes = stable)
  expect_equal(dim(ranked_stable), dim(mat))
  expect_true(isTRUE(attr(ranked_stable, "stable")))
  expect_true(all(ranked_stable >= 0 & ranked_stable <= 1))

  # scoring should work on ranked output
  up <- rownames(mat)[1:2]
  down <- rownames(mat)[5]
  sc <- simpleScore(ranked, upSet = up, downSet = down)
  expect_s3_class(sc, "data.frame")
  expect_equal(nrow(sc), ncol(mat))
  expect_true(all(c("TotalScore", "UpScore", "DownScore") %in% colnames(sc)))
})


test_that("scores are identical for matrix vs HDF5-backed DelayedMatrix", {
  testthat::skip_if_not_installed("HDF5Array")

  # Use package toy data
  emat <- SummarizedExperiment::assay(toy_expr_se)
  # Matrix path
  ranked_mat <- rankGenes(emat)
  sc_mat <- simpleScore(ranked_mat, upSet = toy_gs_up, downSet = toy_gs_dn)

  # Delayed (HDF5-backed) path
  h5_file <- tempfile(fileext = ".h5")
  on.exit(unlink(h5_file), add = TRUE)
  dm <- HDF5Array::writeHDF5Array(emat, filepath = h5_file, name = "expr")
  ranked_dm <- rankGenes(dm)
  sc_dm <- simpleScore(ranked_dm, upSet = toy_gs_up, downSet = toy_gs_dn)

  # Compare
  expect_equal(sc_dm, sc_mat)
})

test_that("stable-genes ranking matches for matrix vs DelayedMatrix", {
  testthat::skip_if_not_installed("HDF5Array")

  emat <- SummarizedExperiment::assay(toy_expr_se)
  st <- rownames(emat)[seq_len(min(10, nrow(emat)))]

  ranked_mat_st <- rankGenes(emat, stableGenes = st)
  sc_mat_st <- simpleScore(ranked_mat_st, upSet = toy_gs_up, downSet = toy_gs_dn)

  h5_file <- tempfile(fileext = ".h5")
  on.exit(unlink(h5_file), add = TRUE)
  dm <- HDF5Array::writeHDF5Array(emat, filepath = h5_file, name = "expr")

  ranked_dm_st <- rankGenes(dm, stableGenes = st)
  sc_dm_st <- simpleScore(ranked_dm_st, upSet = toy_gs_up, downSet = toy_gs_dn)

  expect_equal(sc_dm_st, sc_mat_st)
})


