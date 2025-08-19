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


test_that("rankExpr gives identical results for matrix vs HDF5-backed DelayedMatrix", {
  testthat::skip_if_not_installed("HDF5Array")

  emat <- SummarizedExperiment::assay(toy_expr_se)
  ranked_mat <- rankExpr(emat)

  h5_file <- tempfile(fileext = ".h5")
  on.exit(unlink(h5_file), add = TRUE)
  dm <- HDF5Array::writeHDF5Array(emat, filepath = h5_file, name = "expr")
  ranked_dm <- rankExpr(dm)

  expect_equal(ranked_dm, ranked_mat)
})

test_that("rankExprStable matches rankExprStable_delayed for matrix input", {
  emat <- SummarizedExperiment::assay(toy_expr_se)
  st <- rownames(emat)[seq_len(min(10, nrow(emat)))]

  ranked_old <- rankExprStable(emat, stgenes = st)
  ranked_new <- rankExprStable_delayed(emat, stgenes = st)

  expect_equal(ranked_new, ranked_old)
  expect_true(isTRUE(attr(ranked_new, "stable")))
  expect_true(isTRUE(attr(ranked_old, "stable")))
  expect_true(all(ranked_new >= 0 & ranked_new <= 1))
  expect_true(all(ranked_old >= 0 & ranked_old <= 1))
})


