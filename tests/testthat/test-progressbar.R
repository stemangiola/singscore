library(testthat)
library(HDF5Array)
library(singscore)
library(DelayedArray)
library(DelayedMatrixStats)

test_that("rankGenes with progressbar parameter works on DelayedMatrix", {
  # Create a dummy HDF5-backed matrix
  h5_file <- tempfile(fileext = ".h5")
  h5_path <- "assay"
  mat <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10,
                dimnames = list(paste0("gene", 1:100), paste0("sample", 1:10)))
  h5write(mat, file = h5_file, name = h5_path)
  dm <- HDF5Array(h5_file, h5_path)
  
  # Test with workers = 1 (default)
  stable <- c("gene1", "gene5", "gene10")
  # Ensure stable genes exist in the matrix
  stable <- intersect(stable, rownames(dm))
  ranked_dm_1 <- rankGenes(dm, stableGenes = stable, workers = 1)
  expect_true(is.matrix(ranked_dm_1))
  expect_equal(dim(ranked_dm_1), dim(mat))
  
  # Test with workers = 2 (enables progress bar)
  ranked_dm_2 <- rankGenes(dm, stableGenes = stable, workers = 2)
  expect_true(is.matrix(ranked_dm_2))
  expect_equal(dim(ranked_dm_2), dim(mat))
  
  # Results should be identical regardless of worker count
  expect_equal(ranked_dm_1, ranked_dm_2)
  
  # Clean up
  file.remove(h5_file)
})

test_that("rankGenes with different worker counts produces same results", {
  # Create a dummy HDF5-backed matrix
  h5_file <- tempfile(fileext = ".h5")
  h5_path <- "assay"
  mat <- matrix(rnorm(50 * 5), nrow = 50, ncol = 5,
                dimnames = list(paste0("gene", 1:50), paste0("sample", 1:5)))
  h5write(mat, file = h5_file, name = h5_path)
  dm <- HDF5Array(h5_file, h5_path)
  
  stable <- c("gene1", "gene10", "gene20")
  # Ensure stable genes exist in the matrix
  stable <- intersect(stable, rownames(dm))
  
  # Test with different worker counts
  ranked_1 <- rankGenes(dm, stableGenes = stable, workers = 1)
  ranked_2 <- rankGenes(dm, stableGenes = stable, workers = 2)
  
  # Results should be identical regardless of worker count
  expect_equal(ranked_1, ranked_2)
  
  # Clean up
  file.remove(h5_file)
})
