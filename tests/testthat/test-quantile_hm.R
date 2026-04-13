require(testthat)
require(quantMSImageR)
require(ComplexHeatmap)

context("quantile_hm: heatmap dimensions and orientation")

make_obj <- function(n_features = 3, runs = c("s1", "s2"), n_px = 4,
                     seed = 42) {
  set.seed(seed)
  fdata <- MassDataFrame(mz   = seq_len(n_features),
                         name = paste0("lipid_", seq_len(n_features)))
  n_tot <- length(runs) * n_px
  pdata <- PositionDataFrame(
    run   = rep(runs, each = n_px),
    coord = do.call(rbind, lapply(seq_along(runs),
                                  function(i) expand.grid(x = 1:2, y = 1:2)))
  )
  ints <- matrix(runif(n_features * n_tot, 1, 100), nrow = n_features)
  as(MSImagingExperiment(spectraData = ints, featureData = fdata,
                          pixelData = pdata),
     "quant_MSImagingExperiment")
}

test_that("quantile_hm returns a Heatmap object", {
  obj <- make_obj()
  hm  <- quantile_hm(obj, quant_val = 0.5,
                     heatmap_order = c("s1", "s2"),
                     heatmap_labs  = c("A", "B"))
  expect_s4_class(hm, "Heatmap")
})

test_that("heatmap matrix is n_features x n_samples (not transposed)", {
  n_feat <- 4; runs <- c("s1", "s2", "s3")
  obj <- make_obj(n_features = n_feat, runs = runs)
  hm  <- quantile_hm(obj, quant_val = 0.75,
                     heatmap_order = runs,
                     heatmap_labs  = c("A", "B", "C"))
  mat <- hm@matrix
  expect_equal(nrow(mat), n_feat)
  expect_equal(ncol(mat), length(runs))
})

test_that("z-score values are clipped to [-1, 1]", {
  obj <- make_obj(n_features = 5, runs = c("s1", "s2", "s3", "s4"))
  hm  <- quantile_hm(obj, quant_val = 0.5,
                     heatmap_order = c("s1", "s2", "s3", "s4"),
                     heatmap_labs  = c("A", "A", "B", "B"))
  mat <- hm@matrix
  expect_true(all(mat >= -1 & mat <= 1, na.rm = TRUE))
})

test_that("quantile_hm works without heatmap_order / heatmap_labs", {
  obj <- make_obj()
  expect_s4_class(quantile_hm(obj, quant_val = 0.5), "Heatmap")
})
