require(testthat)
require(quantMSImageR)
require(ggplot2)

context("imageR: ion image generation")

make_img_obj <- function(n_features = 2, seed = 1) {
  set.seed(seed)
  fdata <- MassDataFrame(mz   = seq_len(n_features),
                         name = paste0("lipid_", seq_len(n_features)))
  pdata <- PositionDataFrame(run       = rep("s1", 4),
                              coord     = expand.grid(x = 1:2, y = 1:2),
                              sample_ID = rep("tissue_pixels", 4))
  ints  <- matrix(abs(rnorm(n_features * 4, mean = 100, sd = 20)),
                  nrow = n_features)
  as(MSImagingExperiment(spectraData = ints, featureData = fdata,
                          pixelData = pdata),
     "quant_MSImagingExperiment")
}

test_that("imageR returns a ggplot with scale = 'suppress'", {
  obj <- make_img_obj()
  p   <- imageR(obj, feat_ind = 1, scale = "suppress")
  expect_s3_class(p, "gg")
})

test_that("imageR returns a ggplot with scale = 'sqrt'", {
  obj <- make_img_obj()
  p   <- imageR(obj, feat_ind = 1, scale = "sqrt")
  expect_s3_class(p, "gg")
})

test_that("imageR returns a data.frame / matrix in text_image mode", {
  obj <- make_img_obj()
  p   <- imageR(obj, feat_ind = 1, scale = "suppress", text_image = TRUE)
  expect_true(is.data.frame(p) || is.matrix(p))
})

test_that("imageR selects different features correctly", {
  obj <- make_img_obj(n_features = 3)
  p1  <- imageR(obj, feat_ind = 1, scale = "suppress")
  p2  <- imageR(obj, feat_ind = 2, scale = "suppress")
  expect_false(identical(p1$data$response, p2$data$response))
})
