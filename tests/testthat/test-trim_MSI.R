require(testthat)
require(quantMSImageR)

context("trim_MSI: pure-noise border pixels removed")

make_obj <- function(sample_names_vec, xs, ys) {
  fdata <- MassDataFrame(mz = 1:2, name = c("f1", "f2"))
  coord <- data.frame(x = xs, y = ys)
  pdata <- PositionDataFrame(run         = rep("s1", length(xs)),
                              coord       = coord,
                              sample_name = sample_names_vec)
  ints  <- matrix(rep(1, 2 * length(xs)), nrow = 2)
  as(MSImagingExperiment(spectraData = ints, featureData = fdata,
                          pixelData = pdata),
     "quant_MSImagingExperiment")
}

test_that("columns where all pixels are noise are removed", {
  # x = 1: all noise; x = 2,3: tissue
  xs  <- c(1, 2, 3, 1, 2, 3)
  ys  <- c(1, 1, 1, 2, 2, 2)
  ids <- c("noise_pixels", "tissue_pixels", "tissue_pixels",
           "noise_pixels", "tissue_pixels", "tissue_pixels")
  obj    <- make_obj(ids, xs, ys)
  result <- trim_MSI(obj)

  expect_equal(ncol(result), 4)
  expect_true(all(pData(result)$x != 1))
})

test_that("rows where all pixels are noise are removed", {
  # y = 2: all noise; y = 1: tissue
  xs  <- c(1, 2, 1, 2)
  ys  <- c(1, 1, 2, 2)
  ids <- c("tissue_pixels", "tissue_pixels",
           "noise_pixels",  "noise_pixels")
  obj    <- make_obj(ids, xs, ys)
  result <- trim_MSI(obj)

  expect_equal(ncol(result), 2)
  expect_true(all(pData(result)$y != 2))
})

test_that("object is unchanged when no pure-noise borders exist", {
  xs  <- c(1, 2, 1, 2)
  ys  <- c(1, 1, 2, 2)
  ids <- rep("tissue_pixels", 4)
  obj    <- make_obj(ids, xs, ys)
  result <- trim_MSI(obj)

  expect_equal(ncol(result), 4)
})
