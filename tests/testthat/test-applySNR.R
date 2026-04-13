require(testthat)
require(quantMSImageR)

context("applySNR: intensity set to NA where SNR is NA")

make_obj <- function() {
  fdata <- MassDataFrame(mz = 1:2, name = c("feat_1", "feat_2"))
  pdata <- PositionDataFrame(run   = rep("s1", 4),
                             coord = expand.grid(x = 1:2, y = 1:2))
  ints  <- matrix(c(10, 20, 30, 40, 100, 200, 300, 400), nrow = 2, byrow = TRUE)
  obj   <- as(MSImagingExperiment(spectraData = ints, featureData = fdata,
                                   pixelData = pdata),
              "quant_MSImagingExperiment")
  # Inject an SNR slot with some NAs
  spectra(obj, "snr") <- matrix(c(NA, 2, NA, 4,
                                   5, NA, 7, NA), nrow = 2, byrow = TRUE)
  obj
}

test_that("intensity is NA wherever SNR is NA", {
  obj    <- make_obj()
  result <- applySNR(obj, val_slot = "intensity")
  int_r  <- spectraData(result)[["intensity"]]
  snr_r  <- spectraData(result)[["snr"]]

  expect_true(all(is.na(int_r[is.na(snr_r)])))
})

test_that("intensity is preserved wherever SNR is not NA", {
  obj    <- make_obj()
  result <- applySNR(obj, val_slot = "intensity")
  int_o  <- spectraData(obj)[["intensity"]]
  int_r  <- spectraData(result)[["intensity"]]
  snr_r  <- spectraData(result)[["snr"]]

  expect_equal(int_r[!is.na(snr_r)], int_o[!is.na(snr_r)])
})
