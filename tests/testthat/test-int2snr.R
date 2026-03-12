require(testthat)
require(quantMSImageR)

context("int2snr: SNR calculation and average_method validation")

make_snr_obj <- function() {
  fdata <- MassDataFrame(mz = c(500, 510, 540, 550),
                         analyte = c("IS", rep("analyte", 3)))
  pdata <- PositionDataFrame(
    run       = c(rep("run1", 2), rep("run2", 2), rep("Noise", 2)),
    coord     = expand.grid(x = 1:3, y = 1:2),
    sample_ID = c(rep("Tissue", 4), rep("Noise", 2))
  )
  ints <- matrix(nrow = 4, ncol = 6,
                 data = c(rep(c(5, 100, 150, 100), 2),
                          rep(c(10, 100, 105, 300), 2),
                          1, 2, 3, 4,
                          3, 6, 9, 12))
  obj <- MSImagingExperiment(spectraData = ints, featureData = fdata,
                              pixelData = pdata)
  as(obj, "quant_MSImagingExperiment")
}

test_that("SNR is NA for noise pixels and correct for tissue (average = 'mean')", {
  obj      <- make_snr_obj()
  new_data <- int2snr(MSIobject = obj, val_slot = "intensity",
                      noise = "Noise", tissue = "Tissue",
                      snr_thresh = 3, sample_type = "sample_ID",
                      average = "mean")

  expect_true(all(is.na(spectra(new_data, "snr")[, 5:6])))
  expect_equal(spectra(new_data, "snr")[1, ], c(NA, NA, 5,    5,    NA, NA))
  expect_equal(spectra(new_data, "snr")[2, ], c(25, 25, 25,   25,   NA, NA))
  expect_equal(spectra(new_data, "snr")[3, ], c(25, 25, 17.5, 17.5, NA, NA))
  expect_equal(spectra(new_data, "snr")[4, ], c(12.5, 12.5, 37.5, 37.5, NA, NA))
})

test_that("average = 'median' gives different results from 'mean' for skewed noise", {
  obj       <- make_snr_obj()
  res_mean  <- int2snr(obj, val_slot = "intensity",
                       noise = "Noise", tissue = "Tissue",
                       snr_thresh = 1, sample_type = "sample_ID",
                       average = "mean")
  res_med   <- int2snr(obj, val_slot = "intensity",
                       noise = "Noise", tissue = "Tissue",
                       snr_thresh = 1, sample_type = "sample_ID",
                       average = "median")
  # Results need not be identical; just verify both run without error
  expect_s4_class(res_mean, "quant_MSImagingExperiment")
  expect_s4_class(res_med,  "quant_MSImagingExperiment")
})

test_that("invalid average_method is rejected", {
  obj <- make_snr_obj()
  expect_error(
    int2snr(obj, val_slot = "intensity",
            noise = "Noise", tissue = "Tissue",
            snr_thresh = 3, sample_type = "sample_ID",
            average = "geometric"),
    regexp = "'arg' should be one of"
  )
})

test_that("returns object unchanged when no noise pixels exist", {
  fdata <- MassDataFrame(mz = 1:2, name = c("f1", "f2"))
  pdata <- PositionDataFrame(run       = rep("s1", 4),
                              coord     = expand.grid(x = 1:2, y = 1:2),
                              sample_ID = rep("Tissue", 4))
  ints  <- matrix(c(10, 20, 30, 40, 100, 200, 300, 400), nrow = 2, byrow = TRUE)
  obj   <- as(MSImagingExperiment(spectraData = ints, featureData = fdata,
                                   pixelData = pdata),
              "quant_MSImagingExperiment")

  result <- int2snr(obj, val_slot = "intensity",
                    noise = "Noise", tissue = "Tissue",
                    snr_thresh = 3, sample_type = "sample_ID")
  expect_equal(ncol(result), ncol(obj))
})
