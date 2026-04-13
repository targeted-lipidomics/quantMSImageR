require(testthat)
require(quantMSImageR)

context("generate_txt_images: structure and output_txt flag")

# Paths to the bundled pos04 test acquisition
extdata      <- system.file("extdata", package = "quantMSImageR")
pos04_folder <- extdata                          # data_path: contains pos04_test.raw/
pos04_name   <- "pos04_test"                     # acquisition name (folder = pos04_test.raw/)
lib_ion_path <- file.path(extdata, "ion_library_pos04.csv")

test_that("generate_txt_images returns the expected named list", {
  result <- generate_txt_images(
    fns          = pos04_name,
    data_path    = pos04_folder,
    image_dir    = tempdir(),
    lib_ion_path = lib_ion_path,
    output_txt   = FALSE,
    average_method = "median"
  )

  expect_named(result,
               c("combined", "combined_snr", "combined_FC",
                 "combined_NAbackground"))
  expect_s4_class(result$combined,              "quant_MSImagingExperiment")
  expect_s4_class(result$combined_snr,          "quant_MSImagingExperiment")
  expect_s4_class(result$combined_FC,           "quant_MSImagingExperiment")
  expect_s4_class(result$combined_NAbackground, "quant_MSImagingExperiment")
})

test_that("output_txt = FALSE writes no files", {
  tmp <- file.path(tempdir(), paste0("qmsi_test_", as.integer(Sys.time())))
  dir.create(tmp, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  generate_txt_images(
    fns          = pos04_name,
    data_path    = pos04_folder,
    image_dir    = tmp,
    lib_ion_path = lib_ion_path,
    output_txt   = FALSE
  )

  expect_equal(length(list.files(tmp, recursive = TRUE)), 0L)
})

test_that("output_txt = TRUE creates expected subdirectories", {
  tmp <- file.path(tempdir(), paste0("qmsi_txt_", as.integer(Sys.time())))
  dir.create(tmp, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  result <- generate_txt_images(
    fns          = pos04_name,
    data_path    = pos04_folder,
    image_dir    = tmp,
    lib_ion_path = lib_ion_path,
    snr_thresh   = 1,
    output_txt   = TRUE
  )

  acq_dir <- file.path(tmp, pos04_name)
  expect_true(dir.exists(acq_dir))
  subdirs <- list.dirs(acq_dir, recursive = FALSE, full.names = FALSE)
  expect_true(any(grepl("intensity_SNRfiltered", subdirs)))
  expect_true(any(grepl("intensity_raw",         subdirs)))
})

test_that("combined_snr has fewer non-NA pixels than combined (SNR filter applied)", {
  result <- generate_txt_images(
    fns          = pos04_name,
    data_path    = pos04_folder,
    image_dir    = tempdir(),
    lib_ion_path = lib_ion_path,
    snr_thresh   = 3,
    output_txt   = FALSE
  )

  n_non_na_raw <- sum(!is.na(spectraData(result$combined)[["intensity"]]))
  n_non_na_snr <- sum(!is.na(spectraData(result$combined_snr)[["intensity"]]))

  expect_lte(n_non_na_snr, n_non_na_raw)
})

test_that("invalid average_method is rejected early", {
  expect_error(
    generate_txt_images(
      fns = "dummy", data_path = ".", image_dir = tempdir(),
      lib_ion_path = "dummy.csv", average_method = "geometric"
    ),
    regexp = "'arg' should be one of"
  )
})
