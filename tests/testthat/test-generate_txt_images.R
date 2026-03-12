require(testthat)
require(quantMSImageR)

context("generate_txt_images: structure and output_txt flag")

test_that("generate_txt_images returns the expected named list", {
  skip_if_not(
    file.exists("D:/__STUDIES/021_HDM_mice/data/DESI-MRM/raw"),
    "Live acquisition data not available — skipping integration test"
  )

  result <- generate_txt_images(
    fns          = "20260309_021_M2L_9-LHS_pos04",
    data_path    = "D:/__STUDIES/021_HDM_mice/data/DESI-MRM/raw",
    image_dir    = tempdir(),
    lib_ion_path = "D:/__STUDIES/ion_library20260309.csv",
    output_txt   = FALSE
  )

  expect_named(result,
               c("combined", "combined_snr", "combined_FC",
                 "combined_NAbackground"))
  expect_s4_class(result$combined,              "quant_MSImagingExperiment")
  expect_s4_class(result$combined_snr,          "quant_MSImagingExperiment")
  expect_s4_class(result$combined_FC,           "quant_MSImagingExperiment")
  expect_s4_class(result$combined_NAbackground, "quant_MSImagingExperiment")
})

test_that("no files are written when output_txt = FALSE", {
  skip_if_not(
    file.exists("D:/__STUDIES/021_HDM_mice/data/DESI-MRM/raw"),
    "Live acquisition data not available — skipping integration test"
  )

  tmp <- file.path(tempdir(), paste0("qmsi_test_", Sys.time()))
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  generate_txt_images(
    fns          = "20260309_021_M2L_9-LHS_pos04",
    data_path    = "D:/__STUDIES/021_HDM_mice/data/DESI-MRM/raw",
    image_dir    = tmp,
    lib_ion_path = "D:/__STUDIES/ion_library20260309.csv",
    output_txt   = FALSE
  )

  expect_equal(length(list.files(tmp, recursive = TRUE)), 0)
})

test_that("average_method must be 'mean' or 'median'", {
  expect_error(
    generate_txt_images(
      fns = "dummy", data_path = ".", image_dir = tempdir(),
      lib_ion_path = "dummy.csv", average_method = "geometric"
    ),
    regexp = "'arg' should be one of"
  )
})
