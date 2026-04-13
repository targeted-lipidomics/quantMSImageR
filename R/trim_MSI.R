#' Remove pure-noise border rows and columns from an MSI object
#'
#' Drops any x-column or y-row in which **every** pixel is labelled as
#' `"noise_pixels"` in `pData(MSI_data)$sample_name`. This trims empty
#' scan borders that arise from the Waters acquisition geometry without
#' affecting tissue or mixed-content rows/columns.
#'
#' @import Cardinal
#' @import dplyr
#'
#' @param MSI_data A `quant_MSImagingExperiment` object whose `pData()` contains
#'   columns `x`, `y`, and `sample_name` (populated by `makeFactor()`).
#'
#' @return A subset `quant_MSImagingExperiment` with pure-noise border pixels
#'   removed.
#'
#' @export
trim_MSI = function(MSI_data){
  pd <- as.data.frame(pData(MSI_data))

  # Identify x columns that are entirely noise
  bad_x <- pd |>
    dplyr::group_by(x) |>
    dplyr::summarise(all_noise = all(sample_name == "noise_pixels")) |>
    dplyr::filter(all_noise) |>
    dplyr::pull(x)

  # Identify y rows that are entirely noise
  bad_y <- pd |>
    dplyr::group_by(y) |>
    dplyr::summarise(all_noise = all(sample_name == "noise_pixels")) |>
    dplyr::filter(all_noise) |>
    dplyr::pull(y)

  # Pixels to KEEP
  keep_pixels <- which(
    !(pd$x %in% bad_x | pd$y %in% bad_y)
  )

  return(MSI_data[, keep_pixels])
}
