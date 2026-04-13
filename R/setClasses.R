#' calibrationInfo
#'
#' Class to store information about the calibration spots
#'
#' @name calibrationInfo
#' @export
calibrationInfo = setClass("calibrationInfo",
         slots = c(
           cal_metadata = "data.frame",
           cal_list = "list",
           cal_response_data = "data.frame",
           r2_df = "data.frame"
         )
)

#' tissueInfo
#'
#' Class to store information about the calibration spots
#'
#' @name tissueInfo
#' @export
tissueInfo = setClass("tissueInfo",
                           slots = c(
                             roi_average_matrix = "data.frame",
                             all_pixel_matrix = "data.frame",
                             sample_metadata = "data.frame",
                             feature_metadata = "data.frame"
                           )
)


#' quant_MSImagingExperiment
#'
#' Class containting calibration metadata and MSImaging experiment
#'
#' @import Cardinal
#' @name quant_MSImagingExperiment
#' @export
quant_MSImagingExperiment = setClass("quant_MSImagingExperiment",
         contains = 'MSImagingExperiment',
         slots = c(
           calibrationInfo = "calibrationInfo",
           tissueInfo = "tissueInfo"
         )
)
