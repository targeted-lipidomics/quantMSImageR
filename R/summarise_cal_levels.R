library(Cardinal)
library(dplyr)

setGeneric("summarise_cal_levels", function(MSIobject, ...) standardGeneric("summarise_cal_levels"))

#' Function to calculate the mean response or intensity per pixel for the ROI at each calibration level across all calibration replicates (ng/pixel).
#' @import Cardinal
#' @import dplyr
#' @include setClasses.R
#'
#' @param MSIobject MSI object from Cardinal
#' @param cal_metadata dataframe containing calibration metdata info - including "lipid" = feature name in fData(), "identifier" header to map pData, and "amount_pg" relating to amount of std at each spot.
#' @param val_slot character defining slot name to normalise - takes "intensity" as default
#' @param cal_header Header in pixel metadata to select calibration data from. Default = "sample_type"
#' @param cal_label Label in pixel metadata under the `cal_header` which corresponds to calibration data. Default = "Cal".
#' @param id header in calibration metadata and pData to map (defaults to "identifier") and label unique ROIs
#' @return MSIobject with slots updated for i) matrix of average ng/pixel of m/z (rows = m/z and cols = cal level) ii) list of pixel counts per cal level
#'
#' @export
setMethod("summarise_cal_levels", "quant_MSImagingExperiment",
          function(MSIobject, cal_metadata, val_slot = "response", cal_header = "sample_type", cal_label = "Cal", id = "identifier"){

            MSIobject@calibrationInfo@cal_metadata = cal_metadata

            # create pixel data to associate pixel indices to cal levels
            pixel_data = data.frame(pData(MSIobject)) |>
              dplyr::mutate(pixel_ind = 1:nrow(pData(MSIobject))) |>
              subset(sample_type == cal_label) |>
              subset(!is.na(.[[id]]))

            # Create output response df
            response_df = tibble::tibble(cal_spot = unique(pixel_data[[id]]),
                                 response_perpixel = NA,
                                 pixels = NA) |>
              dplyr::left_join(MSIobject@calibrationInfo@cal_metadata, by=c("cal_spot" = id)) |>
              dplyr::select(dplyr::any_of(c("cal_spot", "response_perpixel", "pixels", "level", "lipid", "amount_pg")))

            for(i in 1:nrow(response_df)){

              # Select feature
              lipid_n = response_df$lipid[i]
              lipid_ind = which(fData(MSIobject)$name == lipid_n)

              # select pixels
              cal_n = response_df$cal_spot[i]
              inds = which(pixel_data[[id]] == cal_n)
              pixels = as.numeric(pixel_data$pixel_ind[inds])


              response_df$pixels[i] = length(pixels)

                ints = spectraData(MSIobject)[[val_slot]][lipid_ind, pixels]
              ints = replace(ints, ints ==0, NA)

              response_df$response_perpixel[i] = mean(ints, na.rm=T)

            }

            response_df = mutate(response_df, pg_perpixel = amount_pg / pixels)


            MSIobject@calibrationInfo@cal_response_data = response_df


            return(MSIobject)
          })
