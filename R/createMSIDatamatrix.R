library(Cardinal)
library(dplyr)
library(tibble)

setGeneric("createMSIDatamatrix", function(MSIobject, ...) standardGeneric("createMSIDatamatrix"))

#' Function to create data matrix from MSI object
#' @import Cardinal
#' @import dplyr
#' @import tibble
#' @include setClasses.R
#'
#' @param MSIobject MSI object from Cardinal, pData() to include sample_ID.
#' @param val_slot character defining slot name to normalise - takes "intensity" as default
#' @param inputNA Whether to convert 0's in matrix to NA (default = T)
#' @param roi_header Header in pData pertaining to ROIs to average. Set to NA to skip generating average df
#' @return MSIobject with slots updated for i) matrix of average ng/pixel of m/z (rows = m/z and cols = cal level) in tissue ROIs ii) sample/ROI metadata
#'
#' @export
setMethod("createMSIDatamatrix", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "intensity", inputNA = T, roi_header = NA){

            # Subset pixels in ROIs only (based on roi_header)
            if(is.na(roi_header)){
              pData(MSIobject)$ROI = 1:nrow(pData(MSIobject))
            } else{
              MSIobject = MSIobject[, which(!is.na(pData(MSIobject)[[roi_header]]))]
            }

            # Update pixel data
            pixel_df = data.frame(pData(MSIobject)) |>
              subset(!is.na(ROI)) |>
              tibble::rownames_to_column("pixel_ind") |>
              dplyr::mutate(pixel_ind = sprintf("pixel_%s", pixel_ind))


            # All pixel df
            all_pixel_df = data.frame( sapply(array(1:nrow(fData(MSIobject))), FUN = function(x) spectraData(MSIobject)[[val_slot]][x, ] ) ) |>
              dplyr::mutate(pixel_ind = pixel_df$pixel_ind)
            colnames(all_pixel_df) = c(fData(MSIobject)$name, "pixel_ind")

            all_pixel_df = dplyr::left_join(all_pixel_df, pixel_df, by = "pixel_ind") |> dplyr::select(dplyr::any_of(c(fData(MSIobject)$name, "pixel_ind", roi_header)))


            if(inputNA){
              all_pixel_df <- replace(all_pixel_df, all_pixel_df==0, NA)
            }


            # Create average ROI df
            if(!is.na(roi_header)){

              #### THIS NEEDS FIXING TO SUMMARISE EACH FEATURE INDEPENDENTLY!!!!
              ave_df = all_pixel_df |>
                dplyr::group_by(dplyr::across(dplyr::all_of(roi_header))) |>
                dplyr::summarise(dplyr::across(dplyr::any_of(c(fData(MSIobject)$name)), \(x) mean(x, na.rm=T))) |>
                tibble::column_to_rownames(roi_header)

              if(inputNA){
                ave_df <- replace(ave_df, ave_df==0, NA)
              }

              MSIobject@tissueInfo@roi_average_matrix = ave_df
            }

            all_pixel_df = all_pixel_df |>
              tibble::column_to_rownames("pixel_ind") |>
              dplyr::select(dplyr::any_of(c(fData(MSIobject)$name)))

            MSIobject@tissueInfo@all_pixel_matrix = all_pixel_df
            MSIobject@tissueInfo@sample_metadata = pixel_df

            return(MSIobject)
          })
