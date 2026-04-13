setGeneric("back2NA", function(MSIobject, ...) standardGeneric("back2NA"))

#' Set background pixel intensities to NA
#'
#' Replaces the spectral values of background/noise pixels with `NA` in a
#' chosen intensity slot. Tissue pixels are left unchanged.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` object whose `pData()` contains
#'   a column identifying pixel type.
#' @param val_slot Character. Name of the spectra slot to modify (default `"intensity"`).
#' @param background Character. Value in `pData(MSIobject)[[sample_type]]` that
#'   labels background/noise pixels (default `"Noise"`).
#' @param tissue Character. Value that labels tissue pixels (default `"Tissue"`).
#'   Currently unused but kept for API symmetry with `int2snr`.
#' @param sample_type Character. Column name in `pData()` holding pixel-type
#'   labels (default `"sample_type"`).
#'
#' @return The input `quant_MSImagingExperiment` with background pixels set to
#'   `NA` in `val_slot`. Returns the object unchanged if no background pixels are
#'   found.
#'
#' @export
setMethod("back2NA", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "response", background = "Noise", tissue = "Tissue",sample_type = "sample_type", ...){

            if(!any(pData(MSIobject)[[sample_type]] == background)){
              print("No background pixels. Return same values")
              return(MSIobject)
            }

            #Set background and tissue pixels
            background_pixels = which(pData(MSIobject)[[sample_type]] == background)

            # Iterate over features in study
            for(mz_ind in 1:nrow(fData(MSIobject))){

              # Save background response vector
              spectraData(MSIobject)[[val_slot]][mz_ind, background_pixels] = NA

            }

            return(MSIobject)
          })
