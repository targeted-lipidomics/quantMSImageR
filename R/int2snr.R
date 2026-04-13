setGeneric("int2snr", function(MSIobject, ...) standardGeneric("int2snr"))

#' Function to convert the intensity values to SNR per pixel based on same transitions in noise/background pixels.Run after IS normalization.
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject quant_MSImagingExperiment - including pData(MSIobject)$sample_type == Noise
#' @param noise character in pData(MSIobject)$sample_type which indicates background / noise pixels
#' @param tissue character in pData(MSIobject)$sample_type which indicates tissue pixels to calculate SNR for
#' @param val_slot character defining slot name to normalise - takes "intensity" as default
#' @param snr_thresh value indicating the minimum SNR value to accept (below this value SNR = 0)
#' @return MSIobject with intensity values replaced with SNR values
#'
#' @export
setMethod("int2snr", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "response", noise = "Noise", tissue = "Tissue", snr_thresh = 3,
                   sample_type = "sample_type", average = c("mean", "median"), ...){
            average <- match.arg(average)

            if(!any(pData(MSIobject)[[sample_type]] == noise)){
              print("No noise pixels. Return same values")
              return(MSIobject)
            }

            #Set noise and tissue pixels
            noise_pixels = which(pData(MSIobject)[[sample_type]] == noise)
            tissue_pixels = which(pData(MSIobject)[[sample_type]] == tissue)

            spectra(MSIobject, "snr") = matrix(nrow = nrow(MSIobject), ncol = ncol(MSIobject))

            # Iterate over features in study
            for(mz_ind in 1:nrow(fData(MSIobject))){

              #print(sprintf("mz - %s", mz_ind))

              # Save noise response vector
              noise_vec = spectraData(MSIobject)[[val_slot]][mz_ind, noise_pixels]

              #Deal with 0 values in noise vector - 10% of lowest
              noise_vec[which(noise_vec ==0)] = NA
              mv_impute = min(noise_vec, na.rm = T) / 10
              noise_vec[which(is.na(noise_vec))] = mv_impute

              noise_level <- if (average == "mean") mean(noise_vec) else median(noise_vec)

              # Calculate S/N of tissue pixels
              vec = spectraData(MSIobject)[[val_slot]][mz_ind, ]
              snr =  vec / noise_level
              snr =  ifelse(snr < snr_thresh, NA, snr)
              snr[which(!1:length(vec) %in% tissue_pixels)] = NA
              spectra(MSIobject, "snr")[mz_ind, ] = snr

            }

            return(MSIobject)
          })
