library(Cardinal)

setGeneric("int2conc", function(MSIobject, ...) standardGeneric("int2conc"))

#' Function to update intensity with concentration values
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject quant_MSImagingExperiment object
#' @param pixels Label in pixel metadata the pixels to quantify
#' @param val_slot character defining slot name to normalise - takes "intensity" as default
#' @return MSIobject with intensity values replaced with concentration values (ng/pixel)
#'
#' @export
setMethod("int2conc", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "response", pixel_header = "sample_type", pixels = "Tissue"){

            cal_list = MSIobject@calibrationInfo@cal_list
            pixel_inds = which(pData(MSIobject)[[pixel_header]] %in% pixels)
            MSIobject = MSIobject[, pixel_inds]

            spectra(MSIobject, "conc - pg/pixel") = matrix(nrow = nrow(MSIobject), ncol = ncol(MSIobject))
            spectra(MSIobject, "conc - pg/mm2") = matrix(nrow = nrow(MSIobject), ncol = ncol(MSIobject))

            no_cal_indices = c()

            for(i in 1:nrow(fData(MSIobject))){

              eqn = cal_list[[i]]

              if(typeof(eqn) != "list"){

                print(sprintf("No calibration curve for feature %s", fData(MSIobject)$name[i]))
                no_cal_indices = c(no_cal_indices, i)

              } else{

                pixel_ints = spectraData(MSIobject)[[val_slot]][i,]
                pixel_concs = as.vector(sapply(pixel_ints, function(y) inverse.predict(eqn, y)$Prediction))

                mm2_scalar = (1000 / as.numeric(experimentData(MSIobject)$pixelSize)) ^2
                mm2_conc = pixel_concs * mm2_scalar

                spectra(MSIobject, "conc - pg/pixel")[i,] = pixel_concs
                spectra(MSIobject, "conc - pg/mm2")[i,] = mm2_conc

              }
            }

            if(length(no_cal_indices) > 0){
              MSIobject = MSIobject[-no_cal_indices, ]
            }

            return(MSIobject)
          })
