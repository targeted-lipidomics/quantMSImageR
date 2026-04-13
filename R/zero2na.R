library(Cardinal)

setGeneric("zero2na", function(MSIobject, ...) standardGeneric("zero2na"))

#' Function to convert intensity values from 0 to NA in MSI dataset.
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject MSI object from Cardinal
#' @return MSIobject with intensity values replaced with response
#'
#' @export
setMethod("zero2na", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "intensity"){

            MSIobject = as(MSIobject, "quant_MSImagingExperiment")

            for(i in 1:nrow(fData(MSIobject))){

              #print(i)

              ints = spectraData(MSIobject)[[val_slot]][i,]

              if(all(is.na(ints))){
                print(sprintf("all intensities are NA for m/z %s. Doing nothing.", i))
              } else if(sum(ints, na.rm=T) == 0){
                print(sprintf("all intensities are 0 for m/z %s. Making NA.", i))
                spectra(MSIobject)[i, ] = NA
              } else{
                ints[which(ints == 0)] = NA
                spectraData(MSIobject)[[val_slot]][i, ] <- ints

              }
            }

            return(MSIobject)

          })
