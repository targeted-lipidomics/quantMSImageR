library(Cardinal)

setGeneric("int2response", function(MSIobject, ...) standardGeneric("int2response"))

#' Function to normalise the intensity values to response per pixel if internal standard is present. Currently only works for a single internal standard to normalise all lipids to.
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject MSI MSIobject from Cardinal
#' @param val_slot character defining slot name to normalise - takes "intensity" as default
#' @param IS_name character defining the IS to use in fData(MSIobject) under the analyte header. Currently only works for a single internal standard to normalise all lipids to. "None" use individual lipids to normalise.
#' @param mode Mode iby which to apply normalisation. "sample" = normalise to median intensity of IS per sample, "line" = normalise to median intensity of IS per line, "pixel" = normalise to intensity of IS per pixel.
#' @param remove_IS Logical whether to remove internal standard feature from object
#' @return MSIobject with intensity values replaced with response
#'
#' @export
setMethod("int2response", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "intensity", IS_name = "None", mode = "line", remove_IS = T, ...){

            if(IS_name == "None"){
              IS_ind = NULL
            } else if(!any(fData(MSIobject)$analyte == IS_name)){
              print("No IS in this study so normalise to individual lipids")

              IS_ind = NULL
            } else{
              IS_ind = which(fData(MSIobject)$analyte == IS_name)
            }

            spectra(MSIobject, "response") = matrix(nrow = nrow(MSIobject), ncol = ncol(MSIobject))

            # Iterate over samples in study
            for(sample in unique(pData(MSIobject)$run)){

              # Select IS mz and pixels for specific sample
              sample_pixels = which(pData(MSIobject)$run == sample)

              tempMSIobject = MSIobject[, sample_pixels]

              # Save IS intensity vector
              IS_vec = spectraData(tempMSIobject)[[val_slot]][IS_ind, ]

              for(mz_ind in 1:nrow(fData(tempMSIobject))){

                #print(sprintf("mz - %s", mz_ind))
                ints = spectraData(tempMSIobject)[[val_slot]][mz_ind, ]

                if(mode == "pixel"){
                  if(!is.null(IS_ind)){
                    response = ints / IS_vec
                  } else{
                    response = ints / ints
                  }
                }
                if(mode == "sample"){
                  if(!is.null(IS_ind)){
                    response = ints / median(IS_vec, na.rm=T)
                  } else{
                    response = ints / median(ints, na.rm=T)
                  }
                }
                if(mode == "line"){
                  response = c()
                  for(line in unique(pData(tempMSIobject)$y)){

                    line_pixels = which(pData(tempMSIobject)$y == line)

                    if(!is.null(IS_ind)){
                      response = c(response, (ints[line_pixels] / median(IS_vec[line_pixels], na.rm=T)))
                    } else{
                      response = c(response, (ints[line_pixels] / median(ints[line_pixels], na.rm=T)))
                    }
                  }
                }

                spectra(MSIobject, "response")[mz_ind, sample_pixels] = response

              }
            }

            # Remove IS m/z
            if(remove_IS == T){
              MSIobject = MSIobject[-IS_ind, ]
            }

            return(MSIobject)
          })
