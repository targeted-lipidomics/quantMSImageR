setGeneric("combine_MSIs", function(MSIobject, ...) standardGeneric("combine_MSIs"))

#' Function to combine MSImagingExperiment objects.
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject MSImagingExperiment object from Cardinal
#' @return quant_MSImagingExperiment object with intensity values replaced with response
#' @param ... additional MSImagingExperiment object to combine - must have matching fData() and same columns form pData().
#'
#' @export
setMethod("combine_MSIs", "MSImagingExperiment",
          function(MSIobject, ...){

            objects <- c(as.list(environment()), list(...))

            f_data = fData(MSIobject)

            for(ind in 2:length(objects)){

              MSIobject = as( cbind(MSIobject, objects[[ind]]), 'MSImagingExperiment')

            }

            MSIobject = as(MSIobject, "quant_MSImagingExperiment")


            fData(MSIobject) = f_data

            return(MSIobject)

          })
