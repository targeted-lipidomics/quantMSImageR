library(Cardinal)
library(dplyr)
library(chemCal)

setGeneric("create_cal_curve", function(MSIobject, ...) standardGeneric("create_cal_curve"))

#' Function to create calibration curves (response v concentration, where concentration is pg/pixel)
#'
#' @import Cardinal
#' @import dplyr
#' @import chemCal
#' @include setClasses.R
#'
#' @param response_matrix matrix of average pg/pixel of m/z (rows = m/z and cols = cal level)
#' @param cal_type string of approach to generate claibration curve - 'std_addition' if standards are on tissue and 'cal' if direct onto glass slide.
#' @param level Column header to find background label from 'MSIobject@calibrationInfo@cal_response_data'
#' @param background string referring to background level from "level" label in 'MSIobject@calibrationInfo@cal_response_data'.
#' @return MSIobject with slots updated for i) cal_list - List of linear models for each m/z (response v concentration, where concentration is pg/pixel) and ii) r2 values for each calibration iii) calibration metadata
#'
#' @export
setMethod("create_cal_curve", "quant_MSImagingExperiment",
          function(MSIobject, cal_type = "std_addition", level = "level", background = "background"){

            cal_data = MSIobject@calibrationInfo@cal_response_data

            features = unique(cal_data$lipid)

            # Set outputs
            cal_list = list()
            r2_df = tibble(feature= features, r2 = NA)

            for(i in 1:length(features)){

              feature = features[i]
              cal_subset = subset(cal_data, lipid == feature)

              if(cal_type == "std_addition"){

                eqn = lm(response_perpixel~pg_perpixel, data = cal_subset, na.action = na.exclude)

                # Calculate background conc
                background_conc = inverse.predict(eqn, 0)$Prediction

                # Update conc values (original conc - background)
                cal_subset$pg_perpixel = cal_subset$pg_perpixel - background_conc


                cal_subset = cal_subset[cal_subset[[level]] != background, ]

              }

              cal_subset = mutate(cal_subset, pg_perpixel = ifelse(pg_perpixel == 0, yes= 1e-9, no = pg_perpixel))

              # Update equation
              eqn = lm(response_perpixel~pg_perpixel, data = cal_subset, na.action = na.exclude, weights = (1/pg_perpixel))

              r2_df$r2[i] = summary(eqn)[["r.squared"]]
              cal_list[[feature]] = eqn
            }

            MSIobject@calibrationInfo@cal_list = cal_list
            MSIobject@calibrationInfo@r2_df = r2_df

            return(MSIobject)
          })
