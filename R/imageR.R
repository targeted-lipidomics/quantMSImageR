setGeneric("imageR", function(MSIobject, ...) standardGeneric("imageR"))

#' Function to create ion images (using ggplot)
#'
#' @import Cardinal
#' @import dplyr
#' @import chemCal
#' @import viridis
#' @import ggplot2
#' @import ggthemes
#' @include setClasses.R
#'
#' @param value character stating what values pertain to
#' @param scale suppress, histogram, sqrt
#' @param percentile percentile to suppress colour scale (for suppress.)if scale == "suppress")
#' @param threshold percintile to remove from colour scale (i.e background noise)
#' @param sample_lab character header from pData(MSIobject) to label images (defaults to "sample_ID")
#' @param pixels character from pData(MSIobject)[[sample_lab]] to take pixels to plot (defaults to NA)
#' @param overlay whether to overlay features (not yet implemented!)
#' @param feat_ind Index of feature from fData(MSIobject) to image
#' @param perc_scale Logical to normalise scale to % (TRUE) or displar raw values (FALSE)
#' @return ggplot object
#'
#' @export
setMethod("imageR", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "intensity", value = "response %", scale = "suppress", threshold = 1,
                   sample_lab = "sample_ID", pixels = NA, percentile=99.0, overlay = F,
                   feat_ind = 1, perc_scale = F, blank_back = T, aspect_ratio=1, text_image = F,
                   roi_labels = NA){

            MSIobject = as(MSIobject[feat_ind, ], "quant_MSImagingExperiment")

            if(!is.na(pixels)){
              MSIobject = MSIobject[, which(pData(MSIobject)$sample_type == pixels)]
            }

            # Generate image data matrix
            image_df = tibble(x = pData(MSIobject)@listData[["x"]],
                              y = pData(MSIobject)@listData[["y"]],
                              response = as.numeric(spectraData(MSIobject)[[val_slot]]),
                              sample = pData(MSIobject)[[sample_lab]],
                              feature = fData(MSIobject)$name)


            if(scale == "sqrt"){
              image_df = mutate(image_df, response = sqrt(response))
            }
            if(scale == "suppress"){
              vals = image_df$response
              vals[which(vals < 0)] = 0
              max_val = quantile(vals, percentile /100, na.rm=TRUE)
              if (!is.na(max_val) && max_val > min(vals, na.rm=TRUE)){
                vals[vals > max_val] = max_val
              }

              image_df$response = vals
            }
            if(scale == "histogram"){

              # Cardinal implementation	https://rdrr.io/bioc/Cardinal/src/R/DIP.R

              vals = image_df$response

              unique_percentiles = unique(quantile(vals, seq(from=0, to=1, length.out=100), na.rm=TRUE))
              if (length(unique_percentiles) < 5){
                return("Range of intensity values too narrow for histogram normalization.") }

              vals_interval = cut(vals, unique_percentiles, include.lowest=TRUE)
              vals_new = as.numeric(vals_interval) / length(levels(vals_interval))

              scale = mean(vals, na.rm=TRUE) / mean(vals_new, na.rm=TRUE) # So mean value same as prior to scaling
              vals_new = scale * vals_new

              image_df$response = vals_new

            }

            if(!is.na(threshold)){
              # Remove x% of values - threshold
              vals = image_df$response
              min_val = quantile(vals, threshold /100, na.rm=TRUE)
              if (!is.na(min_val) && min_val < max(vals, na.rm=TRUE)){
                vals[vals < min_val] = 0
              }
              image_df$response = vals
            }

            # Set NA and negative values to 0 for plotting!
            image_df = image_df |>
              dplyr::mutate(response = ifelse(is.na(response), 0, response)) |>
              dplyr::mutate(response = ifelse(response < 0, 0, response))

            if(perc_scale == T){
              vals = image_df$response
              if(max(vals) > 0){
                perc_vals = 100 * (vals / max(vals))

                image_df$response = perc_vals

              }
            }

            if(text_image == T){

              p = image_df |>
                dplyr::select(x, y, response) |>
                dplyr::arrange(as.numeric(x)) |>
                tidyr::pivot_wider(names_from = x, values_from = response) |>
                dplyr::arrange(as.numeric(y)) |>
                tibble::column_to_rownames("y")

            } else if(blank_back == T){
              image_df$response[is.na(image_df$response) | image_df$response <= 0] <- NA

              p = ggplot(data=image_df, aes(x = x, y = -y, fill = response)) +
                geom_tile() +
                theme_minimal() +
                theme(
                  aspect.ratio = aspect_ratio,
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.line = element_blank(),
                  panel.grid = element_blank(),
                  plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
                ) +
                scale_fill_viridis(na.value = "white") +
                labs(fill = value) +
                facet_grid(sample ~ feature)


            } else{
              image_df[is.na(image_df)] <- 0

              if(overlay == T){

                p = ggplot()+
                  geom_tile(data=image_df,aes(x=x,y=-y,fill=response), alpha=0.25) +
                  theme_minimal() +
                  theme(aspect.ratio=aspect_ratio,
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.line = element_blank(),
                        panel.grid = element_blank(),
                        plot.title = element_text(hjust = 0.5, face="bold", size = 15)) +
                  scale_fill_viridis(na.value = "white") +
                  labs(fill=value) +
                  geom_point(data = roi_labels, aes(x = transformed_x, y = -transformed_y, col = Cell_type),
                             shape = 8, na.rm = TRUE, size = 1.5)

              } else{

                p = ggplot(data=image_df,aes(x=x,y=-y,fill=response))+
                  geom_tile() +
                  theme_minimal() +
                  theme(aspect.ratio=aspect_ratio,
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.line = element_blank(),
                        panel.grid = element_blank(),
                        plot.title = element_text(hjust = 0.5, face="bold", size = 15)) +
                  scale_fill_viridis(na.value = "white") +
                  labs(fill=value) +
                  facet_grid(sample~feature)

              }
            }

            return(p)

        })
