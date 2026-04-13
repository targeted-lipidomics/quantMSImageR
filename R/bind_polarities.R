#' Merge positive- and negative-mode acquisitions from the same tissue area
#'
#' Combines two `MSImagingExperiment` objects (one per ionisation polarity)
#' that were acquired over the same pixel grid into a single object whose
#' feature set spans both polarities. Feature names are taken as-is from the
#' ion library (no prefix is added). Pixel metadata (coordinates, run label,
#' tissue/noise classification) is taken from the positive-mode object.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param pos_obj A `MSImagingExperiment` from a positive-mode acquisition
#'   (e.g. via [read_mrm()]).
#' @param neg_obj A `MSImagingExperiment` from a negative-mode acquisition
#'   covering the same pixel grid as `pos_obj`.
#' @param label Character. Value written to `pData(result)$run`; use the
#'   sample display label (e.g. `"A1"`) so the merged object integrates
#'   cleanly with [combine_MSIs()] and [quantile_hm()]. Defaults to the
#'   run label already in `pos_obj`.
#'
#' @return A `quant_MSImagingExperiment` with features from both polarities
#'   and pixel metadata from `pos_obj`.
#'
#' @seealso [read_mrm()], [combine_MSIs()], [generate_txt_images()]
#' @export
bind_polarities <- function(pos_obj, neg_obj, label = NULL) {

  n_pix_pos <- ncol(pos_obj)
  n_pix_neg <- ncol(neg_obj)

  if (n_pix_pos != n_pix_neg)
    stop(sprintf(
      "bind_polarities: pixel counts differ (POS = %d, NEG = %d). ",
      n_pix_pos, n_pix_neg,
      "POS and NEG must cover the same pixel grid."
    ))

  coord_pos <- as.data.frame(coord(pos_obj))
  coord_neg <- as.data.frame(coord(neg_obj))
  if (!identical(coord_pos, coord_neg))
    stop("bind_polarities: pixel coordinates differ between POS and NEG. ",
         "Register acquisitions to the same grid before binding.")

  # ---- stack intensity matrices (features x pixels) ------------------------
  idata_pos <- as.matrix(spectraData(pos_obj)[["intensity"]])
  idata_neg <- as.matrix(spectraData(neg_obj)[["intensity"]])
  combined_idata <- rbind(idata_pos, idata_neg)

  # ---- build combined featureData (mz indices offset so they stay unique) --
  fd_pos <- fData(pos_obj)
  fd_neg <- fData(neg_obj)

  combined_fdata <- MassDataFrame(
    mz           = c(fd_pos$mz, fd_neg$mz + max(fd_pos$mz)),
    analyte      = c(as.character(fd_pos$analyte),      as.character(fd_neg$analyte)),
    precursor_mz = c(as.character(fd_pos$precursor_mz), as.character(fd_neg$precursor_mz)),
    product_mz   = c(as.character(fd_pos$product_mz),   as.character(fd_neg$product_mz)),
    name         = c(fd_pos$name,                        fd_neg$name)
  )

  # ---- pixel metadata from POS (same grid) ---------------------------------
  pdata <- pixelData(pos_obj)
  if (!is.null(label))
    pdata$run <- factor(rep(label, n_pix_pos))

  # ---- assemble and return -------------------------------------------------
  out <- MSImagingExperiment(
    spectraData = combined_idata,
    featureData = combined_fdata,
    pixelData   = pdata
  )
  featureNames(out) <- combined_fdata$name
  as(out, "quant_MSImagingExperiment")
}
