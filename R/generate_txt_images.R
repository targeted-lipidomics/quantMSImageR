#' Load, process and optionally export per-feature text-image matrices
#'
#' Reads one or more DESI-MRM acquisitions, applies SNR filtering and
#' tissue/background separation, and writes per-feature text matrices that
#' can be imported by imaging software (e.g. MassLynx QuanOptimise overlay
#' tools). Setting `output_txt = FALSE` runs all processing but skips file
#' writing, which is useful when the caller only needs the processed MSI
#' objects (e.g. for HTML report rendering).
#'
#' @import Cardinal
#' @import pracma
#' @include setClasses.R
#'
#' @param fns Character vector of acquisition names (without `.raw` suffix)
#'   for single-polarity studies. For dual-polarity studies pass a list where
#'   each element is either a plain string (single acquisition) or a named list
#'   with fields `pos`, `neg`, `label`, and optionally `prefix_pos`/`prefix_neg`
#'   (see [bind_polarities()]). `run_study.R` builds this list automatically
#'   from the YAML `samples` section.
#' @param data_path Path to the folder that contains the `.raw` acquisition
#'   directories.
#' @param image_dir Root output directory. One sub-directory per acquisition
#'   is created inside it.
#' @param lib_ion_path Full path to the MRM ion-library CSV. Must contain
#'   columns: `transition_id`, `precursor_mz`, `product_mz`, `collision_eV`,
#'   `cone_V`, `Polarity`, `Type`.
#' @param snr_thresh Numeric. Minimum signal-to-noise ratio; pixels below this
#'   threshold are set to `NA` (default `3`).
#' @param tiss_fc Numeric. SNR threshold used for the tissue fold-change layer
#'   (default `0.6`).
#' @param thresh Numeric. Cold-spot percentile passed to `imageR()` as
#'   `threshold` (default `20`).
#' @param perc Numeric. Hot-spot percentile passed to `imageR()` as
#'   `percentile` (default `97`).
#' @param rot_clockwise Integer (0–3). Number of 90° clockwise rotations
#'   applied via `pracma::rot90()` (default `0`).
#' @param average_method Character, `"mean"` or `"median"`. Statistic used to
#'   summarise the noise pixel vector in `int2snr()`. Applied consistently to
#'   every acquisition (default `"median"`).
#' @param output_txt Logical. When `FALSE`, processing runs but no files are
#'   written to disk (default `TRUE`).
#' @param exclude Character vector of feature names to drop before processing.
#'   Names must match `fData()$name` exactly (default `NULL` — keep all).
#' @param rename Named character vector or list mapping old feature names to new
#'   display names, e.g. `c("old name" = "new name")`. Applied to all combined
#'   objects after loading; the new names appear in file names and the heatmap
#'   (default `NULL` — no renaming).
#'
#' @return Invisibly returns a named list with four processed
#'   `quant_MSImagingExperiment` objects:
#'   \describe{
#'     \item{`combined`}{Raw combined object (background pixels retained).}
#'     \item{`combined_snr`}{SNR-filtered intensity; sub-threshold pixels `NA`.}
#'     \item{`combined_FC`}{Tissue fold-change layer (`snr` slot).}
#'     \item{`combined_NAbackground`}{Raw intensity with background set to `NA`.}
#'   }
#'
#' @seealso [int2snr()], [applySNR()], [back2NA()], [imageR()]
#' @export
generate_txt_images <- function(
  fns,
  data_path,
  image_dir,
  lib_ion_path,
  snr_thresh     = 3,
  tiss_fc        = 0.6,
  thresh         = 20,
  perc           = 97,
  rot_clockwise  = 0,
  average_method = "median",
  output_txt     = TRUE,
  exclude        = NULL,
  rename         = NULL
) {

  average_method <- match.arg(average_method, c("mean", "median"))

  # ----- Internal helpers ------------------------------------------------

  make_txt_mat <- function(MSIobject, feat_ind, val_slot, value_label,
                            threshold, percentile) {
    # When every pixel is NA (analyte absent in this sample), return a
    # correctly-dimensioned zero matrix so downstream tools (spatialData /
    # tissuUmaps) receive a file with the right pixel grid.
    ivals <- as.numeric(spectraData(MSIobject[feat_ind, ])[[val_slot]])
    if (all(is.na(ivals))) {
      nx <- length(unique(coord(MSIobject)$x))
      ny <- length(unique(coord(MSIobject)$y))
      return(pracma::rot90(matrix(0, nrow = ny, ncol = nx), k = rot_clockwise))
    }

    tryCatch({
      result <- imageR(
        MSIobject  = MSIobject, val_slot   = val_slot, value      = value_label,
        scale      = "suppress", threshold  = threshold, sample_lab = "run",
        pixels     = NA, percentile = percentile, overlay    = FALSE,
        feat_ind   = feat_ind, blank_back  = FALSE, text_image  = TRUE
      )
      pracma::rot90(as.matrix(result), k = rot_clockwise)
    }, error = function(e) {
      message("imageR error (", value_label, "): ", e$message)
      matrix(NA_real_, 0, 0)
    })
  }

  write_if_nonempty <- function(mat, path) {
    if (nrow(mat) > 0)
      write.table(mat, path, sep = "\t", row.names = FALSE, col.names = FALSE)
  }

  normalize_0_100 <- function(mat) {
    if (nrow(mat) == 0) return(mat)
    rng <- range(mat, na.rm = TRUE, finite = TRUE)
    if (!all(is.finite(rng)) || diff(rng) == 0)
      return(matrix(0, nrow(mat), ncol(mat)))
    (mat - rng[1]) / diff(rng) * 100
  }

  safe_feat_name <- function(x) {
    x <- gsub("\\|\\|", "_or_", x)
    x <- gsub(":",      "_",    x)
    x <- gsub("/",      ".",    x)
    x
  }

  apply_renames <- function(obj, rename_map) {
    if (is.null(rename_map) || length(rename_map) == 0) return(obj)
    old_nms <- names(rename_map)
    for (i in seq_along(old_nms)) {
      ind <- fData(obj)$name == old_nms[i]
      if (any(ind)) fData(obj)$name[ind] <- rename_map[[i]]
    }
    featureNames(obj) <- fData(obj)$name
    obj
  }

  `%||%` <- function(a, b) if (is.null(a)) b else a

  # Load a single acquisition, attach tissue/noise labels, trim, exclude
  load_and_prep_acq <- function(fn_name) {
    obj  <- read_mrm(name = fn_name, folder = data_path, lib_ion_path = lib_ion_path)
    tpdf <- read.csv(sprintf("%s/%s.raw/tissue_pixels.csv", data_path, fn_name))
    pData(obj)$sample_name <- makeFactor(
      tissue_pixels = tpdf[["tissue_pixels"]],
      noise_pixels  = tpdf[["noise_pixels"]]
    )
    obj <- as(obj, "quant_MSImagingExperiment")
    obj <- trim_MSI(MSI_data = obj)
    if (!is.null(exclude) && length(exclude) > 0) {
      keep <- !fData(obj)$name %in% exclude
      if (!all(keep)) {
        message(sprintf("  Excluding %d transition(s): %s",
                        sum(!keep),
                        paste(fData(obj)$name[!keep], collapse = ", ")))
        obj <- obj[keep, ]
      }
    }
    obj
  }

  # Load one or more acquisitions of the same polarity, combine, set run label
  load_and_prep_multiple <- function(fns_vec, label) {
    fns_vec <- as.character(unlist(fns_vec))
    objs    <- lapply(fns_vec, load_and_prep_acq)
    obj     <- objs[[1]]
    for (i in seq_along(objs)[-1])
      obj <- combine_MSIs(obj, objs[[i]])
    pData(obj)$run <- factor(rep(label, ncol(obj)))
    obj
  }

  # ----- Normalise fns → fn_list / fn_labels ----------------------------
  # fns: character vector (backward-compat) OR list where each element is
  # a string (single acq) or named list with pos/neg/label fields.
  fn_list   <- as.list(fns)
  fn_labels <- vapply(fn_list, function(e)
    if (is.list(e)) e$label %||% (e$pos[1] %||% e$neg[1]) else as.character(e),
    character(1))

  # ----- Load and process acquisitions -----------------------------------

  combined     <- NULL
  combined_snr <- NULL
  combined_FC  <- NULL

  for (ind in seq_along(fn_list)) {
    fn_entry <- fn_list[[ind]]
    fn_label <- fn_labels[ind]
    message(sprintf("Loading %s (%d/%d)", fn_label, ind, length(fn_list)))

    if (is.list(fn_entry)) {
      pos_fns <- if (!is.null(fn_entry$pos)) unlist(fn_entry$pos) else NULL
      neg_fns <- if (!is.null(fn_entry$neg)) unlist(fn_entry$neg) else NULL

      if (!is.null(pos_fns) && !is.null(neg_fns)) {
        # Both polarities: combine within each polarity then bind across
        pos_obj <- load_and_prep_multiple(pos_fns, label = fn_label)
        neg_obj <- load_and_prep_multiple(neg_fns, label = fn_label)
        tissue  <- bind_polarities(pos_obj, neg_obj, label = fn_label)
      } else {
        # Single polarity (pos: OR neg: only)
        tissue <- load_and_prep_multiple(pos_fns %||% neg_fns, label = fn_label)
      }
    } else {
      # Plain string: backward-compatible single acquisition
      tissue <- load_and_prep_acq(fn_entry)
    }

    tissue_fc <- int2snr(
      MSIobject = tissue, val_slot = "intensity", sample_type = "sample_name",
      noise = "tissue_pixels", tissue = "tissue_pixels",
      snr_thresh = tiss_fc, average = average_method
    )
    tissue_snr <- int2snr(
      MSIobject = tissue, val_slot = "intensity", sample_type = "sample_name",
      noise = "noise_pixels", tissue = "tissue_pixels",
      snr_thresh = snr_thresh, average = average_method
    )
    tissue_snr <- applySNR(MSIobject = tissue_snr, val_slot = "intensity")

    if (is.null(combined)) {
      combined     <- tissue
      combined_snr <- tissue_snr
      combined_FC  <- tissue_fc
    } else {
      combined     <- combine_MSIs(combined,     tissue)
      combined_snr <- combine_MSIs(combined_snr, tissue_snr)
      combined_FC  <- combine_MSIs(combined_FC,  tissue_fc)
    }
  }

  combined_NAbackground <- back2NA(
    combined, val_slot = "intensity",
    background = "noise_pixels", tissue = "tissue_pixels",
    sample_type = "sample_name"
  )

  # Apply display-name overrides to all four objects
  if (!is.null(rename) && length(rename) > 0) {
    combined              <- apply_renames(combined,              rename)
    combined_snr          <- apply_renames(combined_snr,          rename)
    combined_FC           <- apply_renames(combined_FC,           rename)
    combined_NAbackground <- apply_renames(combined_NAbackground, rename)
  }

  out <- list(
    combined              = combined,
    combined_snr          = combined_snr,
    combined_FC           = combined_FC,
    combined_NAbackground = combined_NAbackground
  )

  if (!output_txt) return(invisible(out))

  # ----- Generate text-image files ---------------------------------------

  for (fn_label in fn_labels) {

    combined_snr_tmp  <- combined_snr[,         pData(combined_snr)$run          == fn_label]
    combined_FC_tmp   <- combined_FC[,           pData(combined_FC)$run           == fn_label]
    combined_back_tmp <- combined_NAbackground[, pData(combined_NAbackground)$run == fn_label]

    image_path <- file.path(image_dir, fn_label)
    dirs <- list(
      snr_filt  = file.path(image_path, sprintf("intensity_SNRfiltered%s",                       snr_thresh)),
      snr_norm  = file.path(image_path, sprintf("response_SNRfiltered%s_NORM",                   snr_thresh)),
      tissue_fc = file.path(image_path, sprintf("tissue-FC%s",                                   tiss_fc)),
      raw       = file.path(image_path, "intensity_raw"),
      hs_filt   = file.path(image_path, sprintf("intensity_SNRfiltered%s_hs%s_cs%s_removal",     snr_thresh, perc, thresh)),
      hs_norm   = file.path(image_path, sprintf("response_SNRfiltered%s_hs%s_cs%s_NORM",         snr_thresh, perc, thresh)),
      combined  = file.path(image_path, sprintf("response_SNRfiltered%s_hs%s_cs%s_COMBINED",     snr_thresh, perc, thresh))
    )
    for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

    n_features <- nrow(fData(combined_snr))

    for (feat_ind in seq_len(n_features)) {

      feat_name <- safe_feat_name(fData(combined_snr)$name[feat_ind])

      # SNR-filtered raw intensities
      mat_snr <- make_txt_mat(combined_snr_tmp, feat_ind, "intensity",
                               "DESI-MRM response - S/N filtered", 0, 100)
      write_if_nonempty(mat_snr, file.path(dirs$snr_filt, paste0(feat_name, ".txt")))

      # 0-100 normalised SNR
      scaled_snr <- normalize_0_100(mat_snr)
      write_if_nonempty(scaled_snr, file.path(dirs$snr_norm, paste0(feat_name, ".txt")))

      # COMBINED: embed global max in [1,1] for cross-sample colour scaling
      if (nrow(scaled_snr) > 0) {
        global_max <- max(
          as.numeric(spectraData(combined_snr[feat_ind, ])[["intensity"]]),
          na.rm = TRUE
        )
        scaled_snr_combined       <- scaled_snr
        scaled_snr_combined[1, 1] <- global_max
        write_if_nonempty(scaled_snr_combined,
                          file.path(dirs$combined, paste0(feat_name, ".txt")))
      }

      # Ratio-to-tissue (FC layer)
      mat_fc <- make_txt_mat(combined_FC_tmp, feat_ind, "snr",
                              "Ratio to tissue", 0, 100)
      write_if_nonempty(mat_fc, file.path(dirs$tissue_fc, paste0(feat_name, ".txt")))

      # Raw DESI-MRM response (background → NA)
      mat_raw <- make_txt_mat(combined_back_tmp, feat_ind, "intensity",
                               "DESI-MRM response", 0, 100)
      write_if_nonempty(mat_raw, file.path(dirs$raw, paste0(feat_name, ".txt")))

      # SNR-filtered + hot/cold-spot removal
      mat_hs <- make_txt_mat(combined_snr_tmp, feat_ind, "intensity",
                              "DESI-MRM response - S/N filtered", thresh, perc)
      write_if_nonempty(mat_hs, file.path(dirs$hs_filt, paste0(feat_name, ".txt")))

      # 0-100 normalised hs/cs version
      scaled_hs <- normalize_0_100(mat_hs)
      write_if_nonempty(scaled_hs, file.path(dirs$hs_norm, paste0(feat_name, ".txt")))
    }
  }

  invisible(out)
}
