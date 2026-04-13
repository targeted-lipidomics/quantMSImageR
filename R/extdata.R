#' Dataset tissue_MRM_data
#'
#' This is a dataset containing Waters MRM data of SphingoMyelin 16:00 in Guinea pig lung tissue.
#'
#' @name tissue_MRM_data
#'
#' @section tissue_MRM_data.raw/imaging/Analyte 1.txt:
#'
#' This data is used in read_mrm().
NULL


#' Dataset cal_MRM_data
#'
#' This is a dataset containing Waters MRM data of SphingoMyelin 16:00 in Guinea pig lung tissue.
#'
#' @name cal_MRM_data
#'
#' @section cal_MRM_data.raw/imaging/Analyte 1.txt:
#'
#' This data is used in read_mrm().
NULL


#' Dataset ion_library
#'
#' This is a table containing a library of MRM transitions included in targeted MSI experiments.
#'
#' @name ion_library
#'
#' @section ion_library.txt:
#'
#' This data is used in read_mrm().
NULL


#' Dataset cal_rois
#'
#' This is a table containing a logical data about caliobration levels pixels are associated with.
#'
#' @name cal_rois
#'
#' @section cal_MRM_data.raw/cal_ROIs.csv:
#'
#' This data is used to identify calibration spots.
NULL


#' Dataset calibration_metadata
#'
#' This is a table containing information about concentrations in the calibration spots.
#'
#' @name calibration_metadata
#'
#' @section calibration_metadata.csv:
#'
#' This data is used to create calibration curves in create_cal_curve().
NULL


#' Dataset tissue_pixels
#'
#' This is a table containing a logical data about whether pixel is part of a calibration level.
#'
#' @name tissue_pixels
#'
#' @section tissue_MRM_data.raw/tissue_ROIs.csv:
#'
#' This data is used to identify tissue pixels.
NULL


#' Dataset tissue_rois
#'
#' This is a table containing a logical data about which tissue types pixels are associated with.
#'
#' @name tissue_rois
#'
#' @section tissue_MRM_data.raw/anat_ROIs.csv:
#'
#' This data is used to identify and label tissue types.
NULL

#' Dataset calibration_metadata
#'
#' This is a table containing information about the study sample(s).
#'
#' @name sample_metadata
#'
#' @section sample_metadata.csv:
#'
#' This data is used to update the pData() of MSIobjects.
NULL

#' Dataset combined
#'
#' This is an DESI-MRM dataset (in imzML format) containing tissue and calibration data.
#'
#' @name combined
#'
#' @section tissue_MRM_data.raw/combined.imzML:
#'
#' This data is used to test the createMSIDatamatrix() summarise_cal_levels() and create_cal_curve() functions.
NULL

#' Dataset cal_response_data
#'
#' This is a table containing mean intensity per pixel for the ROI at each calibration level across all calibration replicates.
#'
#' @name cal_response_data
#'
#' @section cal_response_data.csv:
#'
#' This data is to test the create_cal_curve() function.
NULL

#' Dataset cal_curve_MSI
#'
#' This is an DESI-MRM dataset (in imzML format) containing tissue and calibration data as well as the linera models for quantification.
#'
#' @name cal_curve_MSI
#'
#' @section tissue_MRM_data.raw/cal_curve_MSI.RDS:
#'
#' This data is used to test int2conc() function.
NULL

#' Microscope image of H&E stain (.czi)
#'
#' Microscope image of H&E stain (.czi) to extract image coordinates
#'
#' @name Slide12
#'
#' @section Slide12.czi:
#'
#' This data is used to test colocalising images
NULL

#' Shrunk microscope image of H&E stain (.png)
#'
#' Shrunk image of H&E stain (.png) to colocalise with ion image
#'
#' @name Slide12_shrunk
#'
#' @section Slide12_shrunk.png:
#'
#' This data is used to test colocalising images
NULL

#' Ion image (.png)
#'
#' Ion image of 12,13_DiHOME (.png) to colocalise with H&E stain
#'
#' @name snr_image
#'
#' @section snr_image.png:
#'
#' This data is used to test colocalising images
NULL
