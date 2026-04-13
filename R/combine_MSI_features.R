require(Cardinal)

#' Function to create data matrix from MSI object
#' @import Cardinal
#' @import dplyr
#' @include read_mrm.R
#'
#' @param names Vector of file names, of Waters raw imaging files (including the imaging/Analyte 1.txt file after processing) - feature list in order of vector
#' @param folder Location of folder containing all the raw imaging files
#' @param lib_ion_path Full path to library of transitions for targeted MSI experiments. Must include headers: 'transition_id',	'precursor_mz',	'product_mz',	'collision_eV',	'cone_V',	'Polarity',	'Type'. Where 'Type' is "Analyte" for analytes in tissue and "IS" for internal standards.
#' @param run_name Change the run name of MSI object so can merge different panels
#' @param output_fn Folder and file name to save RDS object to (in the same folder as raw data)
#' @param save logical whether to save RDS for future use
#'
#' @return MSIobject with different features merged if imaging same area
#'
#' @export read_mrm
combine_MSI_features = function(names, folder, lib_ion_path, run_name = "combined", output_fn = "cobined", save=T){

  msi_files = list()

  for(fn in names){
    print(fn)
    msi_file = read_mrm(name = fn, folder = data_path, lib_ion_path = lib_ion_path)

    msi_files[[fn]] = msi_file
  }

  # Set count for mz channel
  count = 0
  for(list_ind in 1:length(msi_files) ){

    # Get indexed file from list of objects
    temp_msi = msi_files[[list_ind]]

    # Change run name for combining
    run(temp_msi) = run_name

    # Update mz channel so they are ordered and not duplicated after merging
    mz(temp_msi) = (count +1) : (count + nrow(temp_msi))
    # Update count
    count = max(mz(temp_msi))

    # Update msi list for combining
    msi_files[[list_ind]] = temp_msi

  }

  combined <- Reduce(rbind, msi_files)

  if(save ==T){
    # Create new folder
    output_folder = sprintf("%s/%s", folder, output_fn)
    dir.create(output_folder)

    # Save the combined DESI-MRM object as RDS - enables multiple spectra data channels to be stored in single object
    rdsfile <- tempfile(fileext=".RDS")
    saveRDS(combined, file = sprintf("%s/%s.RDS", output_folder, output_fn))
  }

  return(combined)
}

