# generate_txt_images() has been moved into the package.
# Install / re-install quantMSImageR and use library(quantMSImageR) instead.
#
# This file is kept only so that existing source() calls do not break
# immediately. It will be removed in a future version.

if (!requireNamespace("quantMSImageR", quietly = TRUE)) {
  stop(
    "Please install quantMSImageR: devtools::install()\n",
    "generate_txt_images() is now a package function."
  )
}

# Re-export the package function into the calling environment
generate_txt_images <- quantMSImageR::generate_txt_images
