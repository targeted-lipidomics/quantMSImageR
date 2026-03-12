#!/usr/bin/env Rscript
# =============================================================================
# run_study.R  —  quantMSImageR DESI-MRM study runner
# =============================================================================
#
# Usage (command line):
#   Rscript run_study.R path/to/config.yaml
#
# Usage (interactive):
#   CONFIG_FILE <- "path/to/config.yaml"
#   source(system.file("run_study.R", package = "quantMSImageR"))
#
# The YAML must follow the structure in inst/config_template.yaml.
# =============================================================================

library(quantMSImageR)
library(Cardinal)
library(yaml)

# ---------------------------------------------------------------------------
# Resolve config path
# ---------------------------------------------------------------------------
if (!exists("CONFIG_FILE")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0)
    stop("Provide a YAML config path:\n  Rscript run_study.R config.yaml")
  CONFIG_FILE <- args[1]
}

if (!file.exists(CONFIG_FILE))
  stop("Config file not found: ", CONFIG_FILE)

cfg <- yaml::read_yaml(CONFIG_FILE)

# ---------------------------------------------------------------------------
# Unpack config
# ---------------------------------------------------------------------------
data_path     <- cfg$paths$data_path
out_path      <- cfg$paths$out_path
image_dir     <- cfg$paths$image_dir
lib_ion_path  <- cfg$paths$lib_ion_path
markdown_rmd  <- cfg$paths$markdown_rmd

fns           <- vapply(cfg$samples, `[[`, character(1), "name")
heatmap_labs  <- vapply(cfg$samples, `[[`, character(1), "label")
heatmap_order <- fns

snr_thresh     <- cfg$parameters$snr_thresh     %||% 3
tiss_fc        <- cfg$parameters$tiss_fc        %||% 0.6
thresh         <- cfg$parameters$thresh         %||% 20
perc           <- cfg$parameters$perc           %||% 97
rot_clockwise  <- cfg$parameters$rot_clockwise  %||% 0
average_method <- cfg$parameters$average_method %||% "median"
baseline_label <- cfg$parameters$baseline_label %||% heatmap_labs[1]

render_report  <- cfg$output$render_report %||% TRUE
output_txt     <- cfg$output$output_txt    %||% TRUE
report_fn      <- cfg$output$report_fn     %||% paste0(cfg$study, "_report")

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Process acquisitions (always runs; controls output via flags)
# ---------------------------------------------------------------------------
result <- generate_txt_images(
  fns            = fns,
  data_path      = data_path,
  image_dir      = image_dir,
  lib_ion_path   = lib_ion_path,
  snr_thresh     = snr_thresh,
  tiss_fc        = tiss_fc,
  thresh         = thresh,
  perc           = perc,
  rot_clockwise  = rot_clockwise,
  average_method = average_method,
  output_txt     = output_txt
)

# ---------------------------------------------------------------------------
# Render HTML report
# ---------------------------------------------------------------------------
if (render_report) {

  if (is.null(markdown_rmd) || !nzchar(markdown_rmd))
    stop("render_report is TRUE but paths$markdown_rmd is not set in the YAML.")
  if (!file.exists(markdown_rmd))
    stop("Rmd not found: ", markdown_rmd)

  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  # Expose objects and metadata expected by the Rmd
  combined       <- result$combined_snr
  # heatmap_order, heatmap_labs, baseline_label already in scope

  rmarkdown::render(
    markdown_rmd,
    output_file = file.path(
      out_path,
      paste0(report_fn, "_response_SNRfiltered.html")
    )
  )

  message("Report written to: ",
          file.path(out_path, paste0(report_fn, "_response_SNRfiltered.html")))
}

message("Done.")
