> **This repository has been archived and is no longer maintained.**
> Active development has moved to **[MJS-708/quantMSImageR](https://github.com/MJS-708/quantMSImageR)**.
> Please update any bookmarks, citations, or `remotes::install_github()` calls accordingly.

[![](https://badgen.net/static/Publication/10.1021.acs.analchem.4c02350/green?.svg)](https://doi.org/10.1021/acs.analchem.4c02350) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10807654.svg)](https://doi.org/10.5281/zenodo.10807654)

# quantMSImageR

Software tools for processing and quantifying targeted multiple reaction monitoring (MRM) mass spectrometry imaging (MSI) data acquired by DESI-MRM.

Built on [Cardinal](https://bioconductor.org/packages/Cardinal/).

## Publication

**Development of a Desorption Electrospray Ionization–Multiple-Reaction-Monitoring Mass Spectrometry (DESI-MRM) Workflow for Spatially Mapping Oxylipins in Pulmonary Tissue**

Matthew J. Smith, Mu Nie, Mikael Adner, Jesper Säfholm, Craig E. Wheelock

*Analytical Chemistry* (2024). DOI: [10.1021/acs.analchem.4c02350](https://doi.org/10.1021/acs.analchem.4c02350)

## Requirements

- R ≥ 4.4.1
- [Cardinal](https://bioconductor.org/packages/Cardinal/) ≥ 3.6.2 (Bioconductor)
- [ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap/) (Bioconductor)

Install Bioconductor dependencies first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Cardinal", "ComplexHeatmap"))
```

## Installation

Install the latest release from GitHub:

```r
# install.packages("remotes")
remotes::install_github("MJS-708/quantMSImageR", ref = "main")
```

For the development version:

```r
remotes::install_github("MJS-708/quantMSImageR", ref = "dev")
```

Local development install:

```r
devtools::install("path/to/quantMSImageR")
```

## Usage

Analysis is configured via a YAML file and run with a single command:

```r
CONFIG_FILE <- "path/to/config.yaml"
source(system.file("run_study.R", package = "quantMSImageR"))
```

This will:
1. Load and process all acquisitions defined in the YAML
2. Apply SNR filtering, tissue fold-change masking, and IS normalisation
3. Export per-feature ion image `.txt` files for each sample
4. Render an HTML report with heatmaps, comparison plots, and ion images

A template YAML config is provided:

```r
file.copy(system.file("config_template.yaml", package = "quantMSImageR"), "config.yaml")
```

### YAML structure

```yaml
study: "study_name"

paths:
  data_path:    "path/to/raw/data"
  out_path:     "path/to/output"
  image_dir:    "path/to/ion/images"
  lib_ion_path: "path/to/ion_library.csv"
  markdown_rmd: "path/to/quantMSImageR___general_heatmap.Rmd"

samples:
  - neg: "acquisition_filename"   # without .raw extension
    label: "GroupA"
  - neg: "acquisition_filename2"
    label: "GroupB"
  - neg: "acquisition_filename3"
    label: "GroupA"               # duplicate labels = biological replicates

parameters:
  snr_thresh:  3
  baseline_label: "GroupA"

output:
  render_report: true
  report_fn: "my_study_SNR3"
```

## Contact

For questions or comments please contact Matthew J. Smith (matthew.smith@ki.se).
