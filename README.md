# Stevens.MSAnalyze v1.0

Processing and Analysis of Metabolomics Datasets

## Description

Utilizes various R packages to perform processing and analysis of metabolomics and lipidomics datasets.    The methods included in this package provide a seamless workflow for standard data processing, statistical methods, and visualization of metabolomics data.    Most analyses can be run in parallel to expedite time-consuming analyses.    The package is compatible with all operating systems. However, parallel processing is only available on Mac and Linux OS.

## Getting Started

### Dependencies
* Windows 10-11, WSL Ubuntu 22.04 or higher, Linux Ubuntu 22.04 or higher, or macOS 12.7.1 or higher
* R version 4.3.1 or higher (https://cran.r-project.org/)
* (Optional) RStudio version 2023.06.2 or higher (https://posit.co/download/rstudio-desktop/)
* (Optional) Conda installation of Python 3.9 and umap-learn (Only used for implementing UMAP via Python)
* R-packages (downloaded from CRAN or Bioconductor):
    * Suggests: 
        * knitr,
        * rmarkdown,
        * reticulate,
        * BiocManager
    * Imports: 
        * dplyr,
        * ggplot2,
        * ggsci,
        * viridis,
        * readxl,
        * ggpubr,
        * igraph,
        * ggraph,
        * shadowtext,
        * parallel,
        * reshape2,
        * ggrepel,
        * plotly,
        * htmlwidgets,
        * umap,
        * mvnormtest,
        * grid,
        * EnhancedVolcano,
        * circlize,
        * ComplexHeatmap

### Installation
* Run the following in a new R session on the command line or within R-Studio:

```
devtools::install_github(
  "cschasestevens/Stevens.MSAnalyze", 
  ref = "master", 
  build_vignettes = T
  )
```

## Help
* Browse vignettes by running the following:

```
browseVignettes("Stevens.MSAnalyze")
```

* Access function documentation by running the following:

```
# Type function name after package name
?Stevens.MSAnalyze::ms_input
```

## Authors

* Nathanial Chase Stevens, PhD, University of North Carolina at Chapel Hill
* Email: Nathanial_Stevens@med.unc.edu
* Alternate email: cschasestevens@gmail.com
* LinkedIn: https://www.linkedin.com/in/nathanial-chase-stevens-phd-08775180/

## Version History
* 1.1 (In Progress)
    * Added functions for regression and correlation analysis
    * Additional plotting functions for violin plots
* 1.0
    * Initial Release

## License

This project is licensed under the GNU General Public License Version 3 - see the LICENSE.md file for details

## Acknowledgments

* ComplexHeatmap package: Gu Z, Eils R, Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics. <doi:10.1093/bioinformatics/btw313>
* Circlize package: Gu, Z. circlize implements and enhances circular visualization in R. Bioinformatics 2014.
