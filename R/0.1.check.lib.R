#' Install Libraries
#'
#' Installs CRAN and Bioconductor dependencies
#'
#' @return Summary of installed packages
#' @examples
#'
#'  # ms_check_lib()
#'
#' @export
ms_check_lib <- function() {
  list_pkg_CRAN <- c( # nolint
    "dplyr",
    "ggplot2",
    "ggsci",
    "viridis",
    "readxl",
    "ggpubr",
    "igraph",
    "ggraph",
    "shadowtext",
    "reshape2",
    "ggrepel",
    "plotly",
    "htmlwidgets",
    "umap",
    "mvnormtest",
    "grid",
    "circlize",
  )

  list_pkg_BIOC <- c( # nolint
    "parallel", "ComplexHeatmap",
    "EnhancedVolcano"
  )

  if(length(list_pkg_CRAN) > 0) { # nolint
    print("Attempting to install/update the following packages...")
    print(
      paste(
        "CRAN:",
        paste(
          list_pkg_CRAN,
          collapse = ", "
        ),
        sep = " "
      )
    )
    lapply(
      list_pkg_CRAN,
      function(x) {
        tryCatch(
          {
            install.packages(x)
          },
          error = function(e) {
            paste(
              "Latest version of",
              x,
              "already installed... skipping to next package.",
              sep = " "
            )
          }
        )
      }
    )
  }

  if("BiocManager" %in% installed.packages()[,"Package"] == FALSE) { # nolint
    install.packages("BiocManager")
  }

  if(length(list_pkg_BIOC) > 0) { # nolint
    print("Attempting to install/update the following packages...")
    print(
      paste(
        "Bioconductor:",
        paste(
          list_pkg_BIOC,
          collapse = ", "
        ),
        sep = " "
      )
    )
    BiocManager::install(list_pkg_BIOC)
  }
}