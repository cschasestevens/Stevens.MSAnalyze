#' Data Input
#'
#' Creates a data frame of processing parameters
#' used for data processing and analysis.
#'
#' @param f Name of a data file in the "data/" folder.
#' Accepts .xlsx, .csv, and .txt files.
#' @param md_num Number of metadata columns in input data file.
#' @param qc_rep Should a QC report be generated for the selected dataset?
#' (Either TRUE or FALSE)
#' @param norm1 Logical indicating if the data have been normalized.
#' (Either TRUE or FALSE)
#' @param ref1 Reference annotation list for assigning compound details and
#' class information (provided as a .txt file).
#' @return List containing formatted input data and metadata
#' for downstream analysis.
#' @examples
#'
#' # ms_input(
#' #   "example.txt",
#' #   4,
#' #   FALSE,
#' #   TRUE,
#' #   "reference.txt"
#' # )
#'
#' @export
ms_input <- function(
  f,
  md_num,
  qc_rep,
  norm1,
  ref1
) {
  if(tools::file_ext(f) == "xlsx") { # nolint
    d <- readxl::read_excel(
      paste(
        "data/",
        f,
        sep = ""
      ),
      sheet = 1
    )
    a1 <- data.frame(
      "ID" = seq.int(
        1,
        length(names(d[, (md_num + 1):ncol(d)])),
        1
      ),
      "name" = names(d[, (md_num + 1):ncol(d)]),
      "type" = ifelse(
        grepl(
          "iSTD|1_CE 22:1|1_Sphingosine d17:1",
          names(d[, (md_num + 1):ncol(d)])
        ),
        "iSTD",
        "an.comp"
      )
    )
    if(missing(ref1) == TRUE) { # nolint
      r1 <- read.table(
        "ref/0.masterlist.txt",
        header = TRUE,
        sep = "\t"
      )
    }
    if(missing(ref1) == FALSE) { # nolint
      r1 <- ref1
    }
    lpar <- list(
      "data" = d[, -c(1:md_num, a1[a1[["type"]] == "iSTD", "ID"] + md_num)],
      "meta" = d[, 1:md_num],
      "anno" = a1,
      "istd" = d[, a1[a1[["type"]] == "iSTD", "ID"] + md_num],
      "qc" = qc_rep,
      "nrm" = norm1,
      "anno.ref" = r1
    )
  }
  if(tools::file_ext(f) == "csv") { # nolint
    d <- read.csv(
      paste(
        "data/",
        f,
        sep = ""
      ),
      check.names = TRUE,
      header = TRUE
    )
    a1 <- data.frame(
      "ID" = seq.int(
        1,
        length(names(d[, (md_num + 1):ncol(d)])),
        1
      ),
      "name" = names(d[, (md_num + 1):ncol(d)]),
      "type" = ifelse(
        grepl(
          "iSTD|1_CE.22.1|1_Sphingosine.d17.1",
          names(d[, (md_num + 1):ncol(d)])
        ),
        "iSTD",
        "an.comp"
      )
    )
    if(missing(ref1) == TRUE) { # nolint
      r1 <- read.table(
        "ref/0.masterlist.txt",
        header = TRUE,
        sep = "\t"
      )
    }
    if(missing(ref1) == FALSE) { # nolint
      r1 <- ref1
    }
    lpar <- list(
      "data" = d[, -c(1:md_num, a1[a1[["type"]] == "iSTD", "ID"] + md_num)],
      "meta" = d[, 1:md_num],
      "anno" = a1,
      "istd" = d[, a1[a1[["type"]] == "iSTD", "ID"] + md_num],
      "qc" = qc_rep,
      "nrm" = norm1,
      "anno.ref" = r1
    )
  }
  if(tools::file_ext(f) == "txt") { # nolint
    d <- read.table(
      paste(
        "data/",
        f,
        sep = ""
      ),
      check.names = TRUE,
      header = TRUE,
      sep = "\t"
    )
    a1 <- data.frame(
      "ID" = seq.int(
        1,
        length(names(d[, (md_num + 1):ncol(d)])),
        1
      ),
      "name" = names(d[, (md_num + 1):ncol(d)]),
      "type" = ifelse(
        grepl(
          "iSTD|1_CE.22.1|1_Sphingosine.d17.1",
          names(d[, (md_num + 1):ncol(d)])
        ),
        "iSTD",
        "an.comp"
      )
    )
    if(missing(ref1) == TRUE) { # nolint
      r1 <- read.table(
        "ref/0.masterlist.txt",
        header = TRUE,
        sep = "\t"
      )
    }
    if(missing(ref1) == FALSE) { # nolint
      r1 <- ref1
    }
    lpar <- list(
      "data" = d[, -c(1:md_num, a1[a1[["type"]] == "iSTD", "ID"] + md_num)],
      "meta" = d[, 1:md_num],
      "anno" = a1,
      "istd" = d[, a1[a1[["type"]] == "iSTD", "ID"] + md_num],
      "qc" = qc_rep,
      "nrm" = norm1,
      "anno.ref" = r1
    )
  }
  return(lpar)
}
