#' Data Input
#'
#' Creates a data frame of processing parameters used in scRNA-Seq data processing by Seurat
#'
#' @param f Path to the name of a data file in the "analysis/" folder. Accepts .xlsx, .csv, and .txt files.
#' @param md.num Number of metadata columns in input data file.
#' @param qc.rep Should a QC report be generated for the selected dataset? (TRUE/FALSE)
#' @param norm1 Logical indicating if the data have been normalized.
#' @param ref1 Reference annotation list for assigning compound details/class information (provided as a .txt file).
#' @return List containing formatted input data and metadata for downstream analysis.
#' @examples
#'
#' # # Dataset Input
#' # ms.input("example.xlsx",4,FALSE,TRUE)
#'
#' @export
ms.input <- function(
    f,
    md.num,
    qc.rep,
    norm1,
    ref1
    ) {
  if(tools::file_ext(f) == "xlsx"){
    d <- readxl::read_excel(paste("data/","example.xlsx",sep = ""),sheet = 1)
    a1 <- data.frame(
      "ID" = seq.int(
        1,
        length(names(d[,(md.num+1):ncol(d)])),
        1),
      "name" = names(d[,(md.num+1):ncol(d)]),
      "type" = ifelse(
        grepl(
          "iSTD|1_CE 22:1|1_Sphingosine d17:1",
          names(d[,(md.num+1):ncol(d)])),
        "iSTD",
        "an.comp"
        )
      )
    list.params <- list(
      "data" = d[,-c(1:md.num,a1[a1[["type"]] == "iSTD","ID"] + md.num)],
      "meta" = d[,1:md.num],
      "anno" = a1,
      "istd" = d[,a1[a1[["type"]] == "iSTD","ID"] + md.num],
      "qc" = qc.rep,
      "nrm" = norm1
      )
    }
  if(tools::file_ext(f) == "csv"){
    d <- read_csv(
      paste("data/",f,sep = ""),
      check.names = T,
      header = T
      )
    a1 <- data.frame(
      "ID" = seq.int(
        1,
        length(names(d[,(md.num+1):ncol(d)])),
        1),
      "name" = names(d[,(md.num+1):ncol(d)]),
      "type" = ifelse(
        grepl(
          "iSTD|1_CE.22.1|1_Sphingosine.d17.1",
          names(d[,(md.num+1):ncol(d)])),
        "iSTD",
        "an.comp"
      )
    )
    list.params <- list(
      "data" = d[,-c(1:md.num,a1[a1[["type"]] == "iSTD","ID"] + md.num)],
      "meta" = d[,1:md.num],
      "anno" = a1,
      "istd" = d[,a1[a1[["type"]] == "iSTD","ID"] + md.num],
      "qc" = qc.rep,
      "nrm" = norm1
    )
  }
  if(tools::file_ext(f) == "txt"){
    d <- read_table(
      paste("data/",f,sep = ""),
      check.names = T,
      header = T,
      sep = "\t"
      )
    a1 <- data.frame(
      "ID" = seq.int(
        1,
        length(names(d[,(md.num+1):ncol(d)])),
        1),
      "name" = names(d[,(md.num+1):ncol(d)]),
      "type" = ifelse(
        grepl(
          "iSTD|1_CE.22.1|1_Sphingosine.d17.1",
          names(d[,(md.num+1):ncol(d)])),
        "iSTD",
        "an.comp"
      )
    )
    list.params <- list(
      "data" = d[,-c(1:md.num,a1[a1[["type"]] == "iSTD","ID"] + md.num)],
      "meta" = d[,1:md.num],
      "anno" = a1,
      "istd" = d[,a1[a1[["type"]] == "iSTD","ID"] + md.num],
      "qc" = qc.rep,
      "nrm" = norm1
    )
  }
  return(list.params)
}




