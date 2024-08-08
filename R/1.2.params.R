#' Data Input
#'
#' Creates a data frame of processing parameters used in scRNA-Seq data processing by Seurat
#'
#' @param f Path to the name of a data file in the "analysis/" folder. Accepts .xlsx, .csv, and .txt files.
#' @param md.num Number of metadata columns in input data file.
#' @param qc.rep Should a QC report be generated for the selected dataset? (TRUE/FALSE)
#' @return List containing formatted input data and metadata for downstream analysis.
#' @examples
#'
#' # # Dataset Input
#' # ms.input("example.xlsx",4,FALSE)
#'
#' @export
ms.input <- function(
    f,
    md.num
    ) {
  if(tools::file_ext(f) == "xlsx"){
    d <- readxl::read_excel(f,sheet = 1)
    list.params <- list(
      "data" = d[,(md.num+1):ncol(d)],
      "meta" = d[,1:md.num],
      "qc" = qc.rep
      )
    }
  if(tools::file_ext(f) == "csv"){
    d <- read_csv(
      f,
      check.names = T,
      header = T
      )
    list.params <- list(
      "data" = d[,(md.num+1):ncol(d)],
      "meta" = d[,1:md.num],
      "qc" = qc.rep
    )
  }
  if(tools::file_ext(f) == "txt"){
    d <- read_table(
      f,
      check.names = T,
      header = T,
      sep = "\t"
      )
    list.params <- list(
      "data" = d[,(md.num+1):ncol(d)],
      "meta" = d[,1:md.num],
      "qc" = qc.rep
    )
  }



  list.params <- data.frame(
    # Universal columns
    data.frame(
      # Sample Number
      Sample.No = seq(1:length(
        basename(
          list.files(d.path)
          )
        )
      ),
      # Individual file names (uses CellRanger folder name by default)
      File.ID = basename(list.files(d.path)),
      # Data file paths (Location of CellRanger files: 'Data/' by default)
      Path = paste(d.path,list.files(d.path),sep = ""),
      # Path to feature files
      Path.feat = paste(d.path,list.files(d.path),"/filtered_feature_bc_matrix/features.tsv.gz",sep = "")
      ),

    # Dataset-specific columns
    study.md
    )
  }




