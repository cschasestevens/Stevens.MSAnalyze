% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4.1.dimrd.R
\name{ms_dim_rd}
\alias{ms_dim_rd}
\title{Dimension Reduction}
\usage{
ms_dim_rd(
  mat1,
  md,
  md_var,
  dim_type,
  ret_umap = FALSE,
  dim1,
  p_type,
  p_lab = TRUE
)
}
\arguments{
\item{mat1}{Pareto-scaled data matrix returned by ms_data_check().}

\item{md}{Data frame containing sample information.}

\item{md_var}{Vector of either a single variable or multiple variables
to plot. For single dimension reduction plots, only a single variable
should be listed.}

\item{dim_type}{Dimension reduction technique to use for plotting.
Currently available options are "PCA", "UMAP", or "PLSDA".}

\item{ret_umap}{Only used if dim_type = "UMAP"; Use Python umap-learn
for calculating UMAP (TRUE or FALSE).}

\item{dim1}{Number of dimensions to plot data (Either "2D" or "3D").
Only 2D plots are used if plotting multiple panels.}

\item{p_type}{Plot either a single variable or a list of metadata
variables. Enter either "single" or "list" to control output.}

\item{p_lab}{Label groups on plot panels (Only available if "2D").}
}
\value{
A single plot or series of dimension reduction plots.
}
\description{
Performs dimension reduction for a selected dataset using one
of the following methods: PCA, UMAP, PLSDA.
}
\examples{

# ms_dim_rd(
#   mat1 = d[["data"]][["data.scale"]],
#   md = d[["data"]][["meta"]],
#   md_var = c("Organ"),
#   dim_type = "PCA",
#   ret_umap = TRUE,
#   dim1 = "3D",
#   p_type = "single",
#   p_lab = FALSE
# )

}
