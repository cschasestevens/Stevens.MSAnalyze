% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2.2.data.trf.sum.R
\name{ms_data_check}
\alias{ms_data_check}
\title{Data Transformation, Scaling and Quality Check}
\usage{
ms_data_check(ld, cl_var, samp_id, sc_d = TRUE, sc_meth = "Pareto")
}
\arguments{
\item{ld}{Data list generated from ms_input().}

\item{cl_var}{Variable for annotating sample correlation heatmap.}

\item{samp_id}{Variable containing sample IDs.}

\item{sc_d}{Transform and scale data?}

\item{sc_meth}{Data scaling method (performed after log2 transformation);
Currently available methods are "Median" or "Pareto" (default).}
}
\value{
List containing formatted input data and metadata
for downstream analysis.
}
\description{
Performs data imputation, log2-transformation, scaling, and
quality checks on an untargeted dataset.
}
\examples{

# ms_data_check(d)

}
