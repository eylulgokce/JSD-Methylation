% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimators.R
\name{js_divergence}
\alias{js_divergence}
\title{Jensen-Shannon Divergence}
\usage{
js_divergence(indata, weights = NULL)
}
\arguments{
\item{indata}{A data frame with hierarchical column structure (sampling_unit, feature)}

\item{weights}{Optional weights for computing weighted averages}
}
\value{
A data frame with JSD values and associated statistics
}
\description{
Computes Jensen-Shannon divergence for genomic data with quality control filters.
}
