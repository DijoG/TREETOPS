% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TREETOPS.R
\name{get_LEVELz}
\alias{get_LEVELz}
\title{Getting height levels II}
\usage{
get_LEVELz(CHM_g, height_bin)
}
\arguments{
\item{CHM_g}{chm from las, SpatRaster}

\item{height_bin}{numeric height bin vector returned by 'get_HB()'}
}
\value{
stacked SpatRaster
}
\description{
Function for getting height levels (stacked SpatRaster), uses output of 'get_HB()'.
}
