% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TREETOPS.R
\name{get_TREETOPS}
\alias{get_TREETOPS}
\title{Getting treetops}
\usage{
get_TREETOPS(CHM_g, min_H, level_increment = 0.2)
}
\arguments{
\item{CHM_g}{SpatRaster object, raster chm derived from las}

\item{min_H}{numeric, minimum (Z) value of chm, higher values can be treetops}

\item{level_increment}{numeric, level cutting increment (default = 0.2m)}
}
\value{
sf object
}
\description{
MAIN FUNCTION for obtaining treetops.
}
