% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dotplot_markers.R
\name{dotplot_markers}
\alias{dotplot_markers}
\title{Dot Plot of Markers}
\usage{
dotplot_markers(seurat_object, markers)
}
\arguments{
\item{seurat_object}{the object to retrieve the data from.}

\item{markers}{a genes to plot.}
}
\description{
This function plots the expression of the markers across clusters in the Seurat object. It's analogous to the \code{Seurat::DotPlot()} but it doesn't change the order of markers.
If the named list of markers provided, the function use the names as titles for plots.
}
\examples{
dotplot_markers()
}
\keyword{dotplot}
