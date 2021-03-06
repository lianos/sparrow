% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{SparrowResult-class}
\alias{SparrowResult-class}
\alias{.SparrowResult}
\alias{SparrowResult}
\title{A SparrowResult object holds the results from a sparrow::seas() call.}
\description{
A call to \code{\link[=seas]{seas()}} will produce analyses for an
arbitrary number of GSEA methods, the results of which will be stored and
accessible here using the \code{\link[=result]{result()}}, \code{\link[=results]{results()}}, and \code{\link[=resultNames]{resultNames()}}.

In addition, the \code{\link[=GeneSetDb]{GeneSetDb()}} used for the analysis is accessible
via \link[=GeneSetDb]{geneSetDb()}, and the results from the differential
expression analysis is available via \code{\link[=logFC]{logFC()}}.

Visualizing results of a geneset based analysis also are functions that
operate over a \code{SparrowResult} object, for instance see the
\code{\link[=iplot]{iplot()}} and the \code{sparrow.shiny} package.
}
\section{Slots}{

\describe{
\item{\code{gsd}}{The \code{\link[=GeneSetDb]{GeneSetDb()}} this analysis was run over}

\item{\code{results}}{The list of individual results generated by each of the
GSEA methods that were run.}

\item{\code{logFC}}{The differential expression statistics for each individual
feature measured in the experiment.}
}}

