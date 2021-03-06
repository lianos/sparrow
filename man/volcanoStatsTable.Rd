% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/volcano_plot.R
\name{volcanoStatsTable}
\alias{volcanoStatsTable}
\title{Extracts x and y axis values from objects to create input for volcano plot}
\usage{
volcanoStatsTable(
  x,
  stats = "dge",
  xaxis = "logFC",
  yaxis = "pval",
  idx = "idx",
  xtfrm = identity,
  ytfrm = function(vals) -log10(vals)
)
}
\arguments{
\item{x}{A \code{SparrowResult} object, or a \code{data.frame}}

\item{stats}{One of \code{"dge"} or \code{resultNames(x)}}

\item{xaxis, yaxis}{the column of the the provided (or extracted)
\code{data.frame} to use for the xaxis and yaxis of the volcano}

\item{idx}{The column of the \code{data.frame} to use as the identifier
for the element in the row. You probably don't want to mess with this}

\item{xtfrm}{A function that transforms the \code{xaxis} column to an
appropriate scale for the x-axis. This is the \code{identity} function
by default, because most often the logFC is plotted as is.}

\item{ytfrm}{A function that transforms the \code{yaxis} column to an
appropriate scale for the y-axis. This is the \code{-log10(yval)} function
by default, because this is how we most often plot the y-axis.}
}
\value{
a \code{data.frame} with \code{.xv}, \code{.xy}, \code{.xvt} and
\code{.xvy} columns that represent the xvalues, yvalues, transformed
xvalues, and transformed yvalues, respectively
}
\description{
You can, in theory, create a volcano plot from a number of different parts
of a \code{\link[=SparrowResult]{SparrowResult()}} object. Most often you want to create a volcano
plot from the differential expressino results, but you could imagine
building a volcan plot where each point is a geneset. In this case, you
would extract the pvalues from the method you like in the
\code{\link[=SparrowResult]{SparrowResult()}} object using the \code{stats} parameter.
}
\details{
Like the \code{\link[=volcanoPlot]{volcanoPlot()}} function, this is mostly used by the
\emph{sparrow.shiny} package.
}
\examples{
mg <- exampleSparrowResult()
v.dge <- volcanoStatsTable(mg)
v.camera <- volcanoStatsTable(mg, 'camera')
}
