% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SparrowResult-methods.R
\name{combine,SparrowResult,SparrowResult-method}
\alias{combine,SparrowResult,SparrowResult-method}
\title{Combines two SparrowResult objects together.}
\usage{
\S4method{combine}{SparrowResult,SparrowResult}(x, y, rename.x = NULL, rename.y = NULL, ...)
}
\arguments{
\item{x}{A \code{SparrowResult} object}

\item{y}{A \code{SparrowResult} object}

\item{rename.x}{A named vector that used to match resultNames(x) and remane
them to something different. \code{names(rename.x)} should match whatever you
want to change in \code{resultNames(x)}, and the values are the new names of
the result.}

\item{rename.y}{Same as \code{rename.x}, but for the results in \code{y}.}

\item{...}{more things}
}
\value{
A combined \code{SparrowResult} object
}
\description{
This would be useful when you want to add a GSEA result to an already
existing one. \code{append} would be more appropriate, but ...
}
\details{
When would you want to do that? Imagine a shiny app that drives sparrow.
You might want to present the results of each analysis as they come "online",
so you would run them independently and make them available to the user
immediately after they each finish (ie. in combination with the promises
package).
}
\examples{
mg1 <- exampleSparrowResult()
mg2 <- exampleSparrowResult()
mgc <- combine(mg1, mg2)
}
