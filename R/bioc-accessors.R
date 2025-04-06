#' Wrapper for easy fetching of rowData/fData, etc.
#'
#' @noRd
fdata <- function(x, as.df = FALSE, ...) {
  UseMethod("fdata", x)
}

#' @noRd
#' @export
fdata.default <- function(x, ...) NULL

`fdata<-` <- function(x, value) {
  UseMethod("fdata<-", x)
}

#' @noRd
pdata <- function(x, as.df = FALSE, ...) {
  UseMethod("pdata", x)
}

#' @noRd
#' @export
pdata.default <- function(x, ...) NULL

`pdata<-` <- function(x, value) {
  UseMethod("pdata<-", x)
}

# DGEList ======================================================================

#' @noRd
fdata.DGEList <- function(x, as.df = FALSE, ...) {
  out <- x$genes
  if (!is.data.frame(out)) {
    warning("No `genes` data.frame found in (DG)EList", immediate. = TRUE)
  }
  out
}

#' @noRd
`fdata<-.DGEList` <- function(x, value) {
  x$genes <- value
  x
}

#' @noRd
pdata.DGEList <- function(x, ...) {
  x$samples
}

#' @noRd
`pdata<-.DGEList` <- function(x, value) {
  x$samples <- value
  x
}

# EList ========================================================================

#' @noRd
fdata.EList <- fdata.DGEList
`fdata<-.EList` <- `fdata<-.DGEList`

#' @noRd
pdata.EList <- function(x, ...) {
  x$targets
}

#' @noRd
`pdata<-.EList` <- function(x, value) {
  x$targets <- value
  x
}

# eSet =========================================================================

#' @noRd
fdata.eSet <- function(x, ...) {
  ns <- tryCatch(loadNamespace("Biobase"), error = function(e) NULL)
  if (is.null(ns)) stop("Biobase required")
  ns$fData(x)
}

#' @noRd
`fdata<-.eSet` <- function(x, value) {
  ns <- tryCatch(loadNamespace("Biobase"), error = function(e) NULL)
  if (is.null(ns)) stop("Biobase required")
  x <- ns$`fData<-`(x, value)
  x
}

#' @noRd
pdata.eSet <- function(x, ...) {
  ns <- tryCatch(loadNamespace("Biobase"), error = function(e) NULL)
  if (is.null(ns)) stop("Biobase required")
  ns$pData(x)
}

#' @noRd
`pdata<-.eSet` <- function(x, value) {
  ns <- tryCatch(loadNamespace("Biobase"), error = function(e) NULL)
  if (is.null(ns)) stop("Biobase required")
  x <- ns$`pData<-`(x, value)
  x
}

# SummarizedExperiment =========================================================

#' @noRd
fdata.SummarizedExperiment <- function(x, as.df = FALSE, ...) {
  if (!requireNamespace("SummarizedExperiment")) {
    stop("SummarizedExperiment package required")
  }
  out <- SummarizedExperiment::rowData(x)
  if (as.df) out <- SummarizedExperiment::as.data.frame(out)
  out
}

#' @noRd
`fdata<-.SummarizedExperiment` <- function(x, value) {
  if (!requireNamespace("SummarizedExperiment")) {
    stop("SummarizedExperiment package required")
  }
  SummarizedExperiment::rowData(x) <- value
  x
}

#' @noRd
pdata.SummarizedExperiment <- function(x, as.df = FALSE, ...) {
  if (!requireNamespace("SummarizedExperiment")) {
    stop("SummarizedExperiment package required")
  }
  out <- SummarizedExperiment::colData(x)
  if (as.df) out <- SummarizedExperiment::as.data.frame(out)
  out
}

#' @noRd
`pdata<-.SummarizedExperiment` <- function(x, value) {
  if (!requireNamespace("SummarizedExperiment")) {
    stop("SummarizedExperiment package required")
  }
  SummarizedExperiment::colData(x) <- value
  x
}
