# correlation colors for corplot, red is high, blue is low.
col.pairs <- circlize::colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red'), 0.5)

#' Plots the correlation among the columns of a numeric matrix.
#'
#' We assume that this is a sample x gene expression matrix, but it can
#' (of course) be any numeric matrix of your choosing. The column names appear
#' in the main diagonal of the plot. Note that you might prefer the corrplot
#' package for similar functionality, and this functionality is intentionally
#' named different from that..
#'
#' TODO: Add with.signature parameter to allow a box to plot the signature
#' score of all genes in E.
#'
#' @export
#' @importFrom graphics smoothScatter text
#' @seealso The [corrplot](http://cran.r-project.org/package=corrplot) package
#'
#' @param E the matrix used to plot a pairs correlation plot. The vectors used
#'   to assess all pairwise correlation should be *in the columns* of the
#'   matrix.
#' @param title The title of the plot
#' @param cluster `logical` indicating whether or not to shuffle genes
#'   around into some clustering.
#' @param col.point the color of the points in the scatterplots
#' @param diag.distro show the distribution of values on the diagnols?
#' @param smooth.scatter boolean to indicate wether to use a normal scatter, or
#'   a [graphics::smoothScatter()]. Defaults to `TRUE` if `nrow(E) > 400`
#' @param ... pass through arguments to internal panel functions
#' @param max.cex.cor the numeric value defining the maximum text size (cor) in the correlation panel.
#'   By default there is no limit on the maximum text size and the text size is calculated with `0.8 / strwidth(text)`.
#'   With `max.cex.cor` defined the text size is calculated as `min(0.8 / strwidth(text), max.cex.cor)`.
#' @return nothing, just creates the plot
#'
#' @examples
#' x <- matrix(rnorm(1000), ncol=5)
#' corplot(x)
corplot <- function(E, title, cluster = FALSE, col.point = '#00000066',
                    diag.distro = TRUE, smooth.scatter = nrow(E) > 400, 
                    max.cex.cor = NULL, ...) {
  E <- as_matrix(E, ...)
  if (missing(title)) {
    title <- 'Pairs Plot'
  }

  if (cluster) {
    cors <- cor(E, method='spearman', use='na.or.complete')
    dists <- as.dist((1 - cors) / 2)
    hc <- hclust(dists)
    dendro <- as.dendrogram(hc)
    idxs <- rev(order.dendrogram(dendro))
    E <- E[, idxs]
  }

  # create a copy of panel.spoints function
  # and update the default values with the ones provided in `corplot`
  # workaround to avoid running `pairs` with `col.point` and `smooth.scatter` params
  # as this resulted in multiple warnings
  c.panel.spoints <- panel.spoints
  formals(c.panel.spoints)$col.point <- col.point
  formals(c.panel.spoints)$smooth.scatter <- smooth.scatter

  # same for the panel.cor
  c.panel.cor <- panel.cor
  if (!is.null(max.cex.cor)) { 
  formals(c.panel.cor)$max.cex.cor <- max.cex.cor
  }
  
  
    pairs(E, ..., 
          lower.panel=c.panel.cor,
          upper.panel=c.panel.spoints,
          diag.panel=if (diag.distro) panel.hist else NULL,
          gap=0.2, pch=16,
          main=title)

  invisible(E)
}

## Helper functions ------------------------------------------------------------

#' Helper function that fills in the lower left part of the corrplot.
#'
#' The correlation between the two vectors is calculated and written in the
#' corresponding panel. The background of the panel is colored in accordance
#' to the strength and direction of the correlation.
#'
#' @noRd
panel.cor <- function(x, y, cex.cor, digits = 2, prefix = "", max.cex.cor = NULL, ...) {
  usr <- par(c("usr"))
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method='spearman', use='na.or.complete')
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if (missing(cex.cor)) {
    cex.cor <- if (is.null(max.cex.cor)) {
      0.8 / strwidth(txt)
    } else {
      min(0.8 / strwidth(txt), max.cex.cor)
    }
  }

  if (!is.na(r)) {
    bg.col <- col.pairs(r)
  } else {
    bg.col <- "#3f3f3f33"
  }

  rect(0, 0, 1, 1, col=bg.col)
  text(0.5, 0.5, txt, cex = cex.cor)
}

#' Helper function to draw the histogram on the diagonal of the corplot.
#'
#' Drawn when `corrplot(..., diag.distro = TRUE)`
#'
#' @noRd
panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr=c(usr[1:2], 0, 1.5))
  h <- hist(x, breaks=20, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  suppressWarnings(rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...))
}

#' Helper function: draws scatter plot in the upperright diagonal of corrplot.
#'
#' @noRd
panel.spoints <- function(x, y, col=par("col"), bg=NA, pch=par("pch"), cex=1,
                          col.smooth="red", span=2/3, iter=3,
                          col.point="#00000022", smooth.scatter=FALSE, ...) {
  ## points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  if (smooth.scatter) {
    smoothScatter(x, y, nrpoints = 0, add=TRUE)
  } else {
    points(x, y, pch=16, col=col.point, bg=bg, cex=cex)
  }

  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) {
    suppressWarnings({
      lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
            col = col.smooth, ...)
    })
  }
}
