#' Visualize gene level behavior of genes within a geneset across a contrast.
#'
#' @description
#' It is informative to look at the individual log fold changes of the genes
#' within a gene set to explore the degree to which they (1) are coherent with
#' respect to each other; and (2) see how the compare to the background
#' distribution of log fold changes of the entire transcriptome.
#'
#' You can visualize this behavior via a `type = "density"` plot, or a
#' `type = "boxplot". It is also common to plot either the individual
#' log fold changes `value = "logFC"` or t-statistics `value = "t"`.
#'
#' @rdname iplot
#' @export
#' @importFrom ggplot2 theme_bw
#'
#' @param x A [SparrowResult()] object
#' @param name the name of the geneset to plot
#' @param value A string indicating the column name for the value of the
#'   gene-level metadata to plot. Default is `"logFC"`. Anoter often used choice
#'   might also be `"t"`, to plot t-statistics (if they're in the result). But
#'   this can be any numeric column found in the data.frame returned by
#'   `geneSet(x, y, j)`. If this is a named string (vector), then the value in
#'   `names(value)` will be used on the axis when plotted.
#' @param type plot the distributions as a `"density"` plot or `"boxplot"`.
#' @param tools the tools to display in the rbokeh plot
#' @param main A title to display. If not specified, the gene set name
#'   will be used, otherwise you can pass in a custom title, or `NULL`
#'   will disable the title altogether.
#' @param with.legend Draws a legend to map point color to meaning. There are
#'   three levels a point (gene level statistic) can be color as, "notsig",
#'   "psig", and "sig". "notsig" implies that the FDR >= 10%, "psig" means that
#'   FDR <= 10%, but the logFC is "unremarkable" (< 1), and "sig" means
#'   that both the FDR <= 10% and the logFC >= 1
#' @param collection If you have genesets with duplicate names in `x`
#'   (only possible with a `GeneSetDb` object), provide the name of the
#'   collection here to disambiguate (default: `NULL`).
#' @param result_name if not `NULL` (the default), it fetches statistics for the
#'   given GSEA result using `result(x, name = result_name, ...)` which can then
#'   be used for plotting.
#' @param shiny_source the name of this element that is used in shiny callbacks.
#'   Defaults to `"mggenes"`.
#' @param width,height the width and height of the output plotly plot
#' @param ggtheme a ggplot theme, like the thing returned from
#'  `ggplot2::theme_bw()`, for instance.
#' @param trim used to define the upper and lower quantiles to max out the
#'   individual gene statistics in the selected geneset.
#' @param ... pass through parameters to internal boxplot/density/gsea
#'   plotting functions
#' @return the ploty plot object
#' @examples
#' mgr <- exampleSparrowResult()
#' upreg <- "SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP"
#' dnreg <- "TURASHVILI_BREAST_LOBULAR_CARCINOMA_VS_DUCTAL_NORMAL_DN"
#' 
#' iplot(mgr, upreg, value = c("t-statistic" = "t"), type = "density")
#' iplot(mgr, upreg, value = c("t-statistic" = "t"), type = "density", .plot_static = TRUE)
#' iplot(mgr, upreg, value = c("log2FC" = "logFC"), type = "boxplot")
#' iplot(mgr, upreg, value = c("t-statistic" = "t"), type = "gsea")
#' iplot(mgr, dnreg, value = c("t-statistic" = "t"), type = "gsea")
iplot <- function(x, name, value = "logFC",
                  type = c("density", "gsea", "boxplot"),
                  tools = c('wheel_zoom', 'box_select', 'reset', 'save'),
                  main = NULL, with.legend = TRUE,
                  collection = NULL,
                  result_name = NULL,
                  stats.report = c("NES", "logFC", "pval", "padj", "n"),
                  shiny_source = 'mggenes', width = NULL, height = NULL,
                  ggtheme = ggplot2::theme_bw(), trim = 0.005, 
                  interactive = TRUE,
                  ...) {
  if (FALSE) {
    x <- xmg; y <- 'H'; j <- 'HALLMARK_E2F_TARGETS'; value <- 'logFC';
    main <- NULL; type <- 'boxplot'; with.legend <- TRUE
  }
  assert_class(x, "SparrowResult")
  assert_string(name)
  assert_string(collection, null.ok = TRUE)
  assert_string(main, null.ok = TRUE)
  assert_character(stats.report, null.ok = TRUE)
  if (missing(result_name) && type == "gsea") {
    result_name <- "fgsea"
  }
  type <- match.arg(type)
  lfc <- copy(logFC(x, as.dt = TRUE)[, group := "bg"])
  match.arg(value, colnames(lfc))
  stopifnot(is.numeric(lfc[[value]]))

  if (is.null(main)) {
    main <- name
    if (!is.null(collection)) main <- sprintf("%s (%s)", main, collection)
  }

  gset <- geneSet(x, name = name, collection = collection, as.dt = TRUE)
  gcoll <- gset$collection[1L]
  gname <- gset$name[1L]
  gstats <- NULL

  if (test_character(stats.report, min.len = 1L) && test_string(result_name) ) {
    if (!result_name %in% resultNames(x)) {
      warning("No GSEA method result for `'", result_name, "'`")
    } else {
      gstats <- result(x, result_name, as.dt = TRUE)[J(gcoll, gname)]
      if (is.na(gstats$N)) {
        # How could we not found geneset?
        warning(
          "geneset not found for stat retrieval: ",
          paste(gcoll, gname, collapse = ",")
        )
        gstats <- NULL
      } else {
        stats.report <- intersect(stats.report, colnames(gstats))
        gstats <- as.list(gstats)[stats.report]
      }
    }
  }
  
  if (type == "gsea") {
    ret <- iplot.gsea.plot(
      lfc, 
      gset,
      rank_by = value,
      title = main,
      gstats = gstats,
      spr = x,
      result_name = result_name,
      .plot_static = !interactive,
      ...
    )
    return(ret)
  }

  # Avoid notes in R CMD check from data.table NSE mojo
  val <- NULL
  lfc[, val := lfc[[value]]]
  gset[, val := gset[[value]]]
  dat <- local({
    gset[, group := "geneset"]
    kcols <- intersect(names(gset), names(lfc))

    out <- rbind(lfc[, kcols, with = FALSE], gset[, kcols, with = FALSE])
    out[, significant := {
      tmp <- ifelse(significant, 'sig', 'notsig')
      ifelse(padj < 0.10 & tmp == 'notsig', 'psig', tmp)
    }]
  })

  if (type == 'density') {
    out <- iplot.density.plotly(x, dat, value, main,
                                with.legend=with.legend, tools=tools,
                                shiny_source=shiny_source,
                                ggtheme=ggtheme, trim = trim,
                                result_name = result_name,
                                gstats = gstats,
                                .plot_static = !interactive,
                                
                                ...)
  } else if (type == 'boxplot') {
    out <- iplot.boxplot.plotly(x, dat, value, main,
                                with.legend=with.legend, tools=tools,
                                shiny_source=shiny_source,
                                width=width, height=height, ggtheme=ggtheme,
                                trim=trim, gstats = gstats,
                                .plot_static = !interactive,
                                ...)
  } else if (type == 'volcano') {
  }

  out
}

## plotly ======================================================================

#' Reimplementation of fgsea::plotEnrichment so we can add more interactive
#' bits, and also to highlight genes on the leading edge, and what not.
#'
#' Lots of code here is copied from fgsea
#'
#' @noRd
#' @examples
#' 
#' 
iplot.gsea.plot <- function(lfc, geneset, rank_by, title, spr, gseaParam = 1,
                            ticksSize = 0.2, result_name = "fgsea", 
                            # gstats is the geneset statistics for this gs
                            gstats = NULL, ...,
                            .plot_default = FALSE,
                            .plot_static = FALSE) {
  checkmate::expect_class(spr, "SparrowResult")
  if (!requireNamespace("fgsea")) stop("'fgsea' package required")
  gcoll <- geneset$collection[1]
  gname <- geneset$name[1]
  # Setup params so we can just copy and paste fgsea::plotEnrichment code
  pathway <- geneset[["feature_id"]]
  stats <- setNames(lfc[[rank_by]], lfc[["feature_id"]])

  if (.plot_default) {
    return(fgsea::plotEnrichment(pathway, stats, gseaParam, ticksSize))
  }
  
  stats.anno <- NULL  
  if (!is.null(gstats)) {
    stats.anno <- sapply(gstats, function(val) {
      if (is.character(val)) return(val)
      if (test_int(val)) return(sprintf("%d", as.integer(val)))
      if (is.numeric(val)) return(sprintf("%0.3f", val))
      if (is.logical(val)) return(tolower(as.character(val)))
      as.character(val)
    })
    # stats.anno <- sprintf("%s: %s\n", names(stats.anno), unname(stats.anno))
    stats.anno <- paste(
      names(stats.anno),
      unname(stats.anno),
      sep = ": ",
      collapse = "\n"
    )
  }
  
  # fgsea::plotEnrichment ......................................................
  rnk <- rank(-stats)
  ord <- order(rnk)

  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))

  pathway <- unname(as.vector(stats::na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)

  gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway,
                                  returnAllExtremes = TRUE)

  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops

  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))

  diff <- (max(tops) - min(bottoms)) / 8

  # Getting rid of NOTEs
  x <- y <- NULL

  # Creates a data.frame to use for the line segments drawn for each  # :custom
  # feature. This allows us to add useful hover information.          # :custom
  features <- data.frame(                                             # :custom
    x = pathway,                                                      # :custom
    y = -diff / 2,                                                    # :custom
    xend = pathway,                                                   # :custom
    yend = diff / 2,                                                  # :custom
    feature_id = names(statsAdj)[pathway])                            # :custom

  add.labels <- c("feature_id", "statistic", "statistic_adj")         # :custom
  xref <- match(features[["feature_id"]], geneset[["feature_id"]])    # :custom

  features[["statistic"]] <- geneset[[rank_by]][xref]                 # :custom
  features[["statistic_adj"]] <- unname(statsAdj)[pathway]            # :custom

  label <- intersect(c("symbol", "name"), colnames(geneset))          # :custom
  if (length(label)) {                                                # :custom
    features[["name"]] <- geneset[[label[1L]]][xref]                  # :custom
    add.labels <- c("name", add.labels)                               # :custom
  }
  add.xtra <- c("logFC", "t", "pval", "padj")
  for (add.me in intersect(add.xtra, colnames(geneset))) {
    features[[add.me]] <- geneset[[add.me]]
    add.labels <- c(add.labels, add.me)
  }

  stat.cols <- c("x", "y", "xend", "yend")                            # :custom
  for (cname in setdiff(colnames(features), stat.cols))   {           # :custom
    if (is.numeric(features[[cname]])) {                              # :custom
      features[[cname]] <- sprintf("%0.3f", features[[cname]])        # :custom
    } else {
      features[[cname]] <- as.character(features[[cname]])            # :custom
    }
  }

  features[["label"]] <- sapply(seq_len(nrow(features)), function(i) {# :custom
    f <- features[i, add.labels]                                      # :custom
    paste(names(f), ":", unname(f[1,]), collapse = "<br>")            # :custom
  })

  if (is.null(names(rank_by))) {                                      # :custom
    xlabel <- sprintf("rank\n(by: %s)", rank_by)                      # :custom
  } else {
    xlabel <- sprintf("rank\n(by: %s)", names(rank_by))               # :custom
  }
  xend <- yend <- NULL # Silence NSE note in R CMD check
  g <- ggplot2::ggplot(toPlot, ggplot2::aes(x=x, y=y)) +
    ggplot2::geom_point(color = "green", size = 0.1) +
    ggplot2::geom_hline(yintercept = max(tops), colour = "red",
                        linetype = "dashed") +
    ggplot2::geom_hline(yintercept = min(bottoms), colour = "red",
                        linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, colour = "black") +
    ggplot2::geom_line(color = "green") +
    theme_bw() +
    # ggplot2::geom_segment(
    #   data = data.frame(x = pathway, label = ),
    #   mapping = ggplot2::aes(x = x, y = -diff / 2, xend = x, yend = diff / 2),
    #   size = ticksSize) +
    suppressWarnings({
      ggplot2::geom_segment(                                          # :custom
        data = features,                                              # :custom
        mapping = ggplot2::aes(                                       # :custom
          x = x, y = y, xend = xend, yend = yend, text = label),      # :custom
        size = ticksSize)                                             # :custom
    }) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(
      x = xlabel,
      y = "enrichment score",
      title = title)
  
  if (!is.null(stats.anno)) {
    if (is.null(gstats$NES) || gstats$NES > 0) {
      position <- "topright"
      hjust <- 1
      vjust <- 1
    } else {
      position <- "bottomleft"
      vjust <- 0
      hjust <- 0
    }
    
    if (.plot_static) {
      geom <- "label"
      if  (position == "topright") {
        x <- 0.95
        y <- 0.95
      } else {
        x <- 0.05
        y <- 0.05
      }
    } else {
      # ggplot2 can't do geom_label yet?
      geom <- "text"
      if (position == "topright") {
        x <- tail(toPlot$x, 1) - 100
        y <- max(tops) - 0.1 # 0.75
      } else {
        y <- min(bottoms) + 0.1 # -0.1
        x <- (toPlot$x[2] + toPlot$x[1]) / 2 # 10
      }
    }

    g <- g +
      ggplot2::annotate(
        geom,
        x = I(x),
        y = I(y),
        label = stats.anno,
        hjust = hjust,
        vjust = vjust
      )
  }

  if (!.plot_static) {
    g <- plotly::ggplotly(g, tooltip = "label")
    g <- layout(g, dragmode = "select")
    g <- config(g, displaylogo = FALSE)
  }

  g
}

#' @noRd
#' @importFrom plotly add_markers add_lines config layout plot_ly
iplot.density.plotly <- function(x, dat, value, main, with.legend=TRUE,
                                 with.points=TRUE, shiny_source='mggenes',
                                 legend.pos=c('inside', 'outside'),
                                 height=NULL, width=NULL, trim=0.02,
                                 colors = NULL,
                                 .plot_static = FALSE,
                                 square=TRUE, ...) {
  stopifnot(is(x, 'SparrowResult'))
  legend.pos <- match.arg(legend.pos)
  gs.dat <- setDF(subset(dat, group == 'geneset'))
  cols <- c('bg'='black', 'geneset'='red',
            'notsig'='grey', 'psig'='lightblue', 'sig'='darkblue')

  colors.default <- list(
    zeroline = "lightgrey",
    bg = "darkgrey",
    geneset = "darkorange",
    genes = "darkorange"
  )
  if (is.list(colors)) {
    colors <- c(colors, colors.default)
    colors[!duplicated(names(colors))]
  } else{
    colors <- colors.default  
  }
  
  if (value == 't') {
    gs.dat$y <- 0.0015 + runif(nrow(gs.dat), 0, 0.004)
    jitter <- 0.005
  } else {
    ## gs.dat$y <- c('notsig'=0.1, 'psig'=0.2, 'sig'=0.3)[gs.dat$significant]
    gs.dat$y <- 0.005 + runif(nrow(gs.dat), 0, 0.040)
    jitter <- 0.05
  }

  if (is.null(names(value))) {
    label <- if (value == "t") "t-statistic" else value
  } else {
    label <- names(value)
  }

  bg <- subset(dat, group == 'bg')

  bgd <- density(bg$val)
  gsd <- density(gs.dat$val)
  lmeta <- list(width=3)

  if (is.numeric(trim) && trim != 0) {
    xrange <- quantile(dat$val, c(trim, 1 - trim))
    xrange[1] <- min(xrange[1], min(gs.dat$val), min(gsd$x))
    xrange[2] <- max(xrange[1], max(gs.dat$val), max(gsd$x))
  } else {
    xrange <- range(dat$val)
  }
  if (square) {
    extreme <- max(abs(xrange))
    xrange <- c(-extreme, extreme)
  }

  if (.plot_static) {
    p <- ggplot2::ggplot() +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = colors$zeroline
      ) +
      ggplot2::geom_density(
        ggplot2::aes(x = val),
        color = colors$bg,
        linewidth = 1,
        data = dat
      ) +
      ggplot2::geom_density(
        ggplot2::aes(x = val),
        color = colors$geneset,
        linewidth = 1,
        data = gs.dat
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = val, y = y),
        size = 0.8,
        data = gs.dat,
        color = colors$genes
      ) +
      ggplot2::xlim(xrange) +
      ggplot2::labs(
        title = main,
        y = "Density",
        x = label
      )
  } else {
    p <- plot_ly(source=shiny_source, width=width, height=height)
    p <- add_lines(p, x=bgd$x, y=bgd$y, name='All Genes', hoverinfo='none', line=lmeta)
    p <- add_lines(p, x=gsd$x, y=gsd$y, name='Geneset', hoverinfo='none', line=lmeta)
    p <- layout(p, xaxis=list(title=label, range=xrange),
                yaxis=list(title="Density"),
                showlegend = with.legend, title=main, dragmode="select")
    if ('symbol' %in% names(gs.dat) && with.points) {
      p <- add_markers(p, x=~val, y=~y, key=~feature_id, data=gs.dat, name="Genes",
                       hoverinfo='text',
                       text=~paste0('Symbol: ', symbol, '<br>',
                                    'logFC: ', sprintf('%.3f', logFC), '<br>',
                                    'FDR: ', sprintf('%.3f', padj)))
    } else if (with.points) {
      p <- add_markers(p, x=~val, y=~y, key=~feature_id, data=gs.dat,
                       name="Genes",
                       text=~paste0('feature_id: ', feature_id, '<br>',
                                    'logFC: ', logFC, '<br>',
                                    'FDR: ', padj))
    }
    if (legend.pos == 'inside') {
      p <- layout(p, legend=list(x=0.75, y=1))
    }
    
    config(p, displaylogo=FALSE)    
  }
  
  p
}

#' @noRd
#' @importFrom plotly config ggplotly layout plotly_build
#' @importFrom ggplot2 aes geom_boxplot geom_jitter ggplot labs
iplot.boxplot.plotly <- function(x, dat, value, main, with.legend=TRUE,
                                 with.points=TRUE, shiny_source='mggenes',
                                 height=NULL, width=NULL, ggtheme=theme_bw(),
                                 trim=0.02, ...) {
  is.gs <- dat[['group']] == 'geneset'
  gs <- setDF(subset(dat, is.gs))
  bg <- setDF(subset(dat, !is.gs))
  n.gs <- sum(is.gs)

  if (is.null(names(value))) {
    label <- if (value == "t") "t-statistic" else value
  } else {
    label <- names(value)
  }

  # silence R CMD check NOTEs
  val <- symbol <- NULL
  gg <- ggplot(bg, aes(group, val)) +
    geom_boxplot(data=bg) +
    geom_boxplot(outlier.shape=NA, data=gs) +
    labs(
      title = main,
      x = NULL,
      y = label)
  if ('symbol' %in% names(bg) && with.points) {
    gg <- gg +
      suppressWarnings({
        geom_jitter(aes(key=feature_id,
                        text=paste0('Symbol: ', symbol, '<br>',
                                    'logFC: ', sprintf('%.3f', logFC), '<br>',
                                    'FDR: ', sprintf('%.3f', padj))),
                    data=gs, width=0.2)
      })

  } else if (with.points) {
    gg <- gg +
      suppressWarnings({
        geom_jitter(aes(key=feature_id,
                        text=paste(paste0('feature_id: ', feature_id, '<br>',
                                          'logFC: ', logFC, '<br>',
                                          'FDR: ', padj))),
                    data=gs, width=0.2)
      })
  }

  if (is(ggtheme, 'theme')) {
    gg <- gg + ggtheme
  }

  ## ggplotly keeps suggesting to use the github/ggplot2
  p <- suppressMessages({
    ggplotly(gg, width = width, height = height, tooltip = "text")
  })
  p <- layout(p, # yaxis = list(title=label),
              dragmode = "select", showlegend = with.legend)
  p <- plotly_build(p)

  p$x$source <- shiny_source
  ## Hacks to hide hover events on boxplots and remove outliers
  ## https://plot.ly/ggplot2/box-plots/#outliers
  p$x$data[[1]]$hoverinfo <- 'none'
  p$x$data[[1]]$marker <- list(opacity=0)
  p$x$data[[2]]$hoverinfo <- 'none'
  p$x$data[[2]]$marker <- list(opacity=0)
  config(p, displaylogo = FALSE)
}
