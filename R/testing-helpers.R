# The datasets references in this file are generated in:
#   tests/testdata/setup-testdata.R

#' Functions that load data for use in examples and testing.
#'
#' We provide examplar expression data (counts or voomed) as well as exemplar
#' gene sets in different forms.
#'
#' @section exampleExpressionSet:
#' The expression data is a subset of the TCGA BRCA indication. Calling
#' `exampleExpressionSet(do.voom = TRUE)` will return a voomed `EList` version
#' of the data. When `do.voom = FALSE`, you will get a DGEList of the counts
#'
#' @rdname examples
#' @aliases exampleExpressionSet
#'
#' @export
#' @importFrom limma voom
#'
#' @param dataset Character vector indicating what samples wanted, either
#'   \code{"tumor-vs-normal"} for a tumor vs normal dataset from TCGA, or
#'   just the tumor samples from the same annotated with subtype.
#' @param do.voom If TRUE, a voomed EList is returned, otherwise an
#'   ExpressionSet of counts.
#' @examples
#' vm <- exampleExpressionSet()
exampleExpressionSet <- function(dataset=c('tumor-vs-normal', 'tumor-subtype'),
                                 do.voom=TRUE) {
  dataset <- match.arg(dataset)
  fn <- system.file("extdata", "testdata", "TCGA-BRCA-some.DGEList.rds",
                    package = "sparrow", mustWork = TRUE)
  y.all <- readRDS(fn)

  # Two samples seem to be outliers:
  axe.samples <- c("TCGA-A2-A3XV-01A-21R-A239-07", "TCGA-A2-A3XU-01A-12R-A22U-07")
  y.all <- y.all[, !colnames(y.all) %in% axe.samples]

  if (dataset == 'tumor-vs-normal') {
    y <- y.all
    design <- model.matrix(~ Cancer_Status, y$samples)
    colnames(design) <- sub('Cancer_Status', '', colnames(design))
  } else {
    y <- y.all[, y.all$samples$Cancer_Status == 'tumor']
    y$samples <- droplevels(y$samples)
    design <- model.matrix(~ 0 + PAM50subtype, y$samples)
    colnames(design) <- sub('PAM50subtype', '', colnames(design))
  }

  y$design <- design

  if (do.voom) voom(y, y$design, plot = FALSE) else y
}


#' @section exampleGeneSets:
#' Returns gene sets as either a list of feature identifiers. Entrez identifiers
#' are used. If `x` is provided, integers that index into the expression
#' container `x` are used (this is a legacy feature that we should nuke).
#'
#' @rdname examples
#' @aliases exampleGeneSets
#'
#' @export
#'
#' @param x If provided, an expression/matrix object so that the genesets are
#'   returned as (integer) index vectors into the rows of x whose rownames
#'   match the ids in the geneset.
#' @param unlist return the genesets as nested list of lists (default: `TRUE`).
#'   The top level lists corresponds to the collection, and the lists within
#'   each are the inidividual gene sets. If `FALSE`, a single list of genesets
#'   is returned.
#' @return A list of lists of entrezIDs when `as == 'lol'`, or
#'   a list of integers into the rows of `x`.
#' @examples
#' head(exampleGeneSets())
exampleGeneSets <- function(x, unlist = !missing(x)) {
  gsl.fn <- system.file('extdata', 'testdata',
                        'genesets-sparrow-list-of-lists.rds',
                        package='sparrow')
  gsl <- readRDS(gsl.fn)
  if (!missing(x)) {
    gsl <- lapply(gsl, function(col) {
      lapply(col, function(gs) {
        out <- match(gs, rownames(x))
        out[!is.na(out)]
      })
    })
  }
  if (unlist) {
    gsl <- unlist(gsl, recursive=FALSE)
    names(gsl) <- sub('\\.', ';;', names(gsl))
  }
  gsl
}

#' @section exampleGeneSetDb:
#' Returns gene sets as a `GeneSetDb` object
#'
#' @rdname examples
#' @aliases exampleGeneSetDb
#' @export
exampleGeneSetDb <- function() {
  out <- GeneSetDb(exampleGeneSets())
  colls <- unique(collectionMetadata(out, as.dt = TRUE)$collection)
  for (col in colls) {
    geneSetCollectionURLfunction(out, col) <- ".geneSetURL.msigdb"
  }
  out
}

#' @section exampleGeneSetDF:
#' Returns a data.frame of gene set definitions. A data.frame of this form
#' can be passed into the `GeneSetDb()` contructor.
#'
#' @export
#' @importFrom utils read.csv
#' @rdname examples
#' @aliases exampleGeneSetDF
exampleGeneSetDF <- function() {
  gs.df <- system.file('extdata', 'testdata', 'custom-sigs.csv',
                       package='sparrow')
  read.csv(gs.df, stringsAsFactors = FALSE, colClasses = 'character')
}


#' @export
#' @rdname examples
#' @aliases exampleGeneSetDF
#' @param cached If `TRUE` (default), returns a pre-saved SparrowResult object.
#'   Otherwise calculates a fresh one using the `methods` provided
#' @param methods the methods to use to create a new SparrowResult for.
exampleSparrowResult <- function(cached = TRUE,
                                 methods = c("cameraPR", "fry")) {
  if (cached) {
    fn <- system.file('extdata', 'testdata', 'test-SparrowResult.rds',
                      package='sparrow')
    out <- readRDS(fn)
  } else {
    vm <- exampleExpressionSet()
    gdb <- exampleGeneSetDb()
    out <- seas(vm, gdb, methods = methods,
                design = vm$design, contrast = "tumor")
  }
  out
}

#' `exampleDgeResult` returns a data.frame of differential expression results.
#' Currently they are for human ensembl genes. Setting the `induce.bias`
#' parameter to `"effective_length"` or `"AveExpr"` will munge the returned
#' result such that larger "bias" values will be associated to lower pvalues,
#' so we can more easily test biased overrepresentation analysis approaches like
#' [ora()] and [goseq()].
#'
#' @export
#' @param species the species to return the example result from (right now,
#'   only "human")
#' @param id.type the type of identifiers to use (right now, only "ensembl")
#' @param induce.bias We can simulate a bias on the pvalue by the gene's
#'   `"effective_length"` or `"AveExpr"`. These are columns that are included
#'   in the output. If `NULL`, no bias is introduced into the result.
#' @rdname examples
exampleDgeResult <- function(species = "human", id.type = "ensembl",
                             induce.bias = NULL) {
  # we only have human/ensembl for now
  species <- match.arg(species, "human")
  id.type <- match.arg(id.type, "ensembl")
  dge.fn <- system.file(
    "extdata", "testdata",
    "dataframe-input-short.csv.gz",
    package = "sparrow")
  out <- data.table::fread(dge.fn, data.table = FALSE)
  if (is.character(induce.bias)) {
    bias <- match.arg(induce.bias, c("effective_length", "AveExpr"))
    o <- order(out[["pval"]])
    out[[bias]][o] <- sort(out[[bias]], decreasing = TRUE)
  }
  out$significant <- out$selected
  out
}

#' Generates a fake GeneSetDb by sampling from features in a seas input.
#'
#' I wrote this because initial fetching from msigdbr can be slow, and also
#' having some weird crashes in the unit tests of bioc3.14-devel.
#'
#' This is a helper function for development, and shouldn't be used by normal
#' users of this package
#'
#' @export
#' @param x an input container to [seas()]
#' @param n number of genesets
#' @param bias column in `x` to bias the geneset creation by
#' @export
#' @examples
#' gdb.rando <- randomGeneSetDb(exampleDgeResult(), 10, bias = "t")
randomGeneSetDb <- function(x, n = 10, bias = NULL, seed = 10, ...) {
  assert_class(x, "data.frame") # only data.frames for now
  assert_number(n, lower = 2, upper = 100)
  set.seed(seed)
  gsets <- lapply(seq_len(n), function(i) {
    idx <- sample(nrow(x), 10, prob = abs(x$t))
    wtf <- data.frame(
      collection = rep("random", 10),
      name = rep(paste0("geneset", i), 10),
      feature_id = x$feature_id[idx],
      symbol = x$symbol[idx])
  })
  GeneSetDb(do.call(rbind, gsets))
}
