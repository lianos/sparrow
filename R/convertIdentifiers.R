#' Converts internal feature identifiers in a GeneSetDb to a set of new ones.
#'
#' The various GeneSetDb data providers limit the identifier types that they
#' provide. Use this function to map the given identifiers to whichever type
#' you like.
#'
#' For best results, provide your own identifier mapping reference, but we
#' provide a convenience wrapper around the [babelgene::orthologs()] function to
#' change between identifier types and species.
#'
#' When there are multiple target id's for the source id, they will all be
#' returned. When there is no target id for the source id, the soucre feature
#' will be axed.

#' @section Custom Mapping:
#' You need to provide a data.frame that has a column for the current
#' identifiers and another column for the target identifiers.
#'
#' @section Convenience identifier and species mapping:
#' If you don't provide a data.frame, you can provide a species name. We will
#' rely on the {babelgene} package for the conversion, so you will have to
#' provide a species name that it recognizes.
#'
#' @export
#' @param x The GeneSetDb with identifiers to convert
#' @param from,to If you are doing identifier and/orspecies conversion using
#'   babelgene, `to` is the species you want to convert to, and `from` is the
#'   species of `x`. If you are only doing id type conversion within the same
#'   species, specify the current species in `from`.
#'   If you are providing a data.frame map of identifiers in `xref`, `to` is
#'   the name of the column that holds the new identifiers, and `from` is the
#'   name of the column that holds the current identifiers.
#' @param id.type If you are using babelgene conversion, this specifies the
#'   type of identifier you want to convert to. It can be any of `"ensembl"`,
#'   `"entrez"`, or `"symbol"`.
#' @param extra.cols a character vector of columns from `to` to add to the
#'   features of the new GeneSetDb. If you want to keep the original identifiers
#'   of the remapped features, include `"original_id"` as one of the values
#'   here.
#' @param min_support,top Parameters used in the internal call to
#'   [babelgene::orthologs()]
#' @return A new GeneSetDb object with converted identifiers. We try to retain
#'   any metadata in the original object, but no guarantees are given. If
#'   `id_type` was stored previously in the collectionMetadata, that will be
#'   dropped.
#' @examples
#' # You can convert the identifiers within a GeneSetDb to some other type
#' # by providing a "translation" table. Check out the unit tests for more
#' # examples.
#' gdb <- exampleGeneSetDb() # this has no symbols in it
#'
#' # Define a silly conversion table.
#' xref <- data.frame(
#'   current_id = featureIds(gdb),
#'   new_id = paste0(featureIds(gdb), "_symbol"))
#' gdb2 <- convertIdentifiers(gdb, from = "current_id", to = "new_id",
#'                            xref = xref, extra.cols = "original_id")
#' geneSet(gdb2, name = "BIOCARTA_AGPCR_PATHWAY")
convertIdentifiers <- function(x, from = NULL, to = NULL,
                               id.type = c("ensembl", "entrez", "symbol"),
                               xref = NULL, extra.cols = NULL,
                               min_support = 3, top = TRUE) {
  assert_class(x, "GeneSetDb")
  id.type <- match.arg(id.type)
  if (test_string(to) && is.null(xref)) {
    stop("Automated species conversion by babelgene is not yet done",
         "https://github.com/lianos/sparrow/issues/2")
    bres <- .prep_babelgene_table(featureIds(x), to, id.type, is.human,
                                  min_support, top)
    xref <- bres[["table"]]
    from <- bres[["id.col"]]
    to <- bres[["target.col"]]
  }
  assert_multi_class(xref, c("data.frame", "data.table", "tbl"))
  assert_string(from)
  assert_string(to)
  take.cols <- c(from, to)
  if (!is.null(extra.cols) && length(extra.cols) > 0) {
    keep.original <- "original_id" %in% extra.cols
    take.cols <- c(take.cols, setdiff(extra.cols, "original_id"))
  }
  assert_subset(take.cols, colnames(xref))
  if (is.data.table(xref)) {
    xref <- xref[, take.cols, with = FALSE]
  } else {
    xref <- xref[, take.cols]
  }
  if (to == "feature_id") {
    setnames(xref, to, "fid.new")
    to <- "fid.new"
  }
  db <- merge(x@db, xref, by.x = "feature_id", by.y = from,
              suffixes = c(".original", ""))
  setnames(db, "feature_id", "original_id")
  setnames(db, to, "feature_id")

  # handle non std eval NOTE in R CMD check when using `:=` mojo
  N <- n <- active <- name <- NULL

  db <- db[!is.na(feature_id) & nchar(feature_id) > 0]
  db <- unique(db, by = c("collection", "name", "feature_id"))
  gs.dt <- merge(db, geneSets(x, as.dt = TRUE), by = c("collection", "name"))
  gs.dt[, N := NULL]
  gs.dt[, n := NULL]
  gs.dt[, active := NULL]
  if (!keep.original) {
    gs.dt[, original_id := NULL]
  }

  out <- GeneSetDb(gs.dt)
  out@collectionMetadata <- x@collectionMetadata[name != "id_type"]
  out
}

#' Internal helper function to handle bookkeeping tasks invovled to enable
#' species conversion from within convertIdentifiers
#'
#' TODO: Implement species conversion book keeping code for convertIdentifiers
#'
#' @noRd
.prep_babelgene_table <- function(ids, species, id.type, is.human,
                                  min_support, top) {
  if (FALSE) {
    x <- exampleGeneSetDb()
    ids <- featureIds(x)
    species <- "rat"
    id.type <- "ensembl"
    is.human <- TRUE
    min_support <- 3
    top <- TRUE

    ids <- c("P2ry12", "Trem2")
    is.human <- FALSE

    ids <- c("ENSMUSG00000036353", "ENSMUSG00000023992")
    species <- "rat"
    is.human <- FALSE
  }
  # orthologs will always return a data.frame with the first 3 columns
  # being human info: human_symbol, human_entrez, human_ensembl
  if (!is.human && species != "human") {
    # If you're query isn't from or to human, you have to stick human in the
    # middle
    # 1. map query id's to human
    human <- babelgene::orthologs(ids, species, human = FALSE)
    # 2. map human ids to target
    #
    xmap <- babelgene::orthologs(, species, human = FALSE)
  }

  xmap <- babelgene::orthologs(ids, species, human = is.human,
                               min_support = min_support, top = top)
}
