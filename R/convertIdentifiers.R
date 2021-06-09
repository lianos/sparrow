#' Converts the feature identifiers in a GeneSetDb to a set of new ones.
#'
#' The various GeneSetDb data providers limit the identifier types that they
#' provide. Use this function to map the given identifiers to whichever type
#' you like.
#'
#' You need to provide a data.frame that has a column for the current
#' identifiers and another column for the target identifiers.
#'
#' If you don't provide a data.frame we can use to translate identifiers,
#' we'll use hte babelgene package.
#'
#' When there are multiple target id's for the source id, they will all be
#' returned. When there is no target id for the source id, the soure feature
#' will be axed.
#'
#' By default, the first column of the `xref` data.frame is assumed to be
#' the source identifiers, and the second is the target_id you want to remap
#' the identifiers to.
#'
#' @export
#' @param x The GeneSetDb with identifiers to convert
#' @param xref If a species (`"mouse"` or `"human"`) we will load the internal
#'   ensembl <-> entrez data.frame, or you can provide your own.
#' @return a remapped GeneSetDb object
#' @examples
#' \dontrun{
#' # This function needs to be reimplemented using the {babelgene} package,
#' # and these examples are left for posterity to remind us of what once was.
#' gdb.entrez <- exampleGeneSetDb()
#' gdb.ens <- convertIdentifiers(gdb.entrez, "human",
#'                              original_id = "entrezgene_id",
#'                              target_id = "ensembl_gene_id")
#' }
convertIdentifiers <- function(x, xref, original_id = colnames(xref)[1L],
                               target_id = colnames(xref)[2L], ...) {
  stop("Re-implement this using the babelgene package\n",
       "https://github.com/lianos/sparrow/issues/2")
  assert_class(x, "GeneSetDb")
  if (test_string(xref)) xref <- load_id_xref(xref)
  assert_multi_class(xref, c("data.frame", "data.table", "tbl"))
  assert_string(original_id)
  assert_string(target_id)
  assert_subset(c(original_id, target_id), colnames(xref))

  db <- merge(x@db, xref, by.x = "feature_id", by.y = original_id,
              suffixes = c(".original", ""))
  setnames(db, "feature_id", "original_id")
  setnames(db, target_id, "feature_id")

  db <- db[!is.na(feature_id) & nchar(feature_id) > 0]
  db <- unique(db, by = c("collection", "name", "feature_id"))
  gs.dt <- merge(db, geneSets(x, as.dt = TRUE), by = c("collection", "name"))
  gs.dt[, N := NULL]
  gs.dt[, n := NULL]
  gs.dt[, active := NULL]

  out <- GeneSetDb(gs.dt)
  out@collectionMetadata <- x@collectionMetadata
  out
}
