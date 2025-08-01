#' Rename the collections in a GeneSetDb
#'
#' This function remaps names of collections in the database from their current
#' names to ones specified by the user, folows the `dplyr::rename` convenction
#' where `names()` of the rename vector are the new names you want, and its
#' values are the old names it came from.
#'
#' @export
#' @param x A GeneSetDb object
#' @param rename a named character vector. `names(rename)` are the new names
#'   names for the collection, and the values are the current names of the 
#'   collection.
#' @param ... pass it along
#' @return GeneSetDb `x` with renamed `geneSets(x)$collection` values.
#' @examples
#' gdb <- exampleGeneSetDb()
#' ngdb <- renameCollection(gdb, c("MSigDB C2" = "c2", "ImmuneSigDb" = "c7"))
#' all.equal(
#'   unname(geneSetURL(gdb, "c7", "GSE3982_BCELL_VS_TH2_DN")),
#'   unname(geneSetURL(ngdb, "ImmuneSigDb", "GSE3982_BCELL_VS_TH2_DN")))
renameCollection <- function(x, rename = NULL, ...) {
  if (is.null(rename)) return(x)
  
  assert_character(rename, names = "unique")
  found <- rename %in% geneSets(x)[["collection"]]
  if (any(!found)) {
    warning("These collections you want to rename were not found: ",
            paste(names(rename[!found])), collapse = ",")
    rename <- rename[found]
  }
  if (length(rename) == 0L) return(x)
  
  out <- x
  out@db <- data.table::copy(out@db)
  out@collectionMetadata <- data.table::copy(out@collectionMetadata)
  out@table <- data.table::copy(out@table)
  rename <- setNames(names(rename), rename)
  
  for (oname in names(rename)) {
    nname <- rename[[oname]]
    out@db[collection == oname, collection := nname]
    out@table[collection == oname, collection := nname]
    out@collectionMetadata[collection == oname, collection := nname]
  }
  setkeyv(out@db, data.table::key(x@db))
  setkeyv(out@table, data.table::key(x@table))
  setkeyv(out@collectionMetadata, data.table::key(x@collectionMetadata))
  
  out
}

#' @noRd
#' @export
renameCollections <- function(x, rename = NULL, ...) {
  .Deprecated("renameCollection", msg = "remove the last S for savings")
  renameCollection(x, rename = rename, ...)
}
