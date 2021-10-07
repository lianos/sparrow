
#' Fetches gene set collections from the moleular signature database (MSigDB)
#'
#' This provides versioned genesets from gene set collections defined in
#' [MSigDB](http://software.broadinstitute.org/gsea/msigdb). Collections can
#' be retrieved by their collection name, ie `c("H", "C2", "C7")`.
#'
#' @section Species and Identifier types:
#' This function utilizes the functionality from the `{msigdbr}` and
#' `{babelgene}` packages to retrieve gene set definitions from a variety of
#' organisms and identifier types.
#'
#' @section KEGG Gene Sets:
#' Due to the licensing restrictions over the KEGG collections, they are not
#' returned from this function unless they are explicitly asked for. You can
#' ask for them through this function by either (i) querying for the `"c2"`
#' collection while setting `with.kegg = TRUE`; or (ii) explicitly calling with
#' `collection = "kegg"`.
#'
#' @section Citing the Molecular Signatures Database:
#' To cite your use of the Molecular Signatures Database (MSigDB), please
#' reference Subramanian, Tamayo, et al. (2005, PNAS 102, 15545-15550) and one
#' or more of the following as appropriate:
#'
#' * Liberzon, et al. (2011, Bionformatics);
#' * Liberzon, et al. (2015, Cell Systems); and
#' * The source for the gene set as listed on the gene set page.
#'
#' @export
#'
#' @param collection character vector specifying the collections you want
#'   (c1, c2, ..., c7, h). By default we load just the hallmark collecitons.
#'   Setting this to `NULL` loads all collections. Alternative you can also
#'   include named subsets of collections, like `"reactome"`. Refer to the
#'   Details section for more information.
#' @param species `"human"` or `"mouse"`? Really, this is anything available
#'   in the `alias` column of the `sparrow:::species_info()` table (except
#'   cyno).
#' @param id.type do you want the feature id's used in the gene sets to be
#'   `"ensembl"` (default), `"entrez"`, or `"symbol"`.
#' @param with.kegg The Broad distributes the latest versions of the KEGG
#'   genesets as part of the c2 collection. These genesets come with a
#'   restricted license, so by default we do not return them as part of the
#'   GeneSetDb. To include the KEGG gene sets when asking for the c2
#'   collection, set this flag to `TRUE`.
#' @param promote_subcategory_to_collection there are different sources of
#'   genesets for a number of the collections in MSigDB. These are included
#'   in the `gs_subcat` column of `geneSets(this)`. When this is set to `TRUE`,
#'   the collection column for the genesets is appended with the subcatory.
#'   So, instead of having a massive `"C2"` collection, you'll have bunch of
#'   collections like `"C2_CGP"`, `"C2_CP:BIOCARTA"`, etc.
#' @param prefix_collection When `TRUE` (default: `FALSE`), the `"C1"`, `"C2"`,
#'   etc. is prefixed with `"MSigDB_*"`
#' @param ... pass through parameters
#' @return a `BiocSet` of the MSigDB collections
#' @examples
#' \donttest{
#'   # these take a while to load initially, so put them in dontrun blocks.
#'   # you should run these interactively to understand what they return
#'   bcs <- getMSigDbCollection("h", "human", "entrez")
#'   bcs.h.entrez <- getMSigDbCollection(c("h", "c2"), "human", "entrez")
#'   bcs.h.ens <- getMSigDbCollection(c("h", "c2"), "human", "ensembl")
#'   bcs.m.entrez <- getMSigDbCollection(c("h", "c2"), "mouse", "entrez")
#'
#'   gdb <- getMSigGeneSetDb("h", "human", "entrez")
#' }
getMSigDbCollection <- function(collection = NULL,
                                species = "human",
                                id.type = c("ensembl", "entrez", "symbol"),
                                with.kegg = FALSE,
                                promote_subcategory_to_collection = FALSE,
                                prefix_collection = TRUE, ...) {
  id.type <- match.arg(id.type)
  out <- getMSigGeneSetDb(
    collection, species, id.type, with.kegg,
    promote_subcategory_to_collection, prefix_collection, ...)
  as(out, "BiocSet")
}


#' @describeIn getMSigDbCollection retrieval method for a GeneSetDb container
#' @export
getMSigGeneSetDb <- function(collection = NULL,
                             species = "human",
                             id.type = c("ensembl", "entrez", "symbol"),
                             with.kegg = FALSE,
                             promote_subcategory_to_collection = FALSE,
                             prefix_collection = FALSE, ...) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("The msigdbr package is required for this functionality")
  }
  id.type <- match.arg(id.type)
  species.info <- species_info(species)
  valid.cols <- c("H", paste0("C", 1:8))
  if (species.info$alias != "human") valid.cols <- setdiff(valid.cols, "C1")
  if (!is.null(collection)) {
    collection <- assert_subset(toupper(collection), valid.cols)
  }

  # handle non std eval NOTE in R CMD check when using `:=` mojo
  # each of these variables are referenced in some data.table NSE mojo below
  gs_cat <- gs_subcat <- gs_name <- symbol <- ensembl_id <- gs_id <- NULL

  sigs.all <- data.table::copy(.pkgcache[["msigdb"]][[species.info$species]])
  if (is.null(sigs.all)) {
    sigs.all <- msigdbr::msigdbr(species.info$species)
    sigs.all <- as.data.table(sigs.all)
    axe.cols <- c("gs_pmid", "gs_geoid", "gs_url",
                  "gs_description", "species_name", "species_common_name",
                  "ortholog_sources", "num_ortholog_sources")
    axe.cols <- intersect(axe.cols, colnames(sigs.all))
    for (axe in axe.cols) sigs.all[, (axe) := NULL]
    setkeyv(sigs.all, c("gs_cat", "gs_name"))
    .pkgcache[["msigdb"]][[species.info$species]] <- data.table::copy(sigs.all)
  }

  if (!is.null(collection)) {
    out <- sigs.all[gs_cat %in% collection]
  } else {
    out <- sigs.all
  }

  if (!with.kegg) {
    out <- out[gs_subcat != "CP:KEGG"]
  }

  if (prefix_collection) {
    out[, collection := paste0("MSigDB_", out$gs_cat)]
  } else {
    out[, collection := gs_cat]
  }

  if (promote_subcategory_to_collection) {
    out[, collection := {
      ifelse(nchar(out$gs_subcat) == 0L,
             out$collection,
             paste(out$collection, out$gs_subcat, sep = "_"))
    }]
  }

  if (id.type == "ensembl") {
    idtype <- GSEABase::ENSEMBLIdentifier()
    idcol <- "ensembl_gene"
  } else if (id.type == "entrez") {
    idtype <- GSEABase::EntrezIdentifier()
    idcol <- "entrez_gene"
  } else {
    idtype <- GSEABase::SymbolIdentifier()
    idcol <- "gene_symbol"
  }

  ret <- out[, {
    list(collection, name = gs_name, feature_id = as.character(.SD[[idcol]]),
         subcategory = gs_subcat)
  }, .SDcols = c("collection", "gs_name", idcol)]
  if (id.type != "symbol") {
    ret[, symbol := out[["gene_symbol"]]]
  } else {
    ret[, ensembl_id := out[["ensembl_gene"]]]
  }
  ret[, gs_id := {
    ifelse(grepl("GO:", out$gs_exact_source), out$gs_exact_source, out$gs_id)
  }]

  ret <- ret[!is.na(feature_id)]
  ret <- unique(ret, by = c("collection", "name", "feature_id"))
  gdb <- GeneSetDb(ret)

  for (col in unique(geneSets(gdb)$collection)) {
    geneSetCollectionURLfunction(gdb, col) <- ".geneSetURL.msigdb"
    featureIdType(gdb, col) <- idtype
    gdb <- addCollectionMetadata(
      gdb, col, 'source', as.character(packageVersion("msigdbr")))
  }

  # org(gdb) <- gsub(" ", "_", species.info[["species"]])
  gdb@collectionMetadata <- gdb@collectionMetadata[name != "count"]
  gdb
}


#' @noRd
.geneSetURL.msigdb <- function(collection, name, ...) {
  url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"
  sprintf(url, name)
}

