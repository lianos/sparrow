
#' Fetches a `GeneSetDb` from geneset collections defined in MSigDB.
#'
#' This provides versioned genesets from gene set collections defined in
#' [MSigDB](http://software.broadinstitute.org/gsea/msigdb). Collections can
#' be retrieved by their collection name, ie `c("H", "C2", "C7")`.
#'
#' @section KEGG Gene Sets:
#' Due to the licensing restrictions over the KEGG collections, they are not
#' returned from this function unless they are explicitly asked for. You can
#' ask for them through this function by either (i) querying for the `"c2"`
#' collection while setting `with.kegg = TRUE`; or (ii) explicitly calling with
#' `collection = "kegg"`.

#' @section MSigDB Versions:
#' We recently switched to using the msigdbr package as the source of truth for
#' these, so v7 is the earliest version of the MSigDB collections we make
#' available. Version 6 are available in the following (deprecated) packages:
#'
#' * https://github.com/lianos/GeneSetDb.MSigDB.Mmusculus.v61
#' * https://github.com/lianos/GeneSetDb.MSigDB.Hsapiens.v61
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
#' @importFrom msigdb.data msigdb_retrieve
#'
#' @param collection character vector specifying the collections you want
#'   (c1, c2, ..., c7, h). By default we load just the hallmark collecitons.
#'   Setting this to `NULL` loads all collections. Alternative you can also
#'   include named subsets of collections, like `"reactome"`. Refer to the
#'   Details section for more information.
#' @param species human or mouse?
#' @param with.kegg The Broad distributes the latest versions of the KEGG
#'   genesets as part of the c2 collection. These genesets come with a
#'   restricted license, so by default we do not return them as part of the
#'   GeneSetDb. To include the KEGG gene sets when asking for the c2
#'   collection, set this flag to `TRUE`.
#' @param allow_multimap,min_ortho_sources configure how to handle orthology
#'   mapping (allow multimappers, and what type of level of db suport required).
#'   See help in [msigdb.data::msigdb_retrieve()]
#' @param version the version of the MSigDB database to use.
#' @return a `GeneSetDb` object
#' @examples
#' \dontrun{
#'   gdb <- getMSigGeneSetDb(c("h", "reactome"), "human", "entrez")
#'   gdb.h.entrez <- getMSigGeneSetDb(c("h", "c2"), "human", "entrez")
#'   gdb.h.ens <- getMSigGeneSetDb(c("h", "c2"), "human", "ensembl")
#'   gdb.m.entrez <- getMSigGeneSetDb(c("h", "c2"), "mouse", "entrez")
#' }
getMSigGeneSetDb <- function(collection = NULL,
                             species = "human",
                             id.type = c("ensembl", "entrez", "symbol"),
                             with.kegg = FALSE,
                             allow_multimap = TRUE, min_ortho_sources = 2,
                             promote_subcategory_to_collection = FALSE,
                             prefix_collection = FALSE,
                             version = NULL, ...) {
  id.type <- match.arg(id.type)
  species.info <- species_info(species)
  valid.cols <- c("H", paste0("C", 1:8))
  if (species.info$alias != "human") valid.cols <- setdiff(valid.cols, "C1")
  if (!is.null(collection)) {
    collection <- assert_subset(toupper(collection), valid.cols)
  }

  sigs.all <- copy(.pkgcache[["msigdb"]][[species.info$species]])
  if (is.null(sigs.all)) {
    sigs.all <- as.data.table(msigdbr::msigdbr(species.info$species))
    axe.cols <- c("gs_pmid", "gs_geoid", "gs_url",
                  "gs_description", "species_name", "species_common_name",
                  "ortholog_sources", "num_ortholog_sources")
    axe.cols <- intersect(axe.cols, colnames(sigs.all))
    for (axe in axe.cols) sigs.all[, (axe) := NULL]
    setkeyv(sigs.all, c("gs_cat", "gs_name"))
    .pkgcache[["msigdb"]][[species.info$species]] <- copy(sigs.all)
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

  # Beef up collectionMetadata
  url.fn <- function(collection, name, ...) {
    url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"
    sprintf(url, name)
  }

  for (col in unique(geneSets(gdb)$collection)) {
    geneSetCollectionURLfunction(gdb, col) <- url.fn
    featureIdType(gdb, col) <- idtype
    gdb <- addCollectionMetadata(
      gdb, col, 'source', as.character(packageVersion("msigdbr")))
  }

  org(gdb) <- gsub(" ", "_", species.info[["species"]])
  gdb@collectionMetadata <- gdb@collectionMetadata[name != "count"]
  gdb
}

.msigdb.cache <- new.env()
