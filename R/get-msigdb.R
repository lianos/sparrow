
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
#' @param promote.subcollection there are different sources of
#'   genesets for a number of the collections in MSigDB. These are included
#'   in the `gs_subcollection` column of `geneSets(this)`. When this is set to
#'  `TRUE`, the collection column for the genesets is appended with the
#'   subcollection.
#'   So, instead of having a massive `"C2"` collection, you'll have bunch of
#'   collections like `"C2_CGP"`, `"C2_CP:BIOCARTA"`, etc.
#' @param prefix.collection When `TRUE` (default: `FALSE`), the `"C1"`, `"C2"`,
#'   etc. is prefixed with `"MSigDB_*"`
#' @param strip.subcollection.prefix removes the CGP: type prefixes for the
#'   `gs_subcollection` column, except for the C5 GO collection.
#' @param merge.human.into.mouse When `TRUE` (default), the OG human collections
#'   are merged into the newly  minted (as of 2024) M* mouse collections. Set
#'   to `FALSE` to not do that.
#' @param ... pass through parameters
#' @return a `BiocSet` of the MSigDB collections
#' @examples
#' \donttest{
#'   # these take a while to load initially, so put them in dontrun blocks.
#'   # you should run these interactively to understand what they return
#'   bcs <- getMSigCollection("h", "human", "entrez")
#'   bcs.h.entrez <- getMSigCollection(c("h", "c2"), "human", "entrez")
#'   bcs.h.ens <- getMSigCollection(c("h", "c2"), "human", "ensembl")
#'   bcs.m.entrez <- getMSigCollection(c("h", "c2"), "mouse", "entrez")
#'
#'   gdb <- getMSigGeneSetDb("h", "human", "entrez")
#' }
getMSigCollection <- function(collection = NULL,
                              species = "human",
                              id.type =  c("ensembl", "entrez", "symbol", "uniprot"),
                              with.kegg = FALSE,
                              promote.subcollection = FALSE,
                              prefix.collection = FALSE,
                              strip.subcollection.prefix = TRUE, 
                              merge.human.into.mouse = TRUE, ...) {
  id.type <- match.arg(id.type)
  out <- getMSigGeneSetDb(
    collection = collection,
    species = species, 
    id.type = id.type, 
    with.kegg = with.kegg,
    promote.subcollection = promote.subcollection,
    prefix.collection = prefix.collection, 
    strip.subcollection.prefix = strip.subcollection.prefix,
    merge.human.into.mouse = merge.human.into.mouse,
    ...)
  as(out, "BiocSet")
}

#' @describeIn getMSigCollection retrieval method for a GeneSetDb container
#' @export
getMSigGeneSetDb <- function(collection = NULL,
                             species = "human",
                             id.type = c("ensembl", "entrez", "symbol", "uniprot"),
                             with.kegg = FALSE,
                             promote.subcollection = FALSE,
                             prefix.collection = FALSE,
                             strip.subcollection.prefix = TRUE,
                             merge.human.into.mouse = TRUE,
                             refetch = FALSE, ...) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("The msigdbr package is required for this functionality")
  }
  id.type <- match.arg(id.type)
  species.info <- species_info(species)
  valid.cols <- c("H", paste0("C", 1:8))
  
  if (id.type == "uniprot") {
    stop("No uniprot support just yet")
  }
  
  # if (species.info$alias %in% c("mouse")) {
  #   # TODO: I think we should mouse geneset db for rats too, but msigdbr does
  #   # not convert between two other species, when human isn't one of them.
  #   valid.cols <- paste0("M", valid.cols)
  #   db.species <- "MM"
  # } else {
  #   db.species <- "HS"
  # }
  
  if (!species.info$alias %in% c("human", "mouse")) {
    valid.cols <- setdiff(valid.cols, "C1")
  }
  if (!is.null(collection)) {
    collection <- assert_subset(toupper(collection), valid.cols)
    if (species.info$alias == "mouse") paste0("M", collection)
  }
  
  # handle non std eval NOTE in R CMD check when using `:=` mojo
  # each of these variables are referenced in some data.table NSE mojo below
  gs_collection <- gs_subcollection <- gs_cat <- gs_subcat <- NULL
  gs_name <- symbol <- ensembl_id <- gs_id <- NULL
  
  cache.key <- species.info$species
  if (species.info$alias == "mouse" && merge.human.into.mouse) {
    cache.key <- paste0(cache.key, ".with_human")
  }
  sigs.all <- data.table::copy(.pkgcache[["msigdb"]][[cache.key]])
  if (is.null(sigs.all) || isTRUE(refetch)) {
    sigs.all <- setDT(msigdbr::msigdbr(species.info$species, db_species = "HS"))
    # Because I want to merge geneset collections/subcollections between human
    # and mouse, I am stripping gs_subcollections that have a `NNN:` prefix.
    # See the bottom of the test-gdb-msigdb.R, there are some collections that
    # have such a prefix in human db but not mouse, and these prefixes have no
    # semantic meaning that is used anyway
    
    if (strip.subcollection.prefix) {
      sigs.all[, gs_subcollection := {
        ifelse(
          gs_collection == "C5",
          gs_subcollection,
          sub(".*?:", "", gs_subcollection))
      }]
    }
    sigs.all[, ncbi_gene := as.character(ncbi_gene)] # sometimes this is an integer
    sigs.all[, db_species := "HS"]
    sigs.all[, msigdb_collection := gs_collection]
    if (species.info$alias == "mouse") {
      # Let's give mouse genesets precedence
      sm.all <- setDT(msigdbr::msigdbr(species.info$species, db_species = "MM"))
      sm.all[, msigdb_collection := gs_collection]
      sm.all[, db_species := "MM"]
      sm.all[, ncbi_gene := as.character(ncbi_gene)]
      sm.all[, gs_collection := {
        ifelse(
          gs_collection == "MH", 
          "H", 
          sub("^M", "C", gs_collection)
        )
      }]
      if (strip.subcollection.prefix) {
        sm.all[, gs_subcollection := {
          ifelse(
            gs_collection == "C5",
            gs_subcollection,
            sub(".*?:", "", gs_subcollection))
        }]
      }

      if (merge.human.into.mouse) {
        # sigs.hs <- sigs.all[!sm.all, on = c("gs_collection", "gs_subcollection", "gs_name")]
        # We can just antijoin on the gs_name
        
        sigs.hs <- sigs.all[!sm.all, on = "gs_name"]
        # Remove C5 -- the mouse db has all the mouse GO collections 
        sigs.hs <- sigs.hs[collection != "C5"]
        
        sigs.all <- rbindlist(list(sm.all, sigs.hs), fill = TRUE)
        sigs.all <- sigs.all[order(gs_collection, gs_subcollection, gs_name)]
      } else {
        sigs.all <- sm.all
      }
    }
    
    axe.cols <- c("gs_pmid", "gs_geoid",
                  "gs_description", "species_name", "species_common_name",
                  "ortholog_sources", "num_ortholog_sources")
    axe.cols <- intersect(axe.cols, colnames(sigs.all))
    for (axe in axe.cols) sigs.all[, (axe) := NULL]
    setkeyv(sigs.all, c("gs_collection", "gs_name"))
    .pkgcache[["msigdb"]][[cache.key]] <- data.table::copy(sigs.all)
  }
  
  if (!is.null(collection)) {
    out <- sigs.all[gs_collection %in% collection]
  } else {
    out <- sigs.all
  }
  
  if (!with.kegg) {
    # out <- out[gs_subcat != "CP:KEGG"]
    out <- out[!endsWith(gs_subcollection, "KEGG_LEGACY")]
  }
  
  if (prefix.collection) {
    out[, collection := paste0("MSigDB_", out$gs_collection)]
  } else {
    out[, collection := gs_collection]
  }
  
  if (promote.subcollection) {
    out[, collection := {
      ifelse(nchar(out$gs_subcollection) == 0L,
             out$collection,
             paste(out$collection, out$gs_subcollection, sep = "_"))
    }]
  }
  
  if (id.type == "ensembl") {
    # idtype <- GSEABase::ENSEMBLIdentifier()
    idcol <- "ensembl_gene"
  } else if (id.type == "entrez") {
    # idtype <- GSEABase::EntrezIdentifier()
    idcol <- "entrez_gene"
  } else {
    # idtype <- GSEABase::SymbolIdentifier()
    idcol <- "gene_symbol"
  }
  idtype <- id.type
  
  ret <- out[, {
    list(
      collection, 
      name = gs_name, 
      feature_id = as.character(.SD[[idcol]]),
      symbol = gene_symbol,
      subcollection = gs_subcollection,
      collection_name = gs_collection_name,
      geneset_id = gs_id,
      geneset_url = gs_url,
      msigdb_collection,
      db_species)
  }]
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
  
  # browser()
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


# @describeIn getMSigCollection retrieval method for a GeneSetDb container
# @export
getMSigGeneSetDb.old <- function(collection = NULL,
                             species = "human",
                             id.type = c("ensembl", "entrez", "symbol"),
                             with.kegg = FALSE,
                             promote.subcollection.to.collection = FALSE,
                             prefix.collection = FALSE, ...) {
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

  if (prefix.collection) {
    out[, collection := paste0("MSigDB_", out$gs_cat)]
  } else {
    out[, collection := gs_cat]
  }

  if (promote.subcollection.to.collection) {
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
         subcollection = gs_subcat)
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
.geneSetURL.msigdb <- function(collection, name, gdb, ...) {
  default.url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"
  default.url <- sprintf(default.url, name)

  gs <- geneSets(gdb)
  info <- gs[gs$collection == collection & gs$name == name,,drop = FALSE]
  if (nrow(info) == 1) {
    url <- info$geneset_url
  } else {
    url <- NA_character_
  }
  
  if (test_string(url, na.ok = FALSE, min.chars = 2)) url else default.url
}

