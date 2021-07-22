#' Get pathways/GOslim information from PANTHER.db Biocondcutor package.
#'
#' @description
#' This is a convience function that orchestrates the PANTHER.db package to
#' return GeneSetDb objects for either pathway or GOslim information for
#' human or mouse.
#'
#' Note that for some reason the `PANTHER.db` package needs to be
#' installed in a user-writable package location for this to work properly.
#' If you see an error like "Error in resqlite_send_query ... attempt to
#' write a readonly database", this is the problem. Please install another
#' version of the `PANTHER.db` package in a user-writable directory using
#' [BiocManager::install()].
#'
#' @section GOSLIM:
#' [GO Slims](http://geneontology.org/page/go-slim-and-subset-guide) are
#' "cut down" versions of the GO ontology that contain a subset of the terms in
#' the whole GO.
#'
#' PANTHER provides their own set of
#' [GO slims](http://www.pantherdb.org/panther/ontologies.jsp), although
#' it's not clear how often these get updated.
#'
#' @rdname getPantherGeneSetDb
#' @export
#' @importFrom GSEABase EntrezIdentifier
#' @param type "pathway" or, "goslim"
#' @param species "human" or "mouse"
#'
#' @return A wired up GeneSetDb
getPantherGeneSetDb <- function(type=c('pathway', 'goslim'),
                                species=c('human', 'mouse')) {
  species <- match.arg(species)
  type <- match.arg(type)

  if (!requireNamespace('PANTHER.db', quietly=TRUE)) {
    stop("The PANTHER.db bioconductor package is required")
  }
  if (species == 'human') {
    org.pkg <- 'org.Hs.eg.db'
    xorg <- 'Homo_sapiens'
  } else {
    org.pkg <- 'org.Mm.eg.db'
    xorg <- 'Mus_musculus'
  }
  on.exit({
    unloadNamespace('PANTHER.db')
    unloadNamespace(org.pkg)
  })
  if (!requireNamespace(org.pkg, quietly=TRUE)) {
    stop(org.pkg, " bioconductor package required for this species query")
  }

  p.db <- PANTHER.db::PANTHER.db
  PANTHER.db::pthOrganisms(p.db) <- toupper(species)
  org.db <- getFromNamespace(org.pkg, org.pkg)

  out <- switch(type,
                pathway=getPantherPathways(p.db, org.db),
                goslim=getPantherGOSLIM(p.db, org.db))
  mapIds <- getFromNamespace('mapIds', 'AnnotationDbi')
  out@db$symbol <- mapIds(org.db, out@db$feature_id, 'SYMBOL', 'ENTREZID')
  org(out) <- xorg
  out
}

#' @rdname getPantherGeneSetDb
#' @export
getGOslimGeneSetDb <- function(species=c('human', 'mouse')) {
  getPantherGeneSetDb('goslim', species)
}

#' @noRd
getPantherPathways <- function(p.db, org.db) {
  cname <- "PANTHER pathway"
  aselect <- getFromNamespace('select', 'AnnotationDbi')
  p.all <- aselect(p.db, AnnotationDbi::keys(p.db, keytype="PATHWAY_ID"),
                   columns=c("PATHWAY_ID", "PATHWAY_TERM", "UNIPROT"),
                   'PATHWAY_ID')
  ## Map uniprot to entrez
  umap <- aselect(org.db, p.all$UNIPROT, c('UNIPROT', 'ENTREZID'), 'UNIPROT')
  m <- merge(p.all, umap, by='UNIPROT')
  m <- m[!is.na(m[['ENTREZID']]),,drop=FALSE]
  names(m) <- c("uniprot_id", "pathway_id", "name", "feature_id")

  gdb <- GeneSetDb(m, collection = cname)
  geneSetCollectionURLfunction(gdb, cname) <- ".geneSetURL.PANTHERpathway"
  featureIdType(gdb, cname) <- EntrezIdentifier()
  gdb
}

.geneSetURL.PANTHERpathway <- function(coll, gsname, gdb, ...) {
  info <- gdb@table[collection == coll & name == gsname]
  if (nrow(info) == 0) {
    "http://pantherdb.org/panther/prowler.jsp?reset=1&selectedView=5"
  } else {
    paste0("http://pantherdb.org/pathway/pathwayDiagram.jsp?catAccession=",
           info$pathway_id)
  }
}

#' @noRd
getPantherGOSLIM <- function(p.db, org.db) {
  cname <- "GOSLIM"
  if (!requireNamespace("GO.db")) {
    stop("GO.db is required for this functionality")
  }
  aselect <- getFromNamespace('select', 'AnnotationDbi')
  p.all <- aselect(p.db,
                   AnnotationDbi::keys(p.db, keytype='GOSLIM_ID'),
                   columns=c('ENTREZ', 'GOSLIM_ID', 'GOSLIM_TERM'),
                   'GOSLIM_ID')
  p.all <- p.all[!is.na(p.all[['ENTREZ']]),,drop=FALSE]
  p.all <- p.all[order(p.all[['ENTREZ']]),]

  go.all <- aselect(
    GO.db::GO.db,
    unique(p.all[['GOSLIM_ID']]),
    c('GOID', 'TERM'),
    'GOID')
  go <- go.all[!is.na(go.all[['TERM']]),,drop=FALSE]
  go <- go[!duplicated(go$GOID),]

  GO <- merge(p.all, go, by.x = 'GOSLIM_ID', by.y = 'GOID', all.x = TRUE)
  names(GO) <- c("gs_id", "feature_id", "ontology", "name")
  GO <- GO[!is.na(GO$name),]

  if (FALSE) {
    u <- unique(as.data.table(GO), by = c("ontology", "name", "gs_id"))
    sum(duplicated(u$name))
    u[name %in% u$name[duplicated(u$name)]]
  }

  GO$collection <- sprintf("GOSLIM_%s", GO$ontology)
  GO$ontology <- NULL

  gdb <- GeneSetDb(GO)

  for (coll in unique(GO$collection)) {
    geneSetCollectionURLfunction(gdb, coll) <- ".geneSetURL.GOSLIM"
    featureIdType(gdb, coll) <- EntrezIdentifier()
  }

  gdb
}

.geneSetURL.GOSLIM <- function(coll, gsname, gdb = NULL, ...) {
  if (!is(gdb, "GeneSetDb")) {
    return("http://www.pantherdb.org/panther/goSlim.jsp")
  }
  info <- gdb@table[collection == coll & name == gsname]
  if (nrow(info) != 1) {
    "http://www.pantherdb.org/panther/goSlim.jsp"
  } else {
    sprintf("http://amigo.geneontology.org/amigo/term/%s", info$gs_id)
  }
}
