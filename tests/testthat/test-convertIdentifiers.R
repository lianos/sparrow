test_that("1:1 convertIdentifiers works as expected", {
  proto <- .GeneSetDb()
  gdb <- exampleGeneSetDb() # this has no symbols in it

  # 1:1 id conversion
  xref <- data.frame(
    current_id = featureIds(gdb),
    new_id = paste0(featureIds(gdb), "_symbol"))
  gdb2 <- convertIdentifiers(gdb, from = "current_id", to = "new_id",
                             xref = xref, extra.cols = "original_id")
  edb <- transform(gdb@db, feature_id = paste0(feature_id, "_symbol"),
                   original_id = feature_id)
  setkeyv(edb, key(proto@db))
  expect_equal(gdb2@db, edb)
})

test_that("1:many and 0:1 convertIdentifiers works as expected", {
  proto <- .GeneSetDb()
  gdb <- exampleGeneSetDb() # this has no symbols in it

  # 1:1 id conversion
  xref <- data.table(
    current_id = featureIds(gdb),
    new_id = paste0(featureIds(gdb), "_symbol"))

  # find an identifier that appears in more than one geneset to duplicate it
  # and add it to the xref table
  nobs <- gdb@db[, list(n = .N), by = feature_id]
  use.me <- nobs[n > 2][1,]
  if (is.na(use.me$feature_id)) {
    use.me <- nobs[1,]
  }

  # This is the base case for the converted db we expect
  edb <- transform(
    gdb@db, feature_id = paste0(feature_id, "_symbol"),
    original_id = feature_id)

  # xref table with 1:many mapping ---------------------------------------------
  dupe <- data.table(
    current_id = use.me$feature_id[1],
    new_id = "additional")
  xref.dupe <- rbind(dupe, xref)
  res.dupe <- convertIdentifiers(gdb, from = "current_id", to = "new_id",
                                 xref = xref.dupe, extra.cols = "original_id")
  # this is what we expect
  has.dupe <- edb[original_id == dupe$current_id]
  has.dupe[, feature_id := dupe$new_id]
  edb.dupe <- rbind(edb, has.dupe)
  setkeyv(edb.dupe, key(proto@db))
  expect_equal(res.dupe@db, edb.dupe)

  # xref table with dropped identifier (0:many) --------------------------------
  xref.drop <- xref[!current_id %in% use.me$feature_id]
  res.drop <- convertIdentifiers(gdb, from = "current_id", to = "new_id",
                                 xref = xref.drop, extra.cols = "original_id")
  edb.drop <- edb[!original_id %in% use.me$feature_id]
  setkeyv(edb.drop, key(proto@db))
  expect_equal(res.drop@db, edb.drop)
})

test_that("1:0 convertIdentifiers drops features as expected", {
  proto <- .GeneSetDb()
  gdb <- exampleGeneSetDb() # this has no symbols in it

  # 1:1 id conversion
  xref <- data.frame(
    current_id = featureIds(gdb),
    new_id = paste0(featureIds(gdb), "_symbol"))

  # find an identifier that appears in more than one geneset to duplicate it
  # and add it to the xref table
  nobs <- gdb@db[, list(n = .N), by = feature_id]
  dupe.me <- nobs[n > 2][1,]
  if (is.na(dupe.me$feature_id)) {
    dupe.me <- nobs[1,]
  }
  dupe <- data.frame(
    current_id = dupe.me$feature_id[1],
    new_id = "additional")
  xref <- rbind(dupe, xref)

  gdb2 <- convertIdentifiers(gdb, from = "current_id", to = "new_id",
                             xref = xref, extra.cols = "original_id")
  edb <- transform(
    gdb@db, feature_id = paste0(feature_id, "_symbol"),
    original_id = feature_id)

  has.dupe <- edb[original_id == dupe$current_id][, feature_id := dupe$new_id]
  edb <- rbind(edb, has.dupe)
  setkeyv(edb, key(proto@db))
  expect_equal(gdb2@db, edb)
})

# Species Conversion -----------------------------------------------------------
# TODO: Exercise species conversion once implemented, with some of these
#       examples.
# # Assume gdb is a GeneSetDb of human entrz id's
# gdb <- exampleGeneSetDb() # this has no symbols in it
#
# # 1. Convert to human ensembl id's
# gdb.hens <- convertIdentifiers(gdb, "human", id.type = "ensembl",
#                                extra.cols = c("symbol", "original"))
# # 2. Convert to mouse ensembl
# gdb.mens <- convertIdentifiers(gdb, from = "human", to = "mouse",
#                                id.type = "ensembl")
# # 3. Convert to mouse entrez
# gdb.ment <- convertIdentifiers(gdb, from = "human", to = "mouse",
#                                id.type = "entrez")
# # 4. You provide your own idenifiers.
# xref <- data.frame(
#   current_id = featureIds(gdb),
#   new_id = paste0(featureIds(gdb), "_symbol"))
# gdb2 <- convertIdentifiers(gdb, from = "current_id", to = "new_id",
#                            xref = xref, extra.cols = "original")
# geneSet(gdb2, name = "BIOCARTA_AGPCR_PATHWAY")

test_that("entrez id's remapped to ensemble", {
  expect_error(convertIdentifiers(), "babelgene")
  # xref <- load_id_xref("human")
  # gdb.entrez <- exampleGeneSetDb()
  # gdb.ens <- remap_identifiers(gdb.entrez, xref,
  #                             original_id = "entrezgene_id",
  #                             target_id = "ensembl_gene_id")
  #
  # gs.entrez <- geneSet(gdb.entrez, name = "REACTOME_RAF_MAP_KINASE_CASCADE")
  # gs.ens <- geneSet(gdb.ens, name = "REACTOME_RAF_MAP_KINASE_CASCADE")
  # expect_true(all(substr(gs.ens$feature_id, 1, 4) == "ENSG"))
  # expect_set_equal(gs.ens$original_id, gs.entrez$feature_id)
})
