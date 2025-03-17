# The GeneSetDb@collectionMetadata url_functions store everything in their
# environemnts. This can take a lot of runtime memory, but my biggest issue is
# that serializing and loading a GeneSetDb with url_functions is very slow and
# cretest huge objects

test_that("url_functions with .fn.local.vars defined trim their environments", {
  with.url_function.fn <- tempfile("GeneSetDb.with.uf", fileext = ".rds")
  no.url_function.fn <- tempfile("GeneSetDb.no.uf", fileext = ".rds")

  gdb <- getMSigGeneSetDb("C5")
  gdb2 <- gdb
  gdb2@collectionMetadata <-
    data.table::copy(gdb2@collectionMetadata)[name != "url_function"]

  saveRDS(gdb, with.url_function.fn)
  saveRDS(gdb2, no.url_function.fn)

  # size difference not so big
  size.with.fn <- file.info(with.url_function.fn)[["size"]]
  size.no.fn <- file.info(no.url_function.fn)[["size"]]

  # size with function is a bit bigger
  expect_gt(size.with.fn, size.no.fn)

  # but not too big -- approx size in MB's
  expect_lt(size.with.fn / 1024^2 - size.no.fn / 1024^2, 0.2)

  restored <- readRDS(with.url_function.fn)
  # As of sparrow >= 1.13.9 & msigdbr >= 10, we are uising the AMIGO URLs for
  # GO pathways
  gs <- geneSet(gdb, name = "GOBP_2FE_2S_CLUSTER_ASSEMBLY")
  gs.info <- subset(geneSets(gdb), name == gs$name[1])
  
  expect_equal(
    geneSetURL(restored, gs$collection[1], gs$name[1]),
    gs.info$geneset_url,
    check.attributes = FALSE)

  unlink(with.url_function.fn)
  unlink(with.url_function.fn)
})
