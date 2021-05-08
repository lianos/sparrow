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
  expect_equal(
    geneSetURL(restored, "C5", "GOBP_2FE_2S_CLUSTER_ASSEMBLY"),
    c("C5" = "http://www.broadinstitute.org/gsea/msigdb/cards/GOBP_2FE_2S_CLUSTER_ASSEMBLY.html"))

  unlink(with.url_function.fn)
  unlink(with.url_function.fn)
})


