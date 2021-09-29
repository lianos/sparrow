test_that("Renaming collections maintains dataset integreity", {
  check <- data.frame(
    old_collection = c("c2", "c7"),
    new_collection = c("MSigDB C2", "ImmuneSigDb"),
    geneset = c("BIOCARTA_AGPCR_PATHWAY", "GSE3982_BCELL_VS_TH2_DN"))
  unchanged <- list(collection = "c6", geneset = "CAMP_UP.V1_DN")

  rename <- setNames(check$new_collection, check$old_collection)

  gdb <- exampleGeneSetDb()
  ngdb <- renameCollections(gdb, rename)

  for (i in seq_len(nrow(check))) {
    params <- check[i,]
    url.orig <- geneSetURL(gdb, params$old_collection, params$geneset)
    url.new <- geneSetURL(ngdb, params$new_collection, params$geneset)
    expect_string(url.new, min.chars = 5)
    expect_equal(unname(url.new), unname(url.orig),
                 info = paste("Iteration:", i))
    expect_equal(names(url.orig), params$old_collection,
                 info = paste("Iteration:", i))
    expect_equal(names(url.new), params$new_collection,
                 info = paste("Iteration:", i))
  }

  # This collection name should not be changed
  expected <- geneSetURL(gdb, unchanged$collection, unchanged$geneset)
  expect_string(expected, min.chars = 5)
  check.same <- geneSetURL(ngdb, unchanged$collection, unchanged$geneset)
  expect_string(check.same, min.chars = 5)
  expect_equal(check.same, expected)
})
