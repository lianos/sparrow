context("SparrowResult")

test_that("subset.SparrowResult works", {
  ## TODO: Test subset.SparrowResult
})

test_that("GeneSetDb defined with logFC-like column names are kosher", {
  gdb <- exampleGeneSetDb()
  gdb@db$logFC <- rnorm(nrow(gdb@db))
  es <- exampleExpressionSet()
  mg <- seas(es, gdb, "camera", design = es$design)

  gs <- geneSet(mg, "c2", "BIOCARTA_AGPCR_PATHWAY")
  expect_true(all(c("logFC", "logFC.gs") %in% colnames(gs)))

  info <- merge(gs, logFC(mg), all.x = TRUE, by = "feature_id")

  expect_equal(nrow(gs), nrow(info))
  expect_true(sum(is.na(info$logFC.y)) == 0)
  expect_equal(info$logFC.x, info$logFC.y)
})
