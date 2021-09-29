context("data.frame ranks/DGE input to seas")

xdf <- exampleDgeResult()
scores <- setNames(xdf$logFC, xdf$feature_id)

gdb <- randomGeneSetDb(xdf)
gdb <- conform(gdb, xdf$feature_id)
gs.idx <- as.list(gdb, active.only = TRUE, value = "x.idx")

test_that("ranks-based gsea works", {
  mg <- seas(scores, gdb, "cameraPR", xmeta. = xdf)
  mgres <- result(mg)
  mgres$key <- encode_gskey(mgres)

  cpr <- limma::cameraPR(scores, gs.idx, sort = FALSE)

  expect_equal(mgres$key, rownames(cpr))
  expect_equal(mgres$n, cpr$NGenes)
  expect_equal(mgres$Direction, cpr$Direction)
  expect_equal(mgres$pval, cpr$PValue)
  expect_equal(mgres$padj, cpr$FDR)
})

test_that("data.frame input is same as ranked vector input", {
  # This follows up from previous test
  mgv <- seas(scores, gdb, "cameraPR", xmeta. = xdf)
  mgdf <- seas(xdf, gdb, "cameraPR", rank_by = "logFC", rank_order = "ordered")
  expect_equal(result(mgdf), result(mgv))
  expect_equal(logFC(mgdf), logFC(mgdf))
})

test_that("enrichment-based methods work", {
  fbias <- setNames(xdf$effective_length, xdf$feature_id)
  mg <- seas(scores, gdb, "goseq", xmeta. = xdf, feature.bias = fbias)
  mgres <- result(mg, "goseq")
  mgres$key <- encode_gskey(mgres)

  gseq <- expect_warning({
    sparrow::goseq(
      gdb,
      selected = subset(xdf, significant)$feature_id,
      universe = xdf$feature_id,
      feature.bias = fbias)
  }, "initial point")

  expect_equal(mgres$key, gseq$category)
  expect_equal(mgres$pval, gseq$over_represented_pvalue)
  expect_equal(mgres$pval.under, gseq$under_represented_pvalue)
  expect_equal(mgres$n.drawn, gseq$numDEInCat)
  expect_equal(mgres$n, gseq$numInCat)
})
