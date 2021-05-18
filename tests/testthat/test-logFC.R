context("Pass-through logFC GSEA method")

test_that("logFC pass through generates expected gene-set stats", {
  vm <- exampleExpressionSet()
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)

  mgc <- seas(vm, gsd, "camera", design = vm$design)
  mgf <- seas(vm, gsd, design = vm$design)
  expect_equal(geneSets(mgc), geneSets(mgf))
})

test_that("t-stats and logFCs match full design when only stats passed", {
  vm <- exampleExpressionSet()
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)

  mgf <- suppressWarnings(seas(vm, gsd, desgin = vm$design))
  x <- logFC(mgf)

  tstats <- setNames(x$t, x$feature_id)
  lfc <- setNames(x$logFC, x$feature_id)

  mgt <- seas(tstats, gsd)
  mgl <- seas(lfc, gsd)

  ro <- suppressWarnings(results(mgf))
  rt <- suppressWarnings(results(mgt))
  rl <- suppressWarnings(results(mgl))

  expect_equal(ro$mean.t, rt$mean.t)
  expect_equal(ro$mean.t.trim, rt$mean.t.trim)

  expect_equal(ro$mean.logFC, rl$mean.logFC)
  expect_equal(ro$mean.logFC.trim, rl$mean.logFC.trim)
})
