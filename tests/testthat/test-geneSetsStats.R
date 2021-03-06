context("geneSetsStats")

test_that("geneSetsStats", {
  vm <- exampleExpressionSet()
  gsl <- exampleGeneSets()
  gsd <- conform(GeneSetDb(gsl), vm)

  mg <- seas(vm, gsd, 'geneSetTest', design = vm$design, ranks.only = TRUE)

  trim <- 0.10
  min.logFC <- 1
  max.padj <- 0.10
  gs.stats <- geneSetsStats(mg, min.logFC, max.padj, trim, as.dt=TRUE)

  ## calculate expected
  istats <- calculateIndividualLogFC(vm, vm$design, ncol(vm$design), as.dt=TRUE)
  data.table::setkeyv(istats, 'feature_id')

  ## The code below is the same used in do.geneSetScores -- let's think of a
  ## more manual way to do this.
  expected <- geneSets(gsd, as.dt=TRUE)[, {
    ids <- featureIds(gsd, .BY[[1]], .BY[[2]])
    lfc <- istats[list(ids)]
    t.nona <- lfc$t[!is.na(lfc$t)]
    list(mean.logFC=mean(lfc$logFC, na.rm=TRUE),
         mean.logFC.trim=mean(lfc$logFC, na.rm=TRUE, trim=trim),
         mean.t=mean(lfc$t, na.rm=TRUE),
         mean.t.trim=mean(lfc$t, na.rm=TRUE, trim=trim))
  }, by=c('collection', 'name')]

  expect_equal(expected, gs.stats[, names(expected), with=FALSE])
})
