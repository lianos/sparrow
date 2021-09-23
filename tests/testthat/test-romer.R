context("romer")

## TODO: Test why the mixied pvalues are so significant always!!?!

test_that('romer runs equivalently from do.romer vs direct call', {
  y <- exampleExpressionSet(do.voom = FALSE)
  y <- edgeR::estimateDisp(y, y$design)

  gdb <- GeneSetDb(exampleGeneSets())
  gdb <- gdb[c(TRUE, rep(FALSE, nrow(gdb) - 1))]
  gdb <- conform(gdb, y)

  # We have to ensure that the genesets are tested in the same order as they
  # are tested from the GeneSetDb for the pvalues to be equivalent given
  # the same random seed.
  global.idxs <- as.list(gdb, value='x.idx')

  nrot <- 250
  seed <- 123

  # run romer directly ---------------------------------------------------------
  set.seed(seed)
  expected <- limma::romer(y, global.idxs, y$design, ncol(y$design), nrot = nrot)

  # run do.romer ---------------------------------------------------------------
  do <- sparrow:::do.romer(gdb, y, y$design, ncol(y$design), nrot = nrot,
                           .random.seed = seed)

  # run romer through seas() ---------------------------------------------------
  mg <- seas(y, gdb, "romer", design = y$design, contrast = ncol(y$design),
             nrot = nrot, .random.seed = seed)
  res <- result(mg, "romer")
  res$key <- encode_gskey(res)

  # do.romer matches limma::romer ----------------------------------------------
  # Test that inernal call matches direct limma call
  expect_true(setequal(rownames(do), rownames(expected)))
  expected <- expected[rownames(do),,drop=FALSE]
  expect_equal(do, expected, check.attributes=FALSE)

  # order of geneset should be the same as the GeneSetDb
  expect_equal(rownames(do), encode_gskey(geneSets(gdb)))
  expect_equal(do[, 'NGenes'], geneSets(gdb)$n, check.attributes=FALSE)

  # seas call matches limma::romer ---------------------------------------------
  expect_equal(res$key, rownames(expected))
  expect_equal(res$n, expected[, 'NGenes'], check.attributes = FALSE)
  expect_equal(res$pval.up, expected[, 'Up'], check.attributes = FALSE)
  expect_equal(res$pval.down, expected[, 'Down'], check.attributes = FALSE)
  expect_equal(res$pval, expected[, 'Mixed'], check.attributes = FALSE)
})
