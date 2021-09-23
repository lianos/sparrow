context("romer")

## TODO: Test why the mixied pvalues are so significant always!!?!

test_that('romer runs equivalently from do.romer vs direct call', {
  y <- exampleExpressionSet(do.voom = FALSE)
  y <- edgeR::estimateDisp(y, y$design)

  gdb <- GeneSetDb(exampleGeneSets())
  gdb <- conform(gdb, y)

  seed <- 123
  nrot <- 250

  # run limma::romer -----------------------------------------------------------
  set.seed(seed)
  expected <- limma::romer(y, as.list(gdb, value = "x.idx"),
                           y$design, ncol(y$design), nrot = nrot)

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
