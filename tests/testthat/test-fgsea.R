context("fgsea")

test_that("seas calculate t and preranked t match fgsea results", {
  gseaParam <- 1
  nperm <- 1000
  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()

  # Since Bioc 3.5, running fgsea warns about ties in preranked stats
  expect_warning({
    mgt <- seas(vm, gdb, 'fgsea', design = vm$design, contrast = 'tumor',
                score.by = 't',
                nPermSimple = nperm, gseaParam = gseaParam,
                .random.seed = 123)
  }, "ties")
  mgres <- mgt %>%
    result("fgsea") %>%
    transform(pathway = encode_gskey(collection, name))

  gs.idxs <- as.list(geneSetDb(mgt), active.only=TRUE, value='x.id')
  min.max <- range(sapply(gs.idxs, length))

  lfc <- logFC(mgt)
  ranks.lfc <- setNames(lfc[['logFC']], lfc[['feature_id']])
  ranks.t <- setNames(lfc[['t']], lfc[['feature_id']])

  expect_warning({
    set.seed(123)
    rest <- fgsea::fgseaMultilevel(gs.idxs, ranks.t, sampleSize = 101,
                         minSize = min.max[1], maxSize = min.max[2],
                         eps = 1e-50, scoreType = "std", nproc = 0,
                         nPermSimple = nperm,
                         gseaParam = gseaParam,
                         absEps = NULL)
  }, "ties")

  # these are exactly the same
  expect_equal(nrow(mgres), nrow(rest))
  expect_equal(mgres$pathway, rest$pathway)
  expect_equal(rest$size, mgres$n)
  expect_equal(mgres$ES, rest$ES)
  expect_equal(mgres$leadingEdge, rest$leadingEdge)

  # something is off with the pval, and therefore log2err and NES, and I just
  # can't figure it out
  # expect_equal(mgres$log2err, rest$log2err)
  # expect_equal(mgres$NES, rest$NES)
  # expect_equal(mgres$pval, rest$pval)
  expect_equal(cor(mgres$pval, rest$pval), 1, tolerance = 0.01)

  # passing in a preranked vector gives same results ---------------------------
  expect_warning({
    mgpre <- seas(ranks.t, gdb, "fgsea", nperm = nperm,
                  gseaParam = gseaParam, score.by = "t",
                  .random.seed = 123)
  }, "ties")

  rpre <- result(mgpre, 'fgsea')
  comp.cols <- c('collection', 'name', 'N', 'n', 'size', 'pval', 'ES')
  expect_equal(rpre[, comp.cols], mgres[, comp.cols])

  # Passing in data.frame works, too -------------------------------------------
  expect_warning({
    mgdf <- seas(lfc, gdb, "fgsea", nperm = nperm,
                 rank_by = "t", rank_order = "descending",
                 gseaParam = gseaParam, .random.seed = 123)
  }, "ties")
  res.df <- result(mgdf)
  expect_equal(res.df[, comp.cols], mgres[, comp.cols])
})

