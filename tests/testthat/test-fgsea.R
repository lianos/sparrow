context("fgsea")

test_that("seas calculate t and preranked t match fgsea results", {
  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()

  gseaParam <- 1
  nperm <- 1000

  # Since Bioc 3.5, running fgsea warns about ties in preranked stats
  expect_warning({
    mgt <- seas(vm, gdb, 'fgsea', design = vm$design, contrast = 'tumor',
                score.by = 't', nPermSimple = nperm, gseaParam = gseaParam,
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
    rest <- fgsea::fgsea(
      gs.idxs, ranks.t,
      minSize = min.max[1], maxSize = min.max[2],
      nPermSimple = nperm,
      gseaParam = gseaParam)
  }, "ties")

  # these are exactly the same
  expect_equal(nrow(mgres), nrow(rest))
  expect_equal(mgres$pathway, rest$pathway)
  expect_equal(rest$size, mgres$n)
  expect_equal(mgres$ES, rest$ES)
  expect_equal(mgres$leadingEdge, rest$leadingEdge)
  expect_equal(mgres$NES, rest$NES)
  expect_equal(mgres$pval, rest$pval)

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

