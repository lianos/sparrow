context("Single Sample Gene Set Scoring")

test_that("scoreSingleSamples can use genesets of size n = 1 and value is gene", {
  vm <- exampleExpressionSet()
  genes <- c(GZMA = "3001", PRF1 = "5551", TGFB1 = "7040")
  lol <- list(a = genes[1], b = genes[1:2], c = genes[1:3])
  gdbx <- GeneSetDb(lol, collectionName = "custom")
  scores <- scoreSingleSamples(gdbx, vm, min.gs.size = 1, as.matrix = TRUE)
  expect_equal(scores["custom;;a",], vm$E[genes[1],])
})

test_that("do.scoreSingleSamples.gsva produces correct gsva,plage,ssGSEA scores", {
  vm <- exampleExpressionSet()
  gdb <- conform(exampleGeneSetDb(), vm)
  lol <- as.list(gdb)

  E <- vm$E
  
  gparams <- list(
    gsva = GSVA::gsvaParam(vm$E, lol),
    plage = GSVA::plageParam(vm$E, lol),
    ssgsea = GSVA::ssgseaParam(vm$E, lol))
  
  for (method in names(gparams)) {
    ex <- GSVA::gsva(gparams[[method]])
    res <- scoreSingleSamples(gdb, vm, methods = method, as.matrix = TRUE)
    expect_equal(res, ex, info = paste0("GSVA::", method), check.attributes = FALSE)
  }
})

test_that("multiple 'melted' scores are returned in a long data.frame", {
  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()
  scores <- scoreSingleSamples(gdb, vm$E, methods=c('svd', 'ssgsea'))
  expect_is(scores, 'data.frame')
  expect_true(setequal(c('svd', 'ssgsea'), scores$method))
  n.samples <- ncol(vm)
  n.gs <- nrow(geneSets(gdb))
  expect_true(nrow(scores) == n.samples * n.gs * 2)
})

test_that("ssGSEA.normalize returns same normalization as GSVA", {
  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()

  scores <- scoreSingleSamples(gdb, vm$E, methods='ssgsea', parallel.sz=1,
                               verbose=FALSE)
  my.norm <- sparrow:::ssGSEA.normalize(scores$score)
  ssgsea.norm <- scoreSingleSamples(gdb, vm$E, methods='ssgsea', parallel.sz=1,
                                    ssgsea.norm=TRUE)
  expect_equal(my.norm, ssgsea.norm$score)
})

test_that("ssGSEA (raw) scores are not affected by samples included in test", {
  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()
  some <- sample(ncol(vm), 10)
  scores.all <- scoreSingleSamples(gdb, vm$E, methods='ssgsea', parallel.sz=1,
                                   verbose=FALSE)
  scores.some <- scoreSingleSamples(gdb, vm$E[, some],
                                    methods='ssgsea', parallel.sz=1,
                                    verbose=FALSE)
  scores <- merge(scores.all, scores.some, suffixes=c('.all', '.some'),
                  by=c('collection', 'name', 'sample_id'))
  expect_equal(scores$scores.all, scores$scores.some)
})

test_that("eigenWeightedMean with equal weights can be same as zScore", {
  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()
  gdbc <- conform(gdb, vm)
  gs.idxs <- as.list(gdbc, value='x.idx')

  E <- vm$E[gs.idxs[[1]], ]
  ewm.pc1 <- eigenWeightedMean(E, unscale=FALSE, uncenter=FALSE)
  ewm.z <- eigenWeightedMean(E, weights=1, unscale=FALSE, uncenter=FALSE)
  zs <- zScore(E)
  expect_equal(ewm.z$score, zs$score)
  expect_is(all.equal(ewm.pc1$score, zs$score), 'character')
})

test_that("normalization works in eigenWeightedMean", {
  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()
  unorm <- scoreSingleSamples(gdb, vm, method = "ewm", normalize = FALSE)
  norm <- scoreSingleSamples(gdb, vm, method = "ewm", normalize = TRUE)
  expect_true(all(unorm$score >= norm$score))
  expect_true(all(unorm$score >= 0))
  expect_true(any(norm$score < 0))
})

test_that("eigenWeightedMean can handle 0sd features", {
  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()
  # Added to address Issue #20
  # https://github.com/lianos/multiGSEA/issues/20
  gs <- geneSet(conform(gdb, vm), name = "GOZGIT_ESR1_TARGETS_DN")
  E.o <- vm$E[gs$feature_id,]

  # 0 out low variance genes: these will provide minor contributions to the
  # geneset score anyway
  rvars <- DelayedMatrixStats::rowVars(E.o)
  rvorder <- order(rvars)
  nuke.n <- 3
  zero.idxs <- head(rvorder, nuke.n)

  E.0sd <- E.o
  # I picked 10 because its contribution to PC1 score is > 0.2
  E.0sd[zero.idxs,] <- 0

  expected.score <- eigenWeightedMean(E.o)
  test.score <- expect_warning({
    eigenWeightedMean(E.0sd, .add_noise = TRUE)
  }, paste("Found NaN.*", nuke.n, "features with 0-sd"))

  # eigenWeightedMean scoring should be robust to random noise
  cors <- cor(expected.score$score, test.score$score, method = "spearman")
  expect_true(cors > 0.99)
})
