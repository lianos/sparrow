context("overrepresentation analysis ('ora')")

# Manually slimming down the GeneSetDb to only include features in the
# 1000 gene data.frame input test data
gdb. <- local({
  dge.res <- exampleDgeResult("human", "ensembl")
  randomGeneSetDb(dge.res)
})

test_that("'naked' ora call vs seas pipeline are equivalent", {
  dfinput <- exampleDgeResult("human", "ensembl",
                              induce.bias = "effective_length")
  nres <- ora(dfinput, gdb., selected = "selected", groups = "direction",
              feature.bias = "effective_length")
  mres <- seas(setNames(dfinput$t, dfinput$feature_id), gdb., "ora",
               feature.bias = "effective_length",
               xmeta. = dfinput)

  groups <- c(all = "ora", up = "ora.up", down = "ora.down")
  expect_setequal(resultNames(mres), groups)

  for (i in seq(groups)) {
    # Call ora direct
    ename <- names(groups)[i]
    pcol <- paste0("P.", ename)
    cmp <- nres[, c("Pathway", "N", ename, pcol)]

    # Pull out of SparrowResult object
    mname <- groups[i]
    mg <- result(mres, mname)
    mg <- mg[, c("collection", "name", "N", "n", "n.drawn", "pval")]

    expect_equal(mg$name, sub(".*;;", "", nres$Pathway), info = ename)
    expect_equal(mg$n, cmp$N, info = ename)
    expect_equal(mg$n.drawn, cmp[[ename]], info = ename)
    expect_equal(mg$pval, cmp[[pcol]], info = ename)
  }
})

test_that("induced length associattion to significance is accounted for", {
  set.seed(0xBEEF)
  biased <- exampleDgeResult("human", "ensembl",
                             induce.bias = "effective_length")

  nbias <- ora(biased, gdb., selected = "selected", feature.bias = NULL)
  lbias <- ora(biased, gdb., selected = "selected",
               feature.bias = "effective_length")

  # general returned structure is the same
  expect_equal(
    nbias[, c("Pathway", "N", "all")],
    lbias[, c("Pathway", "N", "all")])

  # when accounting for length bias, majority of pvalues should be penalized
  frac.less <- mean(lbias$P.all > nbias$P.all)
  expect_true(frac.less > 0.7)

  if (FALSE) {
    plot_ora_bias(biased, "selected", "effective_length")
    plot(-log10(nbias$P.all), -log10(lbias$P.all))
    abline(0,1,col = "red")
  }

  # When we adjust for a randomized length bias, the results should be similar
  # to those that do not adjust for a length bias at all
  set.seed(0xFEED)
  rando <- transform(biased, effective_length = sample(effective_length))
  rbias <- ora(rando, gdb., selected = "selected",
               feature.bias = "effective_length")
  expect_lt(mean(rbias$P.all - nbias$P.all)**2, 1e-4)

  if (FALSE) {
    plot_ora_bias(rando, "selected", "effective_length")
  }
})

test_that("ora,groups variable accepts column or list", {
  dfinput <- exampleDgeResult("human", "ensembl")
  group.list <- split(dfinput$feature_id, dfinput$direction)
  p.cols <- paste0("P.", c("all", "down", "up"))

  g1 <- ora(dfinput, gdb., selected = "selected", groups = "direction")
  g2 <- ora(dfinput, gdb., selected = "selected", groups = group.list)
  expect_equal(g1$Pathway, g2$Pathway)
  for (pname in p.cols) {
    expect_numeric(g1[[pname]], info = pname)
    expect_equal(g1[[pname]], g2[[pname]], info = pname)
  }

  # This is some metadata that is required for correctly processing these
  # results inside the 'seas' pipeline
  expect_true(attr(g1, "mgunlist"))
  expect_setequal(attr(g1, "groups"), c("all", names(group.list)))
  expect_true(attr(g1, "rawresult"))
})

test_that("ora and goseq give probably approximately correct answers", {
  dfinput <- exampleDgeResult("human", "ensembl",
                              induce.bias = "effective_length")

  # no bias correction
  e1 <- ora(dfinput, gdb., selected = "selected", groups = "direction")
  g1 <- expect_warning({
    sparrow::goseq(
      gdb.,
      subset(dfinput, selected)$feature_id,
      dfinput$feature_id,
      setNames(dfinput$effective_length, dfinput$feature_id),
      method = "Hypergeometric")
  }, "initial point")

  expect_equal(e1$Pathway, g1$category)
  expect_equal(e1$P.all, g1$over_represented_pvalue)

  # length correction
  e2 <- ora(dfinput, gdb., selected = "selected", groups = "direction",
            feature.bias = "effective_length")
  g2 <- expect_warning({
    sparrow::goseq(
      gdb.,
      subset(dfinput, selected)$feature_id,
      dfinput$feature_id,
      setNames(dfinput$effective_length, dfinput$feature_id),
      method = "Wallenius")
  }, "initial point")

  if (FALSE) {
    par(mfrow = c(1, 2))
    plot(-log10(e1$P.all), -log10(e2$P.all),
         main = "Uncorrected vs corrected ora",
         xlab = "Uncorrected", ylab = "Corrected")
    abline(0, 1, col = "red")
    plot(-log10(g1$over_represented_pvalue), -log10(g2$over_represented_pvalue),
         main = "Uncorrected vs corrected goseq",
         xlab = "Uncorrected", ylab = "Corrected")
    abline(0, 1, col = "red")

    par(mfrow = c(1,1))

    # ora is a bit more conservative
    plot(-log10(e2$P.all), -log10(g2$over_represented_pvalue),
         main = "Corrected ora vs goseq",
         xlab = "ora", ylab = "goseq")
    abline(0, 1, col = "red")
  }
  # test that average difference is less than a threshold
  pval.diff <- abs(e2$P.all - g2$over_represented_pvalue)
  expect_lt(mean(pval.diff), 0.025)
})

test_that("ora over ANOVA anaysis works through seas", {
  y <- exampleExpressionSet('tumor-subtype', do.voom=FALSE)
  di <- model.matrix(~ PAM50subtype, data = y$samples)
  vm <- voom(y, di)
  gdb <- exampleGeneSetDb()
  mg <- seas(vm, gdb, "ora", design = vm$design, contrast = 2:3)
  r <- result(mg)
  expect_numeric(r[["pval"]])
  expect_true(sum(r[["pval"]] < 0.002) > 0)
})
