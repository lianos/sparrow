context("geneSetSummaryByGenes")

test_that("geneSetSummaryByGenes,GeneSetDb returns a legit result", {
  set.seed(0xBEEF)
  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()
  features <- sample(featureIds(gdb), 10)

  res <- geneSetSummaryByGenes(gdb, features, with.features = TRUE)

  # 1. Ensure that geneset <-> geneset membership is legit. We do this by
  #    manipulating the result into a long form data.table that looks like
  #    what is stored in a GeneSetDb@db. We then filter the gdb.sub@db
  #    table to only include the queried features, then compare the two.
  gdb.sub <- subsetByFeatures(gdb, features)
  db.expect <- gdb.sub@db %>%
    copy() %>%
    subset(feature_id %in% features) %>%
    setkeyv(c('collection', 'name', 'feature_id'))
  db.result <- res %>%
    dplyr::select(collection, name, starts_with('featureId_')) %>%
    setDT() %>%
    melt(c('collection', 'name')) %>%
    setDF() %>%
    dplyr::rename(feature_id = variable, present=value) %>%
    dplyr::mutate(feature_id = sub('featureId_', '', feature_id)) %>%
    dplyr::filter(present) %>%
    dplyr::select(-present) %>%
    setDT() %>%
    setkeyv(key(db.expect))
  expect_equal(db.result, db.expect)
})

test_that("geneSetSummaryByGenes,SparrowResult returns a legit result", {
  set.seed(0xBEEF)
  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()
  mg <- seas(vm, gdb, "camera", design = vm$design, contrast = ncol(vm$design))
  mgdb <- geneSetDb(mg)
  features <- sample(featureIds(mg), 10)

  res <- geneSetSummaryByGenes(mg, features, with.features=TRUE,
                               feature.rename=FALSE)

  ## Check that logFC for each feature is accurate
  lfc <- res %>%
    dplyr::select({{features}}) %>%
    setDT() %>%
    melt(measure.vars = features, variable.factor = FALSE) %>%
    setDF() %>%
    dplyr::rename(feature_id = variable, logFC = value) %>%
    dplyr::filter(logFC != 0) %>%
    unique(by = "feature_id") %>%
    dplyr::arrange(feature_id)

  lfc.ex <- logFC(mg) %>%
    dplyr::select(feature_id, logFC) %>%
    dplyr::filter(feature_id %in% features) %>%
    dplyr::arrange(feature_id)
  expect_equal(lfc, lfc.ex, check.attributes = FALSE)

  ## check that symbol remapping works, too
  res.s <- geneSetSummaryByGenes(mg, features, with.features=TRUE,
                                 feature.rename='symbol')

  lfc.ex <- logFC(mg) %>%
    dplyr::filter(feature_id %in% features) %>%
    dplyr::transmute(
      renamed = ifelse(
        !is.na(symbol),
        symbol, paste0('featureId_', feature_id)),
      logFC) %>%
    dplyr::arrange(renamed)
  expect_true(all(lfc.ex$renamed %in% colnames(res.s)))

  lfc.s <- res.s %>%
    dplyr::select(!!lfc.ex$renamed) %>%
    setDT() %>%
    melt(measure.vars = colnames(.), variable.factor = FALSE) %>%
    setDF() %>%
    dplyr::rename(renamed = variable, logFC = value) %>%
    dplyr::filter(logFC != 0) %>%
    unique(by = "renamed") %>%
    dplyr::arrange(renamed)
  expect_equal(lfc.s, lfc.ex, check.attributes = FALSE)
})

test_that("geneSetSummary,SparrowResult properly filters significant genesets", {
  set.seed(0xBEEF)
  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()
  mg <- seas(vm, gdb, "camera", design = vm$design, contrast = ncol(vm$design))
  p.thresh <- 0.20
  camera.sig <- dplyr::filter(result(mg, 'camera'), padj <= p.thresh)

  features <- sample(featureIds(mg), 10)
  res.all <- geneSetSummaryByGenes(mg, features, with.features=TRUE,
                                   feature.rename='symbol', as.dt=TRUE)
  res.sig <- geneSetSummaryByGenes(mg, features, with.features=TRUE,
                                   feature.rename='symbol', as.dt=TRUE,
                                   method='camera', max.p=p.thresh)
  expect_true(all(res.sig$name %in% res.all$name))
  expect_true(all(res.sig$name %in% camera.sig$name))
})
