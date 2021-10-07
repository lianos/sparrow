context("retrieve MSigDB collections")

test_that("MSigDB retrieval respects collection subsets", {
  gdb.all <- getMSigGeneSetDb()
  expect_setequal(geneSets(gdb.all)$collection, c("H", paste0("C", 1:8)))
  gdb.sub <- getMSigGeneSetDb(c("H", "C6"))
  expect_setequal(geneSets(gdb.sub)$collection, c("H", "C6"))
})

test_that("with.kegg honors inclusion/exclusion of KEGG gene sets", {
  with.kegg <- getMSigGeneSetDb("c2", with.kegg = TRUE)
  no.kegg <- getMSigGeneSetDb("c2", with.kegg = FALSE)
  gs.kegg <- subset(geneSets(with.kegg), subcategory == "CP:KEGG")
  expect_true(nrow(gs.kegg) > 0L)

  gs.nokegg <- subset(geneSets(no.kegg), subcategory == "CP:KEGG")
  expect_true(nrow(gs.nokegg) == 0L)
})

test_that("url function stored correctly", {
  go.bp.df <- sparrow:::.pkgcache$msigdb$`Homo sapiens`[gs_subcat == "GO:BP"]
  go.mf.df <- sparrow:::.pkgcache$msigdb$`Homo sapiens`[gs_subcat == "GO:MF"]
  go.cc.df <- sparrow:::.pkgcache$msigdb$`Homo sapiens`[gs_subcat == "GO:CC"]

  gdb.pro <- getMSigGeneSetDb("C5", promote.subcategory.to.collection = TRUE)
  gdb.npro <- getMSigGeneSetDb("C5")

  genesets <- c(
    BP = "GOBP_LIVER_REGENERATION",
    MF = "GOMF_ENZYME_ACTIVATOR_ACTIVITY",
    CC = "GOCC_GOLGI_APPARATUS")

  base.url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"

  for (gocat in names(genesets)) {
    goname <- genesets[[gocat]]
    expected.url <- sprintf(base.url, goname)
    pro.url <- geneSetURL(gdb.pro, sprintf("C5_GO:%s", gocat), goname)
    expect_equal(unname(pro.url), expected.url,
                 info = paste("promoted subcat url", gocat, goname, sep = ":"))
    npro.url <- geneSetURL(gdb.npro, "C5", goname)
    expect_equal(unname(pro.url), expected.url,
                 info = paste("no promo subcat url", gocat, goname, sep = ":"))
  }
})
