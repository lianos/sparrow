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
  gs.kegg <- subset(geneSets(with.kegg), subcollection == "CP:KEGG")
  expect_true(nrow(gs.kegg) > 0L)

  gs.nokegg <- subset(geneSets(no.kegg), subcollection == "CP:KEGG")
  expect_true(nrow(gs.nokegg) == 0L)
})

test_that("url function stored correctly", {
  go.bp.df <- sparrow:::.pkgcache$msigdb$`Homo sapiens`[gs_subcat == "GO:BP"]
  go.mf.df <- sparrow:::.pkgcache$msigdb$`Homo sapiens`[gs_subcat == "GO:MF"]
  go.cc.df <- sparrow:::.pkgcache$msigdb$`Homo sapiens`[gs_subcat == "GO:CC"]

  gdb.pro <- getMSigGeneSetDb("C5", promote.subcollection = TRUE)
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

# Tests to address functionality to support updated msigdbr / msigdbdf backend
# for newer MsigDB gene sets

test_that("subcategory prefixes are cleaned (or not)", {
  gs.clean <- getMSigGeneSetDb(
    species = "human",
    strip.subcollection.prefix = TRUE,
    refetch = TRUE)
    
  subcols <- geneSets(gs.clean) |> 
    filter(subcollection != "") |> 
    count(subcollection) |> 
    mutate(with_prefix = grepl(":", subcollection))
  # we should have some subcollections
  expect_gt(nrow(subcols), 5)
  
  # The only ones w/ a prefix should be GO:
  prefixed <- subcols |> 
    filter(with_prefix) |> 
    pull(subcollection)
  expect_true(all(startsWith(prefixed, "GO:")))
  
  gs.prefix <- getMSigGeneSetDb(
    species = "human",
    strip.subcollection.prefix = FALSE,
    refetch = TRUE)
  subcols.pre <- geneSets(gs.prefix) |> 
    filter(subcollection != "") |> 
    count(subcollection) |> 
    mutate(with_prefix = grepl(":", subcollection))
  expect_equal(nrow(subcols.pre), nrow(subcols))  
  
  # there should be more subcollection w/ prefix in the 'notstripped version'
  expect_gt(sum(subcols.pre$with_prefix), sum(subcols$with_prefix))
})

test_that("Mixed mouse/human collection returns mutually exclusive genesets", {
  gdb.hs <- getMSigGeneSetDb(
    species = "human", 
    promote.subcollection = FALSE,
    refetch = TRUE)  
  
  gdb.mm <- getMSigGeneSetDb(
    species = "mouse", 
    promote.subcollection = FALSE,
    refetch = TRUE)
  
  # There should not be duplicate geneset names from different db_species
  gs.hs <- gdb.hs@table[, list(n = .N), by = c("name", "db_species")]
  expect_true(all(gs.hs$n == 1))
  expect_true(all(gs.hs$db_species == "HS"))
  
  gs.mm <- gdb.mm@table[, list(n = .N), by = c("name", "db_species")]
  expect_true(all(gs.mm$n == 1))
  expect_setequal(gs.mm$db_species, c("HS", "MM"))
  
  # now don't count db_species in the uniquenes
  gs.mm.unique <- gs.mm[, list(n = .N), by = "name"]
  expect_true(all(gs.mm.unique$n == 1))
})
  

# Explore relationship between mouse and human geneset collections
# Some geneset names are almost similar gs_subcollections, like
# gs_subcollection == "GTRD" in mouse, but "TFT:GTRD" in human
if (FALSE) {
  hs <- msigdbr::msigdbr("human", db_species = "HS")
  mm <- msigdbr::msigdbr("mouse", db_species = "MM")
  
  shared.mm <- subset(mm, gs_name %in% hs$gs_name)
  shared.hs <- subset(hs, gs_name %in% mm$gs_name)
  
  # There is a gs_subcollection in human called MIR:MIRDB and TFT:GTRD, but
  # in mouse they are called MIRDB and GTRD, respectively
}
