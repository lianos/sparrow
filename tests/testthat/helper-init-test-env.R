# library("magrittr")
# library("data.table")
# library("dtplyr")
# library("dplyr")

# Since loading data from msigdbr can be a bit slow, we load a serious of
# GeneSetDb objects here that will be used throughout the test harness.

# Because we started life in this package way back when, entrez is the default
# gene identifier type returned here
# kGdb <- list(
#   all <- getMSigGeneSetDb("H", "human", "entrez"),
#   hallmark <- getMSigGeneSetDb("H", "human", "entrez"),
#   c6 <- getMSigGeneSetDb("C6", "human", "entrez"),
#   hallmark.ens <- getMSigGeneSetDb("H", "human", "ensembl"))

