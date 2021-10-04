# The dataframe input was created from a differential expression analysis
# done "elsewhere". The actual values don't really matter, we just need a
# standard input data.frame downstream from some DGE to play with.
#
# In the future we will generate this result from a standard dataset available
# as an experiment package or from experiment hub.

# we need to add entrez identifiers to it
set.seed(123)
library(data.table)
orig.dat <- fread("inst/extdata/testdata/dataframe-input-short.csv.gz")

library(org.Hs.eg.db)
orig.dat <- transform(
  orig.dat,
  entrez_id = as.character(
    mapIds(
      org.Hs.eg.db,
      orig.dat$feature_id,
      "ENTREZID",
      "ENSEMBL")))

out <- orig.dat[!is.na(entrez_id)]

write.csv(out, "inst/extdata/testdata/dataframe-input-short.csv", row.names = FALSE)
system("gzip inst/extdata/testdata/dataframe-input-short.csv")
