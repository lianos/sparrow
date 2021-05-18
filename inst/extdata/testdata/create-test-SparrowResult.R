library(sparrow)

vm <- exampleExpressionSet(do.voom = TRUE)
gsl <- exampleGeneSets()
gsd <- conform(GeneSetDb(gsl), vm)
x <- seas(vm, gsd, c("camera", "fry"), design = vm$design)
path <- "inst/extdata/testdata/test-SparrowResult.rds"
saveRDS(x, path)
