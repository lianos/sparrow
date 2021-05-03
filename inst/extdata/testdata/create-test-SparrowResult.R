library(sparrow)

vm <- exampleExpressionSet(do.voom=TRUE)
gsl <- exampleGeneSets()
gsd <- conform(GeneSetDb(gsl), vm)
x <- seas(gsd, vm, vm$design, methods=c('camera', 'fry'))
path <- 'inst/extdata/testdata/test-SparrowResult.rds'
saveRDS(x, path)
