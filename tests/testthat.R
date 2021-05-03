library("sparrow")
library("testthat")
library("data.table")
library("dplyr")

## Remove temporary files that were generated
test.dir <- system.file('tests', package = "sparrow")
pdfs <- dir(test.dir, "\\.pdf$", full.names=TRUE)
if (length(pdfs)) {
  unlink(pdfs)
}
