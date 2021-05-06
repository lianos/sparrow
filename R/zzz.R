
# This looks like the way tiledbr packages does it, thanks Dirk!
.pkgcache <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  ## Setup default option values
  opts <- options()

  # Having this in .onLoad requires that babelgene is in Imports not Suggests
  species <- babelgene::species()[["scientific_name"]]
  .pkgcache[["msigdb"]] <- sapply(species, function(x) NULL, simplify = FALSE)
  invisible()
}
