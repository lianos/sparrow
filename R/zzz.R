
# This looks like the way tiledbr packages does it, thanks Dirk!
.pkgcache <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  ## Setup default option values
  opts <- options()

  species <- msigdbr::msigdbr_species()[["species_name"]]
  .pkgcache[["msigdb"]] <- sapply(species, function(x) NULL, simplify = FALSE)
  invisible()
}
