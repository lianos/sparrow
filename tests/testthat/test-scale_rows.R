M <- list(
  wide = matrix(rnorm(50, mean = 1, sd = 2), nrow = 5,
                dimnames = list(LETTERS[1:5], letters[1:10])),
  long = matrix(rnorm(50, mean = 1, sd = 2), nrow = 10,
                dimnames = list(LETTERS[1:10], letters[1:5])),
  square = matrix(rnorm(100, mean = 1, sd = 2), nrow = 10,
                  dimnames = list(LETTERS[1:10], letters[1:10])))

test_that("scale and center are flags", {
  params <- expand.grid(center = c(TRUE, FALSE), scale = c(TRUE,FALSE))

  for (orientation in names(M)) {
    m <- M[[orientation]]
    for (i in seq(nrow(params))) {
      p <- params[i,,drop = FALSE]
      res <- scale_rows(m, center = p$center, scale = p$scale)
      ex <- t(scale(t(m), center = p$center, scale = p$scale))
      info <- paste(names(p), unname(unlist(p)), sep = ":", collapse = ",")
      expect_equal(res, ex, info = sprintf("%s (%s)", info, orientation),
                   check.attributes = FALSE)
    }
  }
})

test_that("centering on a subgroup of specified columns", {
  ctrl <- c("b", "d", "e")
  params <- expand.grid(scale = c(TRUE, FALSE))
  for (orientation in names(M)) {
    m <- M[[orientation]]
    means <- matrixStats::rowMeans2(m, cols = match(ctrl, colnames(m)))
    sds <- matrixStats::rowSds(m, cols = match(ctrl, colnames(m)))
    for (i in seq(nrow(params))) {
      p <- params[i,,drop = FALSE]
      info <- paste(names(p), unname(unlist(p)), sep = ":", collapse = ",")
      info <- sprintf("%s (%s)", info, orientation)

      if (p$scale) {
        wrap <- function(expr) expect_warning(eval(expr), "weird")
      } else {
        wrap <- identity
      }

      res <- wrap({
        scale_rows(m, center = ctrl, scale = p$scale)
      })
      resi <- wrap({
        scale_rows(m, center = match(ctrl, colnames(m)), scale = p$scale)
      })
      resl <- wrap({
        scale_rows(m, center = colnames(m) %in% ctrl, scale = p$scale)
      })
      ex <- t(scale(t(m), center = means, scale = if (p$scale) sds else FALSE))
      expect_equal(res, ex, info = info, check.attributes = FALSE)
      expect_equal(resi, ex, info = paste(info, "[by index]"),
                   check.attributes = FALSE)
      expect_equal(resl, ex, info = paste(info, "[by logical]"),
                   check.attributes = FALSE)
    }
  }
})

test_that("centering on pre-calculated values", {
  max.rows <- max(sapply(M, nrow))
  means.all <- setNames(rnorm(max.rows), sample(head(LETTERS, max.rows)))

  params <- expand.grid(scale = c(TRUE, FALSE))
  for (orientation in names(M)) {
    m <- M[[orientation]]
    means <- means.all[rownames(m)]
    for (i in seq(nrow(params))) {
      p <- params[i,,drop = FALSE]
      info <- paste(names(p), unname(unlist(p)), sep = ":", collapse = ",")
      info <- sprintf("%s (%s)", info, orientation)

      res <- scale_rows(m, center = means, scale = p$scale,
                        base.attributes = TRUE)
      ex <- t(scale(t(m), center = means, scale = p$scale))

      expect_equal(res, ex, info = info, check.attributes = FALSE)
    }
  }
})
