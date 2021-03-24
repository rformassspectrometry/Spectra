library(testthat)
test_that(".isotope_peaks works", {
  
  ints <- c(1, 2, 1, 1, 1, 2, 1, 0.6, 0)
  x <- cbind(mz = seq_along(ints), intensity = ints)
  isoDef <- cbind(mzd = c(1, 2, 3), maxint = c(0.7, 0.6, 0.4))
  
  res <- .isotope_peaks(x, isotopeDefinition = isoDef)
  expect_equal(res, list(c(2, 3, 4), c(6, 7, 8), c(7, 8)))
  
  res <- .isotope_peaks(x, isotopeDefinition = isoDef, seed_mz = x[c(2, 5), 1])
  expect_equal(res, list(c(2, 3, 4)))
  
  res <- .isotope_peaks(x, isotopeDefinition = isoDef, duplicates = "remove")
  expect_equal(res, list(c(2, 3, 4)))
  
  res <- .isotope_peaks(x, isotopeDefinition = isoDef, duplicates = "closest")
  expect_equal(res, list(c(2, 3, 4), c(6, 7, 8)))
  
  x[6] <- x[6] + 10^(-6)
  res <- .isotope_peaks(x, isotopeDefinition = isoDef)
  expect_equal(res, list(c(2, 3, 4), c(6, 7, 8), c(7, 8)))
  
  res <- .isotope_peaks(x, isotopeDefinition = isoDef, duplicates = "closest")
  expect_equal(res, list(c(2, 3, 4), c(7, 8)))
  
})
