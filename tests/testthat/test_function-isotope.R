library(testthat)
test_that(".isotope_peaks works", {
  
  ints <- c(1, 2, 1, 1, 0, 1, 2, 1, 0.6)
  x <- cbind(mz = seq_along(ints), intensity = ints)
  isoDef <- cbind(mzd = c(1, 2, 3), maxint = c(0.7, 0.6, 0.4))
  
  res <- .isotope_peaks(x, isotopeDefinition = isoDef)
  expect_equal(res, list(c(2, 3, 4), c(7, 8, 9), c(8, 9)))
  
  res <- .isotope_peaks(x, isotopeDefinition = isoDef, seed_mz = x[c(2, 6), 1])
  expect_equal(res, list(c(2, 3, 4)))
  
  res <- .isotope_peaks(x, isotopeDefinition = isoDef, duplicates = "remove")
  expect_equal(res, list(c(2, 3, 4)))
  
  res <- .isotope_peaks(x, isotopeDefinition = isoDef, duplicates = "closest")
  expect_equal(res, list(c(2, 3, 4), c(7, 8, 9)))
  
  x[7, 1] <- x[7, 1] + 10^(-6)
  res <- .isotope_peaks(x, isotopeDefinition = isoDef)
  expect_equal(res, list(c(2, 3, 4), c(7, 8, 9), c(8, 9)))
  
  res <- .isotope_peaks(x, isotopeDefinition = isoDef, duplicates = "closest")
  expect_equal(res, list(c(2, 3, 4), c(8, 9)))
  
})
