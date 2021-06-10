library(testthat)
test_that(".isotope_peaks works", {
  
  ints <- c(1, 6, 3, 5, 15)
  mz <- 1:5
  x <- cbind(mz = mz, intensity = ints)
  # plot(x, type = "h")
  substDef <- cbind(md = c(1, 2), 
                    subst_degree = c(1, 2), 
                    min_slope = c(5, 1/2),
                    max_slope = c(7, 1))
  res <- .isotope_peaks(x, substDef)
  expect_equal(res, list(c(1, 2), c(3, 5)))
  
  res <- .isotope_peaks(x, substDef, seedMz = x[c(3), 1])
  expect_equal(res, list(c(3, 5)))
})
                                                                                                 