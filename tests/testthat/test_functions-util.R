test_that(".as_integer works", {
    expect_true(is.integer(.as_integer(2.3)))
    expect_error(.as_integer("a"))
    myvar <- "2"
    expect_error(.as_integer(myvar), "Argument myvar should be")
})

test_that(".i_to_index works", {
    res <- .i_to_index(3:4, 10)
    expect_equal(res, 3:4)
    expect_equal(.i_to_index(-1, 10), -1)
    expect_error(.i_to_index(4, 3), "has to be between 1 and 3")

    expect_error(.i_to_index("a", 5), "object does not have names")
    expect_error(.i_to_index(c("a", "d"), 3, names = letters[1:3]), "not all names")
    expect_equal(.i_to_index(c("a", "d"), 4, names = letters[1:4]), c(1, 4))

    expect_error(.i_to_index(c(TRUE, FALSE), 3), "has to match the length")
    expect_equal(.i_to_index(rep(FALSE, 3), 3), integer())
    expect_equal(.i_to_index(c(FALSE, TRUE, TRUE), 3), 2:3)
})

## test_that(".is_class works", {
##     expect_true(.is_class(5L, "integer"))
##     expect_true(.is_class("a", "character"))
##     expect_true(.is_class(Rle("a", 5), "character"))
##     expect_true(.is_class(Rle(3.4, 2), "numeric"))
## })

test_that(".class_rle works", {
    expect_equal(.class_rle(1L), "integer")
    expect_equal(.class_rle(Rle(1L, 5)), "integer")
})

test_that(".as_rle works", {
    expect_equal(.as_rle("a"), "a")
    expect_equal(.as_rle(c("a", "a")), rep("a", 2))
    expect_equal(.as_rle(c("a", "a", "a")), Rle("a", 3))
    expect_equal(.as_rle(1:4), 1:4)
    expect_equal(.as_rle(rep(1, 10)), Rle(1, 10))
    expect_equal(.as_rle(TRUE), TRUE)
    expect_equal(.as_rle(c(TRUE, TRUE)), rep(TRUE, 2))
    expect_equal(.as_rle(c(TRUE, TRUE, TRUE)), Rle(TRUE, 3))
    expect_equal(.as_rle(list(a = 1:3, b = 2:3)), list(a = 1:3, b = 2:3))
})

test_that("utils.clean works", {
    x <- list(c(0, 0, 0, 1, 0, 1, 0, 0, 0),
              c(0, 0, 0, 1, 0, 0, 1, 0, 0, 0),
              c(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0),
              c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
              c(1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0),
              c(1, NA, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, NA, 1, 0, 0, 0))
    a <- list(c(FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE),
              c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,
                FALSE),
              c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE,
                FALSE, FALSE), c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE,
                FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
              c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE,
                TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
              c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE,
                TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
    b <- list(c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
              c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
                FALSE),
              c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE,
                FALSE, FALSE),
              c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
                TRUE, TRUE, TRUE, FALSE, FALSE),
              c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE),
              c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE))
    d <- list(c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
              c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
                FALSE),
              c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE,
                FALSE, FALSE),
              c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
                TRUE, TRUE, TRUE, FALSE, FALSE),
              c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE),
              c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE))
    for (i in seq(along=x)) {
      expect_identical(MSnbase:::utils.clean(x[[i]], all=TRUE), a[[i]],
                       label=paste0("all=TRUE, i=", i))
      expect_identical(MSnbase:::utils.clean(x[[i]], all=FALSE), b[[i]],
                       label=paste0("all=FALSE, i=", i))
      expect_identical(MSnbase:::utils.clean(x[[i]], na.rm=TRUE), d[[i]],
                       label=paste0("na.rm=TRUE, i=", i))
    }
})

test_that("utils.enableNeighbours works", {
    expect_error(MSnbase:::utils.enableNeighbours(1:10))
    expect_error(MSnbase:::utils.enableNeighbours(LETTERS[1:10]))
    expect_equal(MSnbase:::utils.enableNeighbours(c(FALSE, TRUE, FALSE)),
                 rep(TRUE, 3))
    expect_equal(MSnbase:::utils.enableNeighbours(c(FALSE, TRUE, FALSE, FALSE)),
                 c(rep(TRUE, 3), FALSE))
})

test_that(".filterSpectraHierarchy works", {
    msLevel <- c(1, 1, 2, 2, 2, 3, 3, 3, 4, 4)
    acquisitionNum <- 1:10
    precursorScanNum <- c(0, 0, 200, 2, 2, 4, 4, 5, 6, 20)

    expect_error(.filterSpectraHierarchy(1:3, 2, 1), "have to be the same")

    expect_equal(which(.filterSpectraHierarchy(
        acquisitionNum, precursorScanNum, 1)), 1)
    expect_equal(which(.filterSpectraHierarchy(
        acquisitionNum, precursorScanNum, 2)), c(2, 4:9))
    expect_equal(which(.filterSpectraHierarchy(
        acquisitionNum, precursorScanNum, 1:2)), c(1:2, 4:9))
    expect_equal(which(.filterSpectraHierarchy(
        acquisitionNum, precursorScanNum, 8)), c(2, 5, 8))
    expect_equal(which(.filterSpectraHierarchy(
        acquisitionNum, precursorScanNum, 9)), c(2, 4, 6, 9))
    expect_equal(which(.filterSpectraHierarchy(
        acquisitionNum, precursorScanNum, 10)), 10)
    expect_equal(which(.filterSpectraHierarchy(
        acquisitionNum, precursorScanNum, 11)), integer())
})

test_that(".rbind_fill works", {
    ## matrix
    a <- matrix(1:9, nrow = 3, ncol = 3)
    colnames(a) <- c("a", "b", "c")
    b <- matrix(1:12, nrow = 3, ncol = 4)
    colnames(b) <- c("b", "a", "d", "e")
    res <- Spectra:::.rbind_fill(a, b)
    expect_equal(colnames(res), c("a", "b", "c", "d", "e"))
    expect_equal(class(res), class(a))
    expect_equal(res[, "a"], c(a[, "a"], b[, "a"]))
    expect_equal(res[, "b"], c(a[, "b"], b[, "b"]))
    expect_equal(res[, "d"], c(NA, NA, NA, b[, "d"]))

    res <- Spectra:::.rbind_fill(a, b[, c("b", "a")])
    expect_equal(colnames(res), c("a", "b", "c"))
    expect_equal(res[, "a"], c(a[, "a"], b[, "a"]))

    ## DataFrame
    a <- DataFrame(a = 1:4, b = FALSE, c = letters[1:4])
    b <- DataFrame(d = 1:4, b = TRUE)
    res <- Spectra:::.rbind_fill(a, b)
    expect_equal(colnames(res), c("a", "b", "c", "d"))
    expect_equal(res$a, c(1:4, NA, NA, NA, NA))
    expect_equal(res$b, rep(c(FALSE, TRUE), each = 4))

    b$e <- Rle(1, 4)
    res <- Spectra:::.rbind_fill(a, b)
    expect_identical(res$e, Rle(c(NA_integer_, 1), c(4, 4)))

    a$e <- Rle(2, 4)
    res <- Spectra:::.rbind_fill(a, b)
    expect_identical(res$e, Rle(c(2, 1), c(4, 4)))

    res <- Spectra:::.rbind_fill(a, b, DataFrame(z = factor(1:5)))
    expect_identical(res$z, factor(c(rep(NA, nrow(a) + nrow(b)), 1:5)))
    expect_identical(res$a, c(1:4, rep(NA, nrow(b) + 5)))
    expect_identical(res$b, c(rep(FALSE, nrow(a)), rep(TRUE, nrow(b)), rep(NA, 5)))

    ## DataFrame containing SimpleList.
    a$mz <- NumericList(1:3, 1:4, 1:2, 1:3)
    res <- Spectra:::.rbind_fill(a, b)
    expect_true(is(res$mz, "NumericList"))
    expect_true(all(unlist(is.na(res$mz[5:8]))))
    expect_true(is.na(res$mz[[5]]))

    ## If the first DataFrame doesn't contain a SimpleList but the second one
    res <- DataFrame(d = c(1L:4L, rep(NA_real_, 4L)),
                     b = rep(c(TRUE, FALSE), each = 4L),
                     e = Rle(1:2, 4),
                     a = c(rep(NA_real_, 4L), 1L:4L),
                     c = c(rep(NA_character_, 4L), letters[1L:4L]),
                     mz = NumericList(NA_real_, NA_real_, NA_real_, NA_real_,
                                      1:3, 1:4, 1:2, 1:3))
    expect_equal(.rbind_fill(b, a), res)

    ## Ensure data types are correct after merging.
    a <- DataFrame(int = c(1L, 1L), int_rle = Rle(c(1L, 1L)),
                   log = c(TRUE, FALSE), log_rle = Rle(c(FALSE, FALSE)))
    b <- DataFrame(real = c(1.2, 1.5), real_rle = Rle(c(1.2, 1.5)))
    res <- Spectra:::.rbind_fill(a, b)
    expect_true(is.integer(res$int))
    expect_true(is.integer(res$int_rle@values))
    expect_true(is.logical(res$log))
    expect_true(is.logical(res$log_rle@values))
    expect_true(is.double(res$real))
    expect_true(is.double(res$real_rle@values))
})

test_that(".fix_breaks works", {
    rtr <- c(1, 12)
    brks <- seq(rtr[1], rtr[2], by = 4)
    ## brks do not include rtr
    expect_true(all(rtr[2] > brks))
    brks_f <- .fix_breaks(brks, rtr)
    expect_true(rtr[2] <= max(brks_f))

    ## Next unit tests taken from breaks_Spectrum.
    ## Issue #191
    ints <- 1:4
    brks <- seq(min(ints), max(ints), by = 1)
    expect_equal(.fix_breaks(brks, range(ints)), 1:5)
    expect_true(any(.fix_breaks(1:2, range(ints)) < 4))
    expect_equal(.fix_breaks(seq(1, 4, by = 2), range(ints)),
                 c(1, 3, 5))

    ## Test with values smaller than 1
    rng <- c(0.1, 0.4)
    brks <- seq(rng[1], rng[2], by = 0.04)
    .fix_breaks(brks, rng)
})

test_that(".bin_values works", {
    set.seed(123)
    vals <- abs(rnorm(20, mean = 40))
    xs <- seq(1:length(vals)) + rnorm(length(vals), sd = 0.001)

    res <- .bin_values(vals, xs, binSize = 1)
    brks <- seq(0, 20, by = 1)
    for (i in 1:(length(brks) - 1)) {
        idx <- which(xs > brks[i] & xs < brks[i +1])
        if (length(idx))
            expect_equal(res$x[i], max(vals[idx]))
        else expect_equal(res$x[i], 0)
    }

    ## Ensure that all values are within.
    xs <- seq(1:length(vals))
    brks <- seq(1, 20, by = 3)
    ## brks does not contain all values.
    expect_true(max(brks) < max(xs))
    res <- .bin_values(vals, xs, binSize = 3, fun = sum)
    ## The largest bin should contain all values larger than max(brks)
    expect_equal(res$x[length(res$x)], sum(vals[xs >= max(brks)]))

    ## Check exceptions
    expect_error(.bin_values(1:3, 1:5))
    expect_error(.bin_values(1:3, 1:5), fun = other)
})
