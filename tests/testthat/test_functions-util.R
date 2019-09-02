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

test_that(".class_rle works", {
    expect_equal(.class_rle(1L), "integer")
    expect_equal(.class_rle(Rle(1L, 5)), "integer")
})

test_that(".is_class works", {
    expect_true(.is_class(3L, "integer"))
    expect_true(.is_class(Rle(1:4), "integer"))
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
      expect_identical(utils.clean(x[[i]], all=TRUE), a[[i]],
                       label=paste0("all=TRUE, i=", i))
      expect_identical(utils.clean(x[[i]], all=FALSE), b[[i]],
                       label=paste0("all=FALSE, i=", i))
      expect_identical(utils.clean(x[[i]], na.rm=TRUE), d[[i]],
                       label=paste0("na.rm=TRUE, i=", i))
    }
})

test_that("utils.enableNeighbours works", {
    expect_error(utils.enableNeighbours(1:10))
    expect_error(utils.enableNeighbours(LETTERS[1:10]))
    expect_equal(utils.enableNeighbours(c(FALSE, TRUE, FALSE)),
                 rep(TRUE, 3))
    expect_equal(utils.enableNeighbours(c(FALSE, TRUE, FALSE, FALSE)),
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
