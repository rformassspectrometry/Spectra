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
