test_that(".valid_ms_backend_data_storage works", {
    expect_match(.valid_ms_backend_data_storage(c("a", NA)), "not allowed")
    expect_null(.valid_ms_backend_data_storage(character()))
    expect_null(.valid_ms_backend_data_storage("b"))
})

test_that(".valid_ms_backend_files_exist", {
    expect_match(.valid_ms_backend_files_exist(c("a", "b")), "a, b not")
    tmpf <- tempfile()
    write("hello", file = tmpf)
    expect_null(.valid_ms_backend_files_exist(tmpf))
    expect_null(.valid_ms_backend_files_exist(character()))
    expect_null(.valid_ms_backend_files_exist(NA_character_))
})

test_that("fillCoreSpectraVariables works", {
    a <- data.frame()
    res <- fillCoreSpectraVariables(a)
    expect_true(is.data.frame(res))
    expect_true(ncol(res) == length(coreSpectraVariables()))
    expect_equal(sort(colnames(res)), sort(names(coreSpectraVariables())))

    a <- data.frame(other_col = 1:3, msLevel = 1L)
    res <- fillCoreSpectraVariables(a)
    expect_true(all(c(names(coreSpectraVariables()), "other_col") %in%
                    colnames(res)))
    expect_equal(a$other_col, res$other_col)
    expect_equal(a$msLevel, res$msLevel)

    res <- fillCoreSpectraVariables(a, c("rtime", "precursorMz"))
    expect_equal(colnames(res),
                 c("other_col", "msLevel", "rtime", "precursorMz"))
    expect_true(all(is.na(res$rtime)))
})
