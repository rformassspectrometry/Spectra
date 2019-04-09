test_that(".valid_ms_backend_files_exist works", {
    expect_match(.valid_ms_backend_files_exist("some"), "some not found")
    expect_null(.valid_ms_backend_files_exist(character()))
    tmpf <- tempfile()
    write("hello", file = tmpf)
    expect_null(.valid_ms_backend_files_exist(tmpf))
})

test_that(".mzR_header works", {
    sc <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
    expect_error(.mzR_header(sc), "should have length 1")
    expect_error(.mzR_header(), "should have length 1")
    hdr <- .mzR_header(sc[1])
    expect_true(inherits(hdr, "DataFrame"))
    expect_equal(nrow(hdr), 931)
    expect_equal(hdr$scanIndex, 1:931)

    fl <- dir(system.file("cdf", package = "msdata"), full.names = TRUE)
    hdr <- .mzR_header(fl)
    expect_true(inherits(hdr, "DataFrame"))
    expect_equal(nrow(hdr), 1278)
    expect_equal(hdr$scanIndex, 1:1278)
})
