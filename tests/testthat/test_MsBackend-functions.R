test_that(".valid_ms_backend_files works", {
    expect_match(.valid_ms_backend_files(c("a", "b")), "File(s)")
    expect_match(.valid_ms_backend_files(c("a", "a")), "Duplicated")
    tmpf <- tempfile()
    write("hello", file = tmpf)
    expect_null(.valid_ms_backend_files(tmpf))
})

test_that(".valid_ms_backend_mod_count works", {
    expect_null(.valid_ms_backend_mod_count(1, "a"))
    expect_match(.valid_ms_backend_mod_count(integer(), "a"),
                 "Different number of")
})
