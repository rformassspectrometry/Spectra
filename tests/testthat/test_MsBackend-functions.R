test_that(".valid_ms_backend_files works", {
    expect_match(.valid_ms_backend_files(c("a", "a")), "Duplicated")
    expect_null(.valid_ms_backend_files(character()))
    expect_null(.valid_ms_backend_files(NA_character_))
})

test_that(".valid_ms_backend_files_exist", {
    expect_match(.valid_ms_backend_files_exist(c("a", "b")), "a, b not")
    tmpf <- tempfile()
    write("hello", file = tmpf)
    expect_null(.valid_ms_backend_files_exist(tmpf))
    expect_null(.valid_ms_backend_files_exist(character()))
    expect_null(.valid_ms_backend_files_exist(NA_character_))
})

test_that(".valid_ms_backend_mod_count works", {
    expect_null(.valid_ms_backend_mod_count(1, "a"))
    expect_match(.valid_ms_backend_mod_count(integer(), "a"),
                 "Different number of")
})
