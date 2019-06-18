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
