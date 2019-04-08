test_that("backendInitialize,MsBackend works", {
    be <- MsBackend()
    be <- backendInitialize(be)
    expect_true(validObject(be))
    expect_true(length(be@files) == 0)
    be <- backendInitialize(be, files = "a",
                            spectraData = DataFrame(fromFile = 1L))
    expect_true(validObject(be))
    expect_error(backendInitialize(be, spectraData = DataFrame(msLevel = 1L)),
                 "column(s): fromFile")
    expect_error(backendInitialize(MsBackend(),
                                   spectraData = DataFrame(msLevel = 1L,
                                                           fromFile = 1L)))
})

test_that("length,MsBackend works", {
    be <- MsBackend()
    expect_equal(length(be), 0)
    be <- new("MsBackend", spectraData = DataFrame(a = 1:3, fromFile = 1L),
              files = NA_character_, modCount = 0L)
    expect_equal(length(be), 3)
})

test_that("msLevel,MsBackend works", {
    be <- backendInitialize(MsBackend(), "a",
                            DataFrame(msLevel = c(1L, 2L, 1L),
                                      fromFile = 1L))
    expect_equal(msLevel(be), c(1, 2, 1))
    be <- backendInitialize(MsBackend(), "a",
                            DataFrame(scanIndex = 1:4,
                                      fromFile = 1L))
    expect_equal(msLevel(be), rep(NA_integer_, 4))
})
