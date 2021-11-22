## assuming an object "be" exists with an initialises backend
## and a vector fls of files to be read.

test_that("backendInitialize,given MsBackend* works", {

    if (! exists("be") ) {
        stop("Missing a Spectra Backend")
    }
    if (! exists("fls") | length(fls) <1 | class(fls) != "character") {
        stop("no files passed or no filenames")
    }
        
    res1 <- backendInitialize(be, fls[1])
    n1 <- length(res1) ## 3

    ## errors
    expect_error(backendInitialize(be), "'files' is mandatory")
    expect_error(backendInitialize(be, 4), "expected to be a character")
    expect_error(suppressWarnings(backendInitialize(be, "a")), "a not found")
})
