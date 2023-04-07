#' split:
test_that("split", {
    res <- split(be, f = seq_along(be))
    expect_true(is.list(res))
    expect_equal(length(res), length(be))
    for (i in seq_along(be)) {
        expect_s4_class(res[[i]], class(be)[1L])
        expect_true(validObject(res[[i]]))
        expect_true(length(res[[i]]) == 1L)
    }
})

#'[ , any order, duplication.
test_that("[", {
    set.seed(123)
    ## random order
    idx <- sample(seq_along(be))
    res <- be[idx]
    expect_equal(length(res), length(be))
    expect_true(validObject(res))
    for (i in seq_along(idx)) {
        a <- spectraData(res[i])
        rownames(a) <- NULL
        b <- spectraData(be[idx[i]])
        rownames(b) <- NULL
        expect_equal(a, b)
    }

    ## duplication
    res <- be[c(1, 1, 1)]
    expect_equal(length(res), 3L)
    expect_true(validObject(res))
    a <- spectraData(be[1L])
    rownames(a) <- NULL
    b <- spectraData(res[1L])
    rownames(b) <- NULL
    expect_equal(a, b)
    b <- spectraData(res[2L])
    rownames(b) <- NULL
    expect_equal(a, b)
    b <- spectraData(res[3L])
    rownames(b) <- NULL
    expect_equal(a, b)

    ## Out of range should throw error
    expect_error(be[0])

    ## integer(0) should return an empty object
    res <- be[integer()]
    expect_s4_class(res, class(be)[1L])
    expect_true(length(res) == 0L)
})

#' dropNASpectraVariables: only for not read-only
#' core spectra variables don't get removed, even if only NA.
test_that("dropNaSpectraVariables", {
    cv <- coreSpectraVariables()
    if (!isReadOnly(be) || inherits(be, "MsBackendCached") ||
        inherits(be, "MsBackendDataFrame")) {
        ## Add an empty additional variable
        tmp <- be
        tmp$na_var <- rep(NA, length(tmp))
        expect_true(any(spectraVariables(tmp) == "na_var"))
        tmp <- dropNaSpectraVariables(tmp)
        expect_false(any(spectraVariables(tmp) == "na_var"))
        expect_true(all(names(cv) %in% spectraVariables(tmp)))
    }
})

#' selectSpectraVariables: only for not read-only - and MsBackendMzR?
#' core spectra variables get the values removed, additional variables get
#' completely removed.
test_that("selectSpectraVariables", {
    if (!isReadOnly(be) || inherits(be, "MsBackendCached") ||
        inherits(be, "MsBackendDataFrame")) {
        tmp <- be
        res <- selectSpectraVariables(tmp, c("mz", "intensity",
                                             "dataStorage", "scanIndex"))
        expect_true(all(names(coreSpectraVariables()) %in%
                        spectraVariables(res)))
        expect_true(all(is.na(res$msLevel)))
        expect_true(all(is.na(res$rtime)))
        expect_true(all(is.na(res$dataOrigin)))
        expect_true(all(is.na(res$precursorMz)))
    }
})

#' test if any eventually implemented method yields the same result as the
#' default implementation
test_that("filterDataOrigin", {
    ref <- getMethod("filterDataOrigin", "MsBackend")
    org <- unique(dataOrigin(be))[1L]
    if (!is.na(org)) {
        a <- ref(be, org)
        b <- filterDataOrigin(be, org)
        a <- spectraData(a)
        rownames(a) <- NULL
        b <- spectraData(b)
        rownames(b) <- NULL
        expect_equal(a, b)
    }
})
