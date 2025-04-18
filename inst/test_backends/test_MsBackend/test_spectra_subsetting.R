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

    ## logical
    l <- rep(FALSE, length(be))
    l[sample(seq_along(l), floor(length(l) / 2))] <- TRUE
    res <- be[l]
    expect_true(validObject(res))
    expect_true(length(res) == sum(l))
    expect_equal(res, be[which(l)])
})

test_that("cbind2 works", {
    seql <- length(be)
    df <- data.frame(cola = seq_len(seql), colb = "b", colz = "z")
    res <- cbind2(be, df)
    expect_true(validObject(res))
    expect_equal(ncol(spectraData(res)), length(spectraVariables(be)) + 3)
    expect_equal(res$cola, seq_len(seql))
    expect_equal(res$colb, rep("b", seql))
    expect_equal(res$colz, rep("z", seql))
    df2  <- data.frame(cola = seq_len(length(be) / 2), colb = "b", colz = "z")
    expect_error(cbind2(be, df2), "does not match")
    ## with matrix
    m <- matrix(1:seql, ncol = 1, dimnames = list(NULL, "m"))
    res <- cbind2(be, m)
    expect_true(validObject(res))
    expect_equal(ncol(spectraData(res)), length(spectraVariables(be)) + 1)
    expect_equal(res$m, 1:seql)
    ## no replacing
    expect_error(cbind2(be, data.frame(scanIndex = 1:seql)),
                 "are already present")
})

#' extractByIndex. Uses [ if not implemented
test_that("extractByIndex", {
    i <- sample(seq_along(be), floor(length(be) / 2))
    res <- extractByIndex(be, i)
    expect_true(validObject(res))
    expect_equal(length(res), length(i))
    expect_equal(msLevel(res), msLevel(be)[i])
    expect_equal(rtime(res), rtime(be)[i])
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
        res <- selectSpectraVariables(
            tmp, union(c("mz", "intensity", "dataStorage", "scanIndex"),
                       backendRequiredSpectraVariables(be)))
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

#' If `msLevel.` is not provided or if it maches all unique MS levels the
#' filter should be applied to all spectra, otherwise to only the spectra of
#' the MS level.
test_that("filterRt works", {
    ref <- getMethod("filterRt", "MsBackend")
    rtr <- range(rtime(be))
    rtr <- range(c(rtr[1L] * 2, rtr[2L] / 2))
    res_ref <- ref(be, rtr, msLevel. = integer())
    res <- filterRt(be, rtr, msLevel. = integer())
    expect_equal(length(res_ref), length(res))
    expect_equal(rtime(res_ref), rtime(res))

    res_ref <- ref(be, rtr, msLevel. = unique(msLevel(be)))
    res <- filterRt(be, rtr, msLevel. = unique(msLevel(be)))
    expect_equal(length(res_ref), length(res))
    expect_equal(rtime(res_ref), rtime(res))

    res_ref <- ref(be, rtr, msLevel. = 1L)
    res <- filterRt(be, rtr, msLevel. = 1L)
    expect_equal(length(res_ref), length(res))
    expect_equal(rtime(res_ref), rtime(res))

    res_ref <- ref(be, rtr, msLevel. = 2L)
    res <- filterRt(be, rtr, msLevel. = 2L)
    expect_equal(length(res_ref), length(res))
    expect_equal(rtime(res_ref), rtime(res))
})
