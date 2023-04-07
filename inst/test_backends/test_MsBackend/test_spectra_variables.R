#' This test suite runs tests for accessing (and eventually setting) spectra
#' variables on a provided `MsBackend` instance and validates the returned
#' data/results. Method descriptions and expected outputs are described in
#' `MsBackend.R`. This unit tests might be run for a newly developed backend
#' extending `MsBackend` to ensure that it is compliant with it.
#'
#' To run this unit tests from another package:
#'
#' be <- <code to create and initialize the backend>
#' test_suite <- system.file("test_backends", "test_MsBackend",
#'     package = "Spectra")
#' test_dir(test_suite, stop_on_failure = TRUE)
#'
#' The unit tests in this suite expect a variable `be` to be defined, which
#' has to represent an **already initialized** backend instance.

test_that("backend is valid", {
    expect_true(validObject(be))
})

test_that("acquisitionNum", {
    res <- acquisitionNum(be)
    expect_type(res, "integer")
    expect_identical(length(res), length(be))
})

test_that("centroided", {
    res <- centroided(be)
    expect_type(res, "logical")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        centroided(tmp) <- rep(FALSE, length(be))
        expect_false(any(centroided(tmp)))
        centroided(tmp) <- rep(TRUE, length(be))
        expect_true(all(centroided(tmp)))
    }
})

test_that("collisionEnergy", {
    res <- collisionEnergy(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- abs(rnorm(length(be)))
        collisionEnergy(tmp) <- vals
        expect_equal(unname(collisionEnergy(tmp)), unname(vals))
    }
})

test_that("dataOrigin", {
    res <- dataOrigin(be)
    expect_type(res, "character")
    expect_identical(length(res), length(be))
})

test_that("dataStorage", {
    res <- dataStorage(be)
    expect_type(res, "character")
    expect_identical(length(res), length(be))
})

test_that("intensity", {
    res <- intensity(be)
    expect_true(is.list(res) || is(res, "NumericList"))
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- lapply(lengths(res), function(z) abs(rnorm(z)))
        intensity(tmp) <- vals
        res <- intensity(tmp)
        expect_equal(as.list(res), vals)
    }
})

test_that("isCentroided", {
    res <- isCentroided(be)
    expect_type(res, "logical")
    expect_identical(length(res), length(be))
})

test_that("isEmpty", {
    res <- isEmpty(be)
    expect_type(res, "logical")
    expect_identical(length(res), length(be))
    expect_identical(lengths(res) == 0, res)
})

test_that("isolationWindowLowerMz", {
    res <- isolationWindowLowerMz(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- abs(rnorm(length(be)))
        isolationWindowLowerMz(tmp) <- vals
        expect_equal(unname(isolationWindowLowerMz(tmp)), unname(vals))
    }
})

test_that("isolationWindowTargetMz", {
    res <- isolationWindowTargetMz(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- abs(rnorm(length(be)))
        isolationWindowTargetMz(tmp) <- vals
        expect_equal(unname(isolationWindowTargetMz(tmp)), unname(vals))
    }
})

test_that("isolationWindowUpperMz", {
    res <- isolationWindowUpperMz(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- abs(rnorm(length(be)))
        isolationWindowUpperMz(tmp) <- vals
        expect_equal(unname(isolationWindowUpperMz(tmp)), unname(vals))
    }
})

test_that("msLevel", {
    res <- msLevel(be)
    expect_type(res, "integer")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- sample(1:5, length(res), replace = TRUE)
        msLevel(tmp) <- vals
        expect_equal(msLevel(tmp), vals)
    }
})

test_that("mz", {
    res <- mz(be)
    expect_true(is.list(res) || is(res, "NumericList"))
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- lapply(lengths(res), function(z) abs(rnorm(z)))
        ## m/z needs to be sorted
        expect_error(mz(be) <- vals, "sorted")
        vals <- lapply(lengths(res), function(z) sort(abs(rnorm(z))))
        mz(tmp) <- vals
        res <- mz(tmp)
        expect_equal(as.list(res), vals)
    }
})

test_that("lengths", {
    res <- lengths(be)
    expect_type(res, "integer")
    expect_identical(length(res), length(be))
    expect_identical(res, lengths(mz(be)))
})

test_that("polarity", {
    ## values 0 (negative), 1 (positive), -1 (missing) NA (missing)
    res <- polarity(be)
    expect_type(res, "integer")
    expect_identical(length(res), length(be))
    expect_true(all(res %in% c(0L, 1L, -1L, NA_integer_)))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- sample(c(0L, 1L), length(res), replace = TRUE)
        polarity(tmp) <- vals
        expect_equal(unname(polarity(tmp)), unname(vals))
    }
})

test_that("precScanNum", {
    res <- precScanNum(be)
    expect_type(res, "integer")
    expect_identical(length(res), length(be))
})

test_that("precursorCharge", {
    res <- precursorCharge(be)
    expect_type(res, "integer")
    expect_identical(length(res), length(be))
})

test_that("precursorIntensity", {
    res <- precursorIntensity(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
})

test_that("precursorMz", {
    res <- precursorMz(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
})

test_that("peaksData", {
    res <- peaksData(be)
    expect_true(is.list(res))
    expect_identical(length(res), length(be))
    m <- res[[length(res)]]
    expect_true(is.matrix(m))
    expect_equal(colnames(m), c("mz", "intensity"))
    ## m/z values need to be ordered!
    expect_false(is.unsorted(m[, "mz"]))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- lapply(lengths(res) / 2, function(z)
            cbind(mz = sort(abs(rnorm(z))),
                  intensity = abs(rnorm(z))))
        peaksData(tmp) <- vals
        expect_equal(peaksData(tmp), vals)
    }
})

test_that("rtime", {
    res <- rtime(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- rnorm(length(be))
        rtime(tmp) <- vals
        expect_equal(unname(rtime(tmp)), unname(vals))
    }
})

test_that("scanIndex", {
    res <- scanIndex(be)
    expect_type(res, "integer")
    expect_identical(length(res), length(be))
})

test_that("smoothed", {
    res <- smoothed(be)
    expect_type(res, "logical")
    expect_identical(length(res), length(be))
    ## set data, if not readonly
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- sample(c(FALSE, TRUE), length(be), replace = TRUE)
        smoothed(tmp) <- vals
        expect_equal(unname(vals), unname(smoothed(tmp)))
    }
})

test_that("spectraData", {
    res <- spectraData(be)
    ## mz and intensity are required columns
    expect_s4_class(res, "DataFrame")
    expect_true(all(c("mz", "intensity") %in% colnames(res)))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- res[sample(seq_len(nrow(res))), ]
        spectraData(tmp) <- vals
        expect_equal(vals, spectraData(tmp))
    }
})

test_that("spectraNames", {
    res <- spectraNames(be)
    if (!is.null(res)) {
        expect_type(res, "character")
        expect_identical(length(res), length(be))
    }
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- as.character(sample(seq_len(length(be))))
        spectraNames(tmp) <- vals
        expect_equal(spectraNames(tmp), vals)
    }
})

test_that("spectraVariables", {
    res <- spectraVariables(be)
    expect_type(res, "character")
    expect_true(all(c("mz", "intensity") %in% res))
})

test_that("tic", {
    res <- tic(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
})

test_that("$ works", {
    res <- be$msLevel
    expect_true(is.integer(res))
    expect_equal(length(res), length(be))

    ## Error if spectra variable not available
    expect_error(be$doesnt_exist)
})
