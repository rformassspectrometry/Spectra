test_that("Spectra,DataFrame works", {
    df <- DataFrame()
    res <- Spectra(df)
    expect_true(validObject(res))
    expect_true(length(res) == 0)

    df <- DataFrame(msLevel = c(1L, 2L))
    res <- Spectra(df)
    expect_identical(msLevel(res), c(1L, 2L))
    expect_true(length(res) == 2)

    df$polarity <- "NEG"
    expect_error(Spectra(df), "wrong data type: polarity")
})

test_that("Spectra,missing works", {
    res <- Spectra()
    expect_true(length(res) == 0)

    be <- backendInitialize(MsBackendDataFrame(), files = NA_character_,
                            spectraData = DataFrame(msLevel = c(1L, 2L),
                                                    fromFile = 1L))
    res <- Spectra(backend = be)
    expect_true(length(res) == 2)
    expect_identical(msLevel(res), c(1L, 2L))
})

test_that("acquisitionNum,Spectra works", {
    sps <- Spectra()
    res <- acquisitionNum(sps)
    expect_identical(res, integer())
    df <- DataFrame(msLevel = c(1L, 1L))
    sps <- Spectra(df)
    res <- acquisitionNum(sps)
    expect_identical(res, c(NA_integer_, NA_integer_))
    df$acquisitionNum <- 1:2
    sps <- Spectra(df)
    res <- acquisitionNum(sps)
    expect_identical(res, 1:2)
})

test_that("centroided,centroided<-,Spectra works", {
    sps <- Spectra()
    res <- centroided(sps)
    expect_identical(res, logical())
    df <- DataFrame(msLevel = c(1L, 2L, 1L), centroided = c(TRUE, FALSE, TRUE))
    sps <- Spectra(df)
    res <- centroided(sps)
    expect_identical(res, c(TRUE, FALSE, TRUE))
    centroided(sps) <- c(FALSE, TRUE, FALSE)
    res <- centroided(sps)
    expect_identical(res, c(FALSE, TRUE, FALSE))
    centroided(sps) <- FALSE
    expect_true(all(centroided(sps) == FALSE))
    expect_error(centroided(sps) <- 3, "'logical' of length")
    expect_error(centroided(sps) <- c(TRUE, FALSE, FALSE, TRUE), "length 1 or 3")
})

test_that("collisionEnergy,collisionEnergy<-,Spectra works", {
    sps <- Spectra()
    res <- collisionEnergy(sps)
    expect_identical(res, numeric())

    df <- DataFrame(msLevel = c(1L, 2L, 2L))
    sps <- Spectra(df)
    res <- collisionEnergy(sps)
    expect_true(all(is.na(res)))
    expect_equal(length(res), length(sps))

    collisionEnergy(sps) <- c(1.2, 1.4, 1.7)
    res <- collisionEnergy(sps)
    expect_identical(res, c(1.2, 1.4, 1.7))

    df$collisionEnergy <- c(3.4, 4.3, 2.3)
    sps <- Spectra(df)
    res <- collisionEnergy(sps)
    expect_identical(res, c(3.4, 4.3, 2.3))

    expect_error(collisionEnergy(sps) <- 4, "of length 3")
    expect_error(collisionEnergy(sps) <- c("a", "b", "c"), "'numeric'")
})

test_that("fileNames,Spectra works", {
    sps <- Spectra()
    res <- fileNames(sps)
    expect_identical(res, character())

    df <- DataFrame(msLevel = c(1L, 2L))
    sps <- Spectra(df)
    res <- fileNames(sps)
    expect_identical(res, NA_character_)

    be <- backendInitialize(MsBackendMzR(), file = sciex_file)
    sps <- Spectra(backend = be)
    res <- fileNames(sps)
    expect_identical(res, sciex_file)
})

test_that("fromFile,Spectra works", {
    sps <- Spectra()
    expect_identical(fromFile(sps), integer())

    df <- DataFrame(msLevel = c(1L, 2L))
    sps <- Spectra(df)
    expect_identical(fromFile(sps), c(1L, 1L))
})

test_that("length,Spectra works", {
    sps <- Spectra()
    expect_equal(length(sps), 0)

    df <- DataFrame(msLevel = c(1L, 2L))
    sps <- Spectra(df)
    expect_true(length(sps) == 2)
})

test_that("intensity,Spectra works", {
    sps <- Spectra()
    res <- intensity(sps)
    expect_true(is(res, "SimpleList"))
    expect_true(length(res) == 0)

    df <- DataFrame(msLevel = c(1L, 1L),
                    rtime = c(1.2, 1.4),
                    centroided = TRUE)
    df$intensity <- list(c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                         c(9, 6, 0, 0, 3, 9, 4, 0, 0, 0))
    df$mz <- list(c(1:10), c(1:10))
    sps <- Spectra(df)
    res <- intensity(sps)
    expect_equal(res, SimpleList(
                          c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                          c(9, 6, 0, 0, 3, 9, 4, 0, 0, 0)))
    sps <- clean(sps, all = TRUE)
    res <- intensity(sps)
    expect_equal(res, SimpleList(c(1, 6, 3, 9, 1), c(9, 6, 3, 9, 4)))
})

test_that("ionCount,Spectra works", {
    sps <- Spectra()
    res <- ionCount(sps)
    expect_identical(res, numeric())

    df <- DataFrame(msLevel = c(1L, 2L), centroided = TRUE)
    df$intensity <- list(c(5, 9, 3), c(9, 8, 2))
    df$mz <- list(1:3, 1:3)
    sps <- Spectra(df)

    res <- ionCount(sps)
    expect_identical(res, c(17, 19))

    sps <- removePeaks(sps, t = 4)
    res <- ionCount(sps)
    expect_identical(res, c(14, 17))
})

test_that("isCentroided,Spectra works", {
    sps <- Spectra()
    res <- isCentroided(sps)
    expect_identical(res, logical())

    df <- DataFrame(msLevel = c(1L, 1L))
    df$intensity <- list(c(5, 6, 1), c(5, 3, 1))
    df$mz <- list(1:3, 1:3)
    sps <- Spectra(df)

    res <- isCentroided(sps)
    expect_identical(res, c(NA, NA))

    sps <- Spectra(backendInitialize(MsBackendMzR(), file = sciex_file))
    res <- isCentroided(sps)
    expect_true(length(res) == length(sps))
    expect_true(all(!res))
})

test_that("isEmpty,Spectra works", {
    sps <- Spectra()
    res <- isEmpty(sps)
    expect_identical(res, logical())

    df <- DataFrame(msLevel = c(2L, 2L), centroided = TRUE)
    df$intensity <- list(c(4, 6, 1), c(45, 2))
    df$mz <- list(1:3, 1:2)

    sps <- Spectra(df)
    res <- isEmpty(sps)
    expect_identical(res, c(FALSE, FALSE))

    sps <- removePeaks(sps, t = 100)
    res <- isEmpty(sps)
    expect_identical(res, c(FALSE, FALSE))

    sps <- clean(sps, all = TRUE)
    res <- isEmpty(sps)
    expect_identical(res, c(TRUE, TRUE))
})

test_that("msLevel,Spectra works", {
    sps <- Spectra()
    expect_identical(msLevel(sps), integer())

    df <- DataFrame(msLevel = c(1L, 2L))
    sps <- Spectra(df)
    expect_identical(msLevel(sps), c(1L, 2L))
})

test_that("mz,Spectra works", {
    sps <- Spectra()
    res <- mz(sps)
    expect_true(is(res, "SimpleList"))
    expect_true(length(res) == 0)

    df <- DataFrame(msLevel = c(1L, 1L),
                    rtime = c(1.2, 1.4),
                    centroided = TRUE)
    df$intensity <- list(c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                         c(9, 6, 0, 0, 3, 9, 4, 0, 0, 0))
    df$mz <- list(c(1:10), c(1:10))
    sps <- Spectra(df)
    res <- mz(sps)
    expect_equal(res, SimpleList(1:10, 1:10))
    sps <- clean(sps, all = TRUE)
    res <- mz(sps)
    expect_equal(res, SimpleList(c(3, 4, 5, 8, 9), c(1, 2, 5, 6, 7)))
})

test_that("peaks,Spectra works", {
    df <- DataFrame(msLevel = c(1L, 2L), fromFile = 1L)
    df$mz <- list(1:4, 1:5)
    df$intensity <- list(1:4, 1:5)
    be <- backendInitialize(MsBackendDataFrame(), file = NA_character_, df)
    sps <- Spectra(backend = be)
    res <- peaks(sps)
    expect_true(is(res, "SimpleList"))
    expect_equal(res[[1]][, 1], 1:4)
    expect_equal(res[[2]][, 1], 1:5)

    sps <- Spectra(backend = MsBackendDataFrame())
    res <- peaks(sps)
    expect_true(is(res, "SimpleList"))
    expect_true(length(res) == 0)
})

test_that("peaksCount,Spectra works", {
    sps <- Spectra()
    res <- peaksCount(sps)
    expect_identical(res, integer())

    df <- DataFrame(msLevel = c(2L, 2L), centroided = TRUE)
    df$intensity <- list(c(4, 6, 1), c(45, 2))
    df$mz <- list(1:3, 1:2)

    sps <- Spectra(df)
    res <- peaksCount(sps)
    expect_identical(res, c(3L, 2L))

    sps <- removePeaks(sps, t = 100)
    res <- peaksCount(sps)
    expect_identical(res, c(3L, 2L))

    sps <- clean(sps, all = TRUE)
    res <- peaksCount(sps)
    expect_identical(res, c(0L, 0L))
})

test_that("polarity,polarity<-,Spectra works", {
    sps <- Spectra()
    expect_identical(polarity(sps), integer())

    df <- DataFrame(msLevel = c(1L, 1L))
    sps <- Spectra(df)
    expect_identical(polarity(sps), c(NA_integer_, NA_integer_))
    polarity(sps) <- c(1L, -1L)
    expect_identical(polarity(sps), c(1L, -1L))
    polarity(sps) <- 2L
    expect_identical(polarity(sps), c(2L, 2L))

    expect_error(polarity(sps) <- c(1L, 2L, 1L), "of length 1 or 2")
    expect_error(polarity(sps) <- c("a", "b"), "an 'integer'")
})

test_that("precScanNum,Spectra works", {
    sps <- Spectra()
    expect_identical(precScanNum(sps), integer())

    sps <- Spectra(DataFrame(msLevel = c(2L, 2L)))
    expect_identical(precScanNum(sps), c(NA_integer_, NA_integer_))
    sps <- Spectra(DataFrame(msLevel = c(2L, 2L), precScanNum = c(1L, 2L)))
    expect_identical(precScanNum(sps), c(1L, 2L))
})

test_that("precursorCharge,Spectra works", {
    sps <- Spectra()
    expect_identical(precursorCharge(sps), integer())

    sps <- Spectra(DataFrame(msLevel = c(2L, 2L)))
    expect_identical(precursorCharge(sps), c(NA_integer_, NA_integer_))
    sps <- Spectra(DataFrame(msLevel = c(2L, 2L), precursorCharge = c(1L, 3L)))
    expect_identical(precursorCharge(sps), c(1L, 3L))
})

test_that("precursorIntensity,Spectra works", {
    sps <- Spectra()
    expect_identical(precursorIntensity(sps), numeric())

    sps <- Spectra(DataFrame(msLevel = c(2L, 2L)))
    expect_identical(precursorIntensity(sps), c(NA_real_, NA_real_))
    sps <- Spectra(DataFrame(msLevel = c(2L, 2L),
                             precursorIntensity = c(1.3, 3.2)))
    expect_identical(precursorIntensity(sps), c(1.3, 3.2))
})

test_that("precursorMz,Spectra works", {
    sps <- Spectra()
    expect_identical(precursorMz(sps), numeric())

    sps <- Spectra(DataFrame(msLevel = c(2L, 2L)))
    expect_identical(precursorMz(sps), c(NA_real_, NA_real_))
    sps <- Spectra(DataFrame(msLevel = c(2L, 2L),
                             precursorMz = c(234.2, 668.2)))
    expect_identical(precursorMz(sps), c(234.2, 668.2))
})

test_that("rtime,rtime<-,Spectra works", {
    sps <- Spectra()
    expect_identical(rtime(sps), numeric())

    sps <- Spectra(DataFrame(msLevel = c(1L, 2L)))
    expect_identical(rtime(sps), c(NA_real_, NA_real_))
    sps <- Spectra(DataFrame(msLevel = c(1L, 2L), rtime = c(2.1, 4.2)))
    expect_identical(rtime(sps), c(2.1, 4.2))
    rtime(sps) <- c(1.2, 1.3)
    expect_identical(rtime(sps), c(1.2, 1.3))

    expect_error(rtime(sps) <- c(TRUE, FALSE), "'numeric' of length 2")
    expect_error(rtime(sps) <- 1.3, "'numeric' of length 2")
})

test_that("scanIndex,Spectra works", {
    sps <- Spectra()
    expect_identical(scanIndex(sps), integer())
    sps <- Spectra(DataFrame(msLevel = c(1L, 2L)))
    expect_identical(scanIndex(sps), c(NA_integer_, NA_integer_))

    sps <- Spectra(DataFrame(msLevel = c(1L, 2L), scanIndex = c(1L, 2L)))
    expect_identical(scanIndex(sps), c(1L, 2L))
})

test_that("selectSpectraVariables,Spectra works", {
    sps <- Spectra()
    res <- selectSpectraVariables(sps, c("msLevel", "rtime"))
    expect_equal(sps, res)

    sps <- Spectra(DataFrame(msLevel = c(1L, 2L), rtime = c(1.2, 1.4),
                             other_col = "a"))
    res <- selectSpectraVariables(sps, c("rtime", "other_col"))
    expect_equal(msLevel(res), c(NA_integer_, NA_integer_))
    expect_equal(rtime(res), c(1.2, 1.4))
    expect_equal(spectraData(res, "other_col")[, 1], c("a", "a"))

    res <- selectSpectraVariables(sps, c("msLevel"))
    expect_identical(rtime(res), c(NA_real_, NA_real_))
    expect_identical(msLevel(res), 1:2)
    expect_true(!any(spectraVariables(res) %in% "other_col"))

    expect_error(selectSpectraVariables(sps, c("rtime", "something")), "not")
})

test_that("smoothed,smoothed<-,Spectra works", {
    sps <- Spectra()
    expect_identical(smoothed(sps), logical())

    sps <- Spectra(DataFrame(msLevel = 1:2))
    expect_identical(smoothed(sps), c(NA, NA))
    sps <- Spectra(DataFrame(msLevel = 1:2, smoothed = TRUE))
    expect_identical(smoothed(sps), c(TRUE, TRUE))

    smoothed(sps) <- c(FALSE, TRUE)
    expect_identical(smoothed(sps), c(FALSE, TRUE))
    smoothed(sps) <- FALSE
    expect_identical(smoothed(sps), c(FALSE, FALSE))

    expect_error(smoothed(sps) <- c("a", "b"), "of length 1 or 2")
    expect_error(smoothed(sps) <- c(TRUE, FALSE, TRUE), "of length 1 or 2")
})

test_that("spectraData,Spectra works", {
    sps <- Spectra()
    res <- spectraData(sps)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0)
    expect_true(all(c("msLevel", "rtime", "mz", "intensity") %in% colnames(res)))

    df <- DataFrame(msLevel = c(1L, 1L), rtime = c(1.2, 1.3), centroided = TRUE)
    df$mz <- list(1:10, 1:10)
    df$intensity <- list(c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                         c(9, 6, 0, 0, 3, 0, 0, 0, 3, 2))
    sps <- Spectra(df)
    res <- spectraData(sps)
    expect_true(is(res, "DataFrame"))
    expect_equal(res$msLevel, c(1L, 1L))
    expect_equal(res$rtime, c(1.2, 1.3))
    expect_equal(res$precScanNum, c(NA_integer_, NA_integer_))
    expect_equal(res$mz, SimpleList(1:10, 1:10))
    expect_equal(res$intensity, SimpleList(c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                                           c(9, 6, 0, 0, 3, 0, 0, 0, 3, 2)))
    sps <- clean(sps, all = TRUE)
    res <- spectraData(sps)
    expect_equal(res$mz, SimpleList(c(3, 4, 5, 8, 9), c(1, 2, 5, 9, 10)))
    expect_equal(res$intensity, SimpleList(c(1, 6, 3, 9, 1), c(9, 6, 3, 3, 2)))
})

test_that("spectraData<-,Spectra works", {
    df <- DataFrame(msLevel = c(1L, 1L), rtime = c(1.2, 1.3), centroided = TRUE)
    df$mz <- list(1:10, 1:10)
    df$intensity <- list(c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                         c(9, 6, 0, 0, 3, 0, 0, 0, 3, 2))
    sps <- Spectra(df)

    spectraData(sps)$add_col <- 3
    expect_true(any(spectraVariables(sps) == "add_col"))
    expect_equal(spectraData(sps)$add_col, c(3, 3))

    df$mz <- list(11:20, 11:20)
    df$precursorMz <- c(0, 0)
    df$fromFile <- 1L

    spectraData(sps) <- df
    expect_true(!any(spectraVariables(sps) == "add_col"))
    expect_equal(mz(sps), SimpleList(11:20, 11:20))

    expect_error(spectraData(sps) <- c(1, 2, 4, 5), "has to be 2")

    expect_error(spectraData(sps) <- df[c(1, 1, 2), ], "with 2 rows")

    tmp <- sciex_mzr
    sps <- Spectra(tmp)
    spectraData(sps)$some_col <- "yes"
    expect_true(any(spectraVariables(sps) == "some_col"))
    expect_true(all(spectraData(sps, "some_col")[, 1] == "yes"))
    expect_true(is(sps@backend@spectraData$some_col, "Rle"))

    spd <- spectraData(sps)
    spd$msLevel <- 2L
    spectraData(sps) <- spd
    expect_true(all(msLevel(sps) == 2L))
})

test_that("spectraNames,spectraNames<-,Spectra works", {
    sps <- Spectra()
    expect_identical(spectraNames(sps), NULL)

    sps <- Spectra(DataFrame(msLevel = c(1L, 1L)))
    expect_identical(spectraNames(sps), NULL)

    spectraNames(sps) <- c("a", "b")
    expect_identical(spectraNames(sps), c("a", "b"))

    expect_error(spectraNames(sps) <- 1, "invalid rownames length")
})

test_that("spectraVariables,Spectra works", {
    sps <- Spectra()
    res <- spectraVariables(sps)
    exp_col <- c("msLevel", "rtime", "acquisitionNum", "scanIndex", "mz",
                 "intensity", "fromFile", "centroided", "smoothed",
                 "polarity", "precScanNum", "precursorMz", "precursorIntensity",
                 "precursorCharge", "collisionEnergy")
    expect_true(all(exp_col %in% res))
    df <- DataFrame(other_col = "a", msLevel = c(1L, 1L))
    sps <- Spectra(df)
    res <- spectraVariables(sps)
    expect_true(all(c(exp_col, "other_col") %in% res))
})

test_that("tic,Spectra works", {
    sps <- Spectra()
    expect_identical(tic(sps), numeric())

    sps <- Spectra(DataFrame(msLevel = c(1L, 1L)))
    expect_identical(tic(sps), c(NA_real_, NA_real_))

    df <- DataFrame(msLevel = c(1L, 1L), totIonCurrent = c(5, 3))
    sps <- Spectra(df)
    expect_identical(tic(sps), c(5, 3))
    expect_identical(tic(sps, initial = FALSE), c(0, 0))

    df$intensity <- list(c(3, 3, 1), c(5, 3, 1))
    df$mz <- list(1:3, 1:3)
    sps <- Spectra(df)
    expect_identical(tic(sps, initial = FALSE), c(7, 9))
})

test_that("$, $<-,Spectra works", {
    sps <- Spectra()
    expect_identical(sps$msLevel, integer())

    sps <- Spectra(DataFrame(msLevel = c(1L, 2L), other_column = "a"))
    expect_identical(sps$msLevel, c(1L, 2L))
    expect_identical(sps$other_column, c("a", "a"))

    sps$second_col <- c(1, 4)
    expect_identical(sps$second_col, c(1, 4))
    expect_true(any(spectraVariables(sps) == "second_col"))

    df <- DataFrame(msLevel = c(1L, 2L), other_column = "a")
    df$mz <- list(1:3, 1:4)
    df$intensity <- list(c(3, 6, 3), c(56, 6, 3, 2))
    sps <- Spectra(df)

    sps$intensity <- list(c(4, 4, 4), c(2, 4, 6, 3))
    expect_equal(intensity(sps), SimpleList(c(4, 4, 4), c(2, 4, 6, 3)))

    sps <- Spectra(sciex_mzr)
    sps$add_col <- "something"
    expect_true(all(sps$add_col == "something"))
})

#### ---------------------------------------------------------------------------
##
##                      FILTERING AND SUBSETTING
##
#### ---------------------------------------------------------------------------

test_that("[,Spectra works", {
    df <- DataFrame(msLevel = c(1L, 2L, 1L, 2L),
                    rtime = c(1, 2, 1, 2))
    sps <- Spectra(df)
    res <- sps[c(2, 4), ]
    expect_identical(msLevel(res), c(2L, 2L))
    expect_identical(sps, sps[])

    expect_error(sps[, 4], "columns is not")

    res <- sps[integer()]
    expect_true(is(res, "Spectra"))
    expect_true(length(res) == 0)

    sps <- Spectra(sciex_mzr)
    tmp <- sps[fromFile(sps) == 2L, ]
    expect_true(all(fromFile(tmp) == 1))
    expect_equal(fileNames(tmp), fileNames(sps)[2])
    expect_equal(rtime(tmp), rtime(sps)[fromFile(sps) == 2])
})

#### ---------------------------------------------------------------------------
##
##                      DATA MANIPULATION METHODS
##
#### ---------------------------------------------------------------------------

test_that("clean,Spectra works", {
    sps <- Spectra()
    res <- clean(sps)
    expect_true(length(res@processingQueue) == 1)
    expect_equal(res@processingQueue[[1]],
                 ProcessingStep(.clean_peaks,
                                list(all = FALSE, msLevel. = integer())))

    res <- clean(sps, all = TRUE, msLevel. = 2L)
    expect_true(length(res@processingQueue) == 1)
    expect_equal(res@processingQueue[[1]],
                 ProcessingStep(.clean_peaks,
                                list(all = TRUE, msLevel. = 2L)))
})

test_that("removePeaks,Spectra works", {
    sps <- Spectra()
    res <- removePeaks(sps, t = 10)
    expect_true(length(res@processingQueue) == 1)
    expect_equal(res@processingQueue[[1]],
                 ProcessingStep(.remove_peaks,
                                list(t = 10, msLevel. = integer())))
})
