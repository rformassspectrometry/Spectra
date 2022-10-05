test_df <- DataFrame(msLevel = c(1L, 2L, 2L), scanIndex = 4:6)
test_df$mz <- list(c(1.1, 1.3, 1.5), c(4.1, 5.1), c(1.6, 1.7, 1.8, 1.9))
test_df$intensity <- list(c(45.1, 34, 12), c(234.4, 1333), c(42.1, 34.2, 65, 6))

test_that("backendInitialize,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_true(validObject(be))

    be <- backendInitialize(be, data = DataFrame(msLevel = 2L))
    expect_true(validObject(be))
    expect_identical(dataStorage(be), "<memory>")
    expect_true(length(be@peaksData) == nrow(be@spectraData))
    expect_equal(be@peaksData[[1L]],
                 matrix(numeric(), ncol = 2, nrow = 0,
                        dimnames = list(character(), c("mz", "intensity"))))

    be_2 <- backendInitialize(be, data = data.frame(msLevel = 2L))
    expect_equal(be, be_2)

    expect_error(backendInitialize(be, data = 4), "has to be a")

    df <- test_df
    be <- backendInitialize(be, df)
    expect_true(validObject(be))
    expect_equal(be@spectraData[, c("msLevel", "scanIndex")],
                 as.data.frame(df[, c("msLevel", "scanIndex")]))
    expect_true(length(be@peaksData) == 3)
    expect_equal(be@peaksData[[1L]], cbind(mz = df$mz[[1L]],
                                           intensity = df$intensity[[1L]]))
    expect_equal(be@peaksData[[2L]], cbind(mz = df$mz[[2L]],
                                           intensity = df$intensity[[2L]]))
    expect_equal(be@peaksData[[3L]], cbind(mz = df$mz[[3L]],
                                           intensity = df$intensity[[3L]]))

    df$mz <- NULL
    be <- backendInitialize(be, df)
    expect_true(validObject(be))
    expect_equal(be@peaksData[[1L]], cbind(intensity = df$intensity[[1L]]))

    df <- test_df
    df$other_col <- list(1:3, 1:4, 1:3)
    be <- backendInitialize(be, df)
    expect_equal(be@peaksDataFrame, list())
    expect_true(any(colnames(be@spectraData) == "other_col"))

    df$other_col <- list(1:3, 1:2, 1:4)
    be <- backendInitialize(be, df)
    expect_equal(be@peaksDataFrame, list(data.frame(other_col = 1:3),
                                         data.frame(other_col = 1:2),
                                         data.frame(other_col = 1:4)))

    df$yet_another_col <- list(c("a", "b", "c"), c("a", "b"), letters[3:6])
    be <- backendInitialize(be, df)
    expect_equal(be@peaksDataFrame,
                 list(data.frame(other_col = 1:3,
                                 yet_another_col = c("a", "b", "c")),
                      data.frame(other_col = 1:2,
                                 yet_another_col = c("a", "b")),
                      data.frame(other_col = 1:4,
                                 yet_another_col = letters[3:6])))
})

test_that("show,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_output(show(be), "MsBackendMemory")
})

test_that("dataStorage,dataStorage<-MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_equal(dataStorage(be), character())

    be <- backendInitialize(be, test_df)
    expect_equal(dataStorage(be), rep("<memory>", 3))

    dataStorage(be) <- c("other", "storage",  "mode")
    expect_equal(dataStorage(be), c("other", "storage", "mode"))

    expect_error(dataStorage(be) <- 3, "length 3")
    expect_error(dataStorage(be) <- 1:3, "'character'")
})

test_that("length,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_true(length(be) == 0)
    be <- backendInitialize(be, test_df)
    expect_true(length(be) == 3)
})

test_that("spectraVariables,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    res <- spectraVariables(be)
    expect_equal(res, names(coreSpectraVariables()))

    df <- data.frame(new_col = "a")
    df$mz <- list(1:3)
    be <- backendInitialize(be, df)
    res <- spectraVariables(be)
    expect_equal(res, c(names(coreSpectraVariables()), "new_col"))
})

test_that("peaksVariables,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    res <- peaksVariables(be)
    expect_equal(res, c("mz", "intensity"))

    df <- data.frame(msLevel = 1:2)
    df$mz <- list(1:2, 1:3)
    be <- backendInitialize(be, df)
    res <- peaksVariables(be)
    expect_equal(res, "mz")

    df$peak_ann <- list(c("a", "b"), c("a", "b", "c"))
    be <- backendInitialize(be, df)
    res <- peaksVariables(be)
    expect_equal(res, c("mz", "peak_ann"))
})

test_that("lengths,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_equal(lengths(be), integer())

    be <- backendInitialize(be, test_df)
    expect_equal(lengths(be), c(3L, 2L, 4L))
})

test_that("mz,mz<-,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_equal(mz(be), IRanges::NumericList(compress = FALSE))

    be <- backendInitialize(be, test_df)
    expect_equal(mz(be), IRanges::NumericList(test_df$mz, compress = FALSE))

    ## Replacing
    expect_error(mz(be) <- "a", "list or")
    expect_error(mz(be) <- list(1:3), "length of 'object'")
    expect_error(mz(be) <- list(1:3, 1:3, 1:3), "number of peaks")

    vals <- list(c(1.4, 5.2, 1.3), c(6.2, 5.6), c(1.1, 3.3, 4.4, 5.5))
    mz(be) <- vals
    expect_equal(mz(be), IRanges::NumericList(vals, compress = FALSE))

    tmp <- test_df
    tmp$mz <- NULL
    be <- backendInitialize(be, tmp)
    mz(be) <- vals
    expect_equal(mz(be), IRanges::NumericList(vals, compress = FALSE))
    expect_equal(colnames(be@peaksData[[1L]]), c("intensity", "mz"))
})

test_that("intensity,intensity<-,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_equal(intensity(be), IRanges::NumericList(compress = FALSE))

    be <- backendInitialize(be, test_df)
    expect_equal(intensity(be), IRanges::NumericList(test_df$intensity,
                                                     compress = FALSE))

    ## Replacing
    expect_error(intensity(be) <- "a", "list or")
    expect_error(intensity(be) <- list(1:3), "length of 'object'")
    expect_error(intensity(be) <- list(1:3, 1:3, 1:3), "number of peaks")

    vals <- list(c(1.4, 5.2, 1.3), c(6.2, 5.6), c(1.1, 3.3, 4.4, 5.5))
    intensity(be) <- vals
    expect_equal(intensity(be), IRanges::NumericList(vals, compress = FALSE))

    tmp <- test_df
    tmp$intensity <- NULL
    be <- backendInitialize(be, tmp)
    intensity(be) <- vals
    expect_equal(intensity(be), IRanges::NumericList(vals, compress = FALSE))
    expect_equal(colnames(be@peaksData[[1L]]), c("mz", "intensity"))
})

test_that("spectraData,MsBackendMemory works", {
    be <- new("MsBackendMemory")

    res <- spectraData(be)
    expect_s4_class(res, "DataFrame")
    expect_true(all(colnames(res) %in% names(coreSpectraVariables())))

    be <- backendInitialize(be, test_df)
    res <- spectraData(be)
    expect_true(all(names(coreSpectraVariables()) %in% colnames(res)))

    expect_equal(res$msLevel, test_df$msLevel)
    expect_equal(res$scanIndex, test_df$scanIndex)
    expect_true(all(res$dataStorage == "<memory>"))
    expect_equal(res$mz, IRanges::NumericList(test_df$mz, compress = FALSE))
    expect_equal(res$intensity,
                 IRanges::NumericList(test_df$intensity, compress = FALSE))

    tmp <- test_df
    tmp$pk_anno <- list(c("a", "b", "c"), c("", "d"), letters[12:15])
    be <- backendInitialize(be, tmp)
    expect_true(length(be@peaksDataFrame) == 3)
    res <- spectraData(be)
    expect_equal(res$msLevel, test_df$msLevel)
    expect_equal(res$scanIndex, test_df$scanIndex)
    expect_true(all(res$dataStorage == "<memory>"))
    expect_equal(res$mz, IRanges::NumericList(test_df$mz, compress = FALSE))
    expect_equal(res$intensity,
                 IRanges::NumericList(test_df$intensity, compress = FALSE))
    expect_equal(res$pk_anno, tmp$pk_anno)

    tmp$add_anno <- list(c(1:3), 1:2, 1:4)
    be <- backendInitialize(be, tmp)
    res <- spectraData(be)
    expect_equal(res$pk_anno, tmp$pk_anno)
    expect_equal(res$add_anno, tmp$add_anno)

    res <- spectraData(be, "mz")
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), "mz")
    expect_equal(res$mz, IRanges::NumericList(tmp$mz, compress = FALSE))
    res <- spectraData(be, "msLevel")
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), "msLevel")
    expect_equal(res$msLevel, tmp$msLevel)
    res <- spectraData(be, "rtime")
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), "rtime")
    expect_equal(res$rtime, rep(NA_real_, 3))
    res <- spectraData(be, "pk_anno")
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), "pk_anno")
    expect_equal(res$pk_anno, tmp$pk_anno)
})

test_that("spectraData<-,MsBackendMemory works", {
    be <- backendInitialize(new("MsBackendMemory"), test_df)

    newDF <- test_df
    newDF$rtime <- 1:3
    newDF$mz <- list(1:3, 1:2, 1:4)
    spectraData(be) <- newDF
    expect_equal(be@spectraData$rtime, 1:3)
    expect_true(all(dataStorage(be) == "<memory>"))
    expect_equal(mz(be), IRanges::NumericList(list(1:3, 1:2, 1:4),
                                              compress = FALSE))
    expect_error(spectraData(be) <- newDF[1:2, ], "3 rows")
    expect_error(spectraData(be) <- "a", "DataFrame")
})

test_that("peaksData,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    res <- peaksData(be)
    expect_equal(res, list())

    be <- backendInitialize(be, test_df)
    res <- peaksData(be)
    expect_equal(res[[1L]], cbind(mz = test_df$mz[[1L]],
                                  intensity = test_df$intensity[[1L]]))
    expect_equal(res[[2L]], cbind(mz = test_df$mz[[2L]],
                                  intensity = test_df$intensity[[2L]]))
    expect_equal(res[[3L]], cbind(mz = test_df$mz[[3L]],
                                  intensity = test_df$intensity[[3L]]))
    expect_error(peaksData(be, "other"), "not available")
    res <- peaksData(be, c("intensity", "mz"))
    expect_equal(res[[1L]], cbind(intensity = test_df$intensity[[1L]],
                                  mz = test_df$mz[[1L]]))
    expect_equal(res[[2L]], cbind(intensity = test_df$intensity[[2L]],
                                  mz = test_df$mz[[2L]]))
    expect_equal(res[[3L]], cbind(intensity = test_df$intensity[[3L]],
                                  mz = test_df$mz[[3L]]))
    res <- peaksData(be, c("mz"))
    expect_equal(res[[1L]], cbind(mz = test_df$mz[[1L]]))
    expect_equal(res[[2L]], cbind(mz = test_df$mz[[2L]]))
    expect_equal(res[[3L]], cbind(mz = test_df$mz[[3L]]))

    tmp <- test_df
    tmp$pk_ann <- list(letters[1:3], c("", ""), letters[1:4])
    be <- backendInitialize(be, tmp)
    res <- peaksData(be)
    expect_equal(res[[1L]], cbind(mz = test_df$mz[[1L]],
                                  intensity = test_df$intensity[[1L]]))
    expect_equal(res[[2L]], cbind(mz = test_df$mz[[2L]],
                                  intensity = test_df$intensity[[2L]]))
    expect_equal(res[[3L]], cbind(mz = test_df$mz[[3L]],
                                  intensity = test_df$intensity[[3L]]))
    res <- peaksData(be, c("mz", "pk_ann"))
    expect_equal(res[[1L]], cbind(mz = tmp$mz[[1L]],
                                  pk_ann = tmp$pk_ann[[1L]]))
    expect_equal(res[[2L]], cbind(mz = tmp$mz[[2L]],
                                  pk_ann = tmp$pk_ann[[2L]]))
    expect_equal(res[[3L]], cbind(mz = tmp$mz[[3L]],
                                  pk_ann = tmp$pk_ann[[3L]]))
    res <- peaksData(be, c("pk_ann", "mz"))
    expect_equal(res[[1L]], cbind(pk_ann = tmp$pk_ann[[1L]],
                                  mz = tmp$mz[[1L]]))
    expect_equal(res[[2L]], cbind(pk_ann = tmp$pk_ann[[2L]],
                                  mz = tmp$mz[[2L]]))
    expect_equal(res[[3L]], cbind(pk_ann = tmp$pk_ann[[3L]],
                                  mz = tmp$mz[[3L]]))
    tmp$add_ann <- list(1:3, 1:2, 1:4)
    be <- backendInitialize(be, tmp)
    res <- peaksData(be, c("mz", "pk_ann"))
    expect_equal(res[[1L]], cbind(mz = tmp$mz[[1L]],
                                  pk_ann = tmp$pk_ann[[1L]]))
    expect_equal(res[[2L]], cbind(mz = tmp$mz[[2L]],
                                  pk_ann = tmp$pk_ann[[2L]]))
    expect_equal(res[[3L]], cbind(mz = tmp$mz[[3L]],
                                  pk_ann = tmp$pk_ann[[3L]]))
    res <- peaksData(be, c("add_ann", "pk_ann"))
    expect_equal(res[[1L]], cbind(add_ann = tmp$add_ann[[1L]],
                                  pk_ann = tmp$pk_ann[[1L]]))
    expect_equal(res[[2L]], cbind(add_ann = tmp$add_ann[[2L]],
                                  pk_ann = tmp$pk_ann[[2L]]))
    expect_equal(res[[3L]], cbind(add_ann = tmp$add_ann[[3L]],
                                  pk_ann = tmp$pk_ann[[3L]]))

    res <- peaksData(be, "pk_ann")
    expect_equal(res[[1L]], cbind(pk_ann = tmp$pk_ann[[1L]]))
    expect_equal(res[[2L]], cbind(pk_ann = tmp$pk_ann[[2L]]))
    expect_equal(res[[3L]], cbind(pk_ann = tmp$pk_ann[[3L]]))
})

test_that("peaksData<-,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    peaksData(be) <- list()

    be <- backendInitialize(be, test_df)
    lst <- list(cbind(mz = 1, intensity = 10.1),
                cbind(mz = 1:3, intensity = c(12.1, 12.4, 12.4)),
                cbind(other = 3.1, intensity = 100))
    expect_error(peaksData(be) <- "4", "list-like")
    expect_error(peaksData(be) <- list(1:3), "match length")
    expect_error(peaksData(be) <- list("a", "b", "c"), "'matrix'")
    expect_error(peaksData(be) <- lst, "same column names")
    lst <- list(cbind(mz = 1, intensity = 10.1),
                cbind(mz = 1:3, intensity = c(12.1, 12.4, 12.4)),
                cbind(mz = 3.1, intensity = 100))
    peaksData(be) <- lst
    expect_equal(peaksData(be), lst)

    lst2 <- list(cbind(intensity = 10.1, mz = 1),
                cbind(intensity = c(12.1, 12.4, 12.4), mz = 1:3),
                cbind(intensity = 100, mz = 3.1))
    peaksData(be) <- lst2
    expect_equal(peaksData(be), lst)

    lst2 <- list(cbind(intensity = 10.1, mz = 1, add_col = 15),
                cbind(intensity = c(12.1, 12.4, 12.4), mz = 1:3, add_col = 5:7),
                cbind(intensity = 100, mz = 3.1, add_col = 100))
    peaksData(be) <- lst2
    expect_equal(peaksData(be), lst)
    expect_equal(be@peaksDataFrame, list(data.frame(add_col = 15),
                                         data.frame(add_col = 5:7),
                                         data.frame(add_col = 100)))
})

test_that("$,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_equal(be$msLevel, integer())
    expect_equal(be$rtime, numeric())
    expect_error(be$other, "are not available")

    be <- backendInitialize(be, test_df)
    expect_equal(be$msLevel, test_df$msLevel)
    expect_equal(be$scanIndex, test_df$scanIndex)
    expect_equal(be$mz, IRanges::NumericList(test_df$mz, compress = FALSE))
    expect_equal(be$intensity, IRanges::NumericList(test_df$intensity,
                                                    compress = FALSE))

    tmp <- test_df
    tmp$peak_ann <- list(letters[1:3], letters[1:2], letters[1:4])
    be <- backendInitialize(be, tmp)
    expect_equal(be$peak_ann, tmp$peak_ann)
})

test_that("$<-,MsBackendMemory works", {
    be <- new("MsBackendMemory")

    expect_error(be$rtime <- 1.3, "data has 0")
    expect_error(be$mz <- list(1:3), "has to match")

    be <- backendInitialize(be, test_df)
    expect_error(be$rtime <- 1:2, "data has")

    ## existing spectra variable
    be$msLevel <- 1:3
    expect_equal(be$msLevel, 1:3)

    ## new spectra variable
    be$rtime <- c(4.2, 1.23, 4.23)
    expect_equal(be$rtime, c(4.2, 1.23, 4.23))
    be$new_var <- letters[1:3]
    expect_equal(be$new_var, letters[1:3])

    ## mz
    expect_error(be$mz <- list(1:2, 1:3, 1:4), "number of values")
    be$mz <- list(c(2.2, 2.3, 2.4), c(1.1, 1.2), c(1.1, 1.2, 1.3, 1.4))
    expect_equal(be$mz, IRanges::NumericList(list(c(2.2, 2.3, 2.4),
                                                  c(1.1, 1.2),
                                                  c(1.1, 1.2, 1.3, 1.4)),
                                             compress = FALSE))
    ## intensity
    vals <- list(c(12.2, 12.3, 12.4), c(12.2, 12.3), c(12.2, 12.3, 12.4, 12.5))
    expect_error(be$intensity <- list(1:2, 1:3, 1:4), "number of values")
    be$intensity <- vals
    expect_equal(be$intensity, IRanges::NumericList(vals, compress = FALSE))

    ## new peaks annotation
    vals <- list(letters[1:3], letters[1:2], letters[1:4])
    be$peak_anno <- vals
    expect_true(length(be@peaksDataFrame) == 3)
    expect_equal(be$peak_anno, vals)

    ## add to peaks annotation
    vals <- list(1:3, 1:2, 1:4)
    be$add_anno <- vals
    expect_equal(colnames(be@peaksDataFrame[[1L]]), c("peak_anno", "add_anno"))
    expect_equal(be$add_anno, vals)

    ## replace peaks annotation
    vals <- list(2:4, 2:3, 2:5)
    be$add_anno <- vals
    expect_equal(be$add_anno, vals)
})

test_that("[,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    df <- data.frame(scanIndex = 1:2, a = "a", b = "b")
    be <- backendInitialize(be, df)
    res <- be[1]
    expect_true(validObject(res))
    expect_equal(be@spectraData[1, ], res@spectraData[1, ])
    res <- be[2]
    expect_true(validObject(res))
    expect_equal(be@spectraData[2, ], res@spectraData[1, ])
    res <- be[2:1]
    expect_true(validObject(res))
    expect_equal(be@spectraData[2:1, ], res@spectraData)

    res <- be[c(FALSE, FALSE)]
    expect_true(validObject(res))
    expect_true(length(res) == 0)
    res <- be[c(FALSE, TRUE)]
    expect_true(validObject(res))
    expect_equal(be@spectraData[2, ], res@spectraData[1, ])

    expect_error(be[TRUE], "match the length of")
    expect_error(be["a"], "names")

    df <- data.frame(scanIndex = c(1L, 2L, 1L, 2L),
                     file = c("a", "a", "b", "b"))
    be <- backendInitialize(be, df)
    dataStorage(be) <- c("1", "1", "2", "2")
    res <- be[3]
    expect_true(validObject(res))
    expect_equal(dataStorage(res), "2")
    expect_equal(res@spectraData$file, "b")

    res <- be[c(3, 1)]
    expect_true(validObject(res))
    expect_equal(dataStorage(res), c("2", "1"))
    expect_equal(res@spectraData$file, c("b", "a"))
})

test_that("split,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    be <- backendInitialize(be, test_df)
    f <- factor(c("b", "a", "a"))
    res <- split(be, f)
    expect_true(is.list(res))
    expect_true(length(res) == 2)
    expect_s4_class(res[[1L]], "MsBackendMemory")
    expect_s4_class(res[[2L]], "MsBackendMemory")
    expect_equal(res[[1L]]$scanIndex, c(5, 6))
    expect_equal(res[[2L]]$scanIndex, c(4))

    res <- split(be, factor(c("b", "a", "a"), levels = c("b", "a")))
    expect_true(is.list(res))
    expect_true(length(res) == 2)
    expect_s4_class(res[[1L]], "MsBackendMemory")
    expect_s4_class(res[[2L]], "MsBackendMemory")
    expect_equal(res[[2L]]$scanIndex, c(5, 6))
    expect_equal(res[[1L]]$scanIndex, c(4))
})

test_that("filterAcquisitionNum,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_equal(be, filterAcquisitionNum(be, n = 4))

    df <- data.frame(acquisitionNum = c(1L, 2L, 3L, 2L, 3L, 1L, 2L, 4L),
                     msLevel = 1L)
    be <- backendInitialize(be, df)
    be$dataStorage <- c("1", "1", "1", "2", "2", "3", "3", "3")
    res <- filterAcquisitionNum(be, n = c(2L, 4L))
    expect_equal(length(res), 4)
    expect_equal(dataStorage(res), c("1", "2", "3", "3"))
    expect_equal(acquisitionNum(res), c(2L, 2L, 2L, 4L))

    res <- filterAcquisitionNum(be, n = 2L, dataStorage = "2")
    expect_equal(dataStorage(res), c("1", "1", "1", "2", "3", "3", "3"))
    expect_equal(acquisitionNum(res), c(1L, 2L, 3L, 2L, 1L, 2L, 4L))

    expect_error(filterAcquisitionNum(be, n = "a"), "integer representing")
    expect_equal(filterAcquisitionNum(be), be)
})

test_that("dataOrigin,dataOrigin<-,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_equal(dataOrigin(be), character())

    be <- backendInitialize(be, test_df)
    expect_equal(dataOrigin(be), c(NA_character_, NA_character_, NA_character_))

    dataOrigin(be) <- c("A", "B", "C")
    expect_equal(dataOrigin(be), c("A", "B", "C"))

    expect_error(dataOrigin(be) <- c("c", "d"), "length 3")
    expect_error(dataOrigin(be) <- 1:3, "length 3")

    expect_error(be$dataOrigin <- 1:3, "invalid")
    be$dataOrigin <- as.character(1:3)
    expect_equal(dataOrigin(be), as.character(1:3))
})

test_that("collisionEnergy,collisionEnergy<-,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_equal(collisionEnergy(be), numeric())

    be <- backendInitialize(be, test_df)
    expect_equal(collisionEnergy(be), rep(NA_real_, 3))
    collisionEnergy(be) <- c(1.2, 1.4, 1.3)
    expect_equal(collisionEnergy(be), c(1.2, 1.4, 1.3))
    expect_equal(be$collisionEnergy, c(1.2, 1.4, 1.3))

    expect_error(collisionEnergy(be) <- 3.2, "length 3")
    expect_error(collisionEnergy(be) <- c("a", "b", "c"), "length 3")
})

test_that("acquisitionNum,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_equal(acquisitionNum(be), integer())

    be <- backendInitialize(be, test_df)
    expect_equal(acquisitionNum(be), rep(NA_integer_, 3))
    be$acquisitionNum <- 1:3
    expect_equal(acquisitionNum(be), 1:3)
})

test_that("centroided,centroided<-,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_equal(centroided(be), logical())

    be <- backendInitialize(be, test_df)
    expect_equal(centroided(be), rep(NA, 3))
    centroided(be) <- c(TRUE, FALSE, TRUE)
    expect_equal(centroided(be), c(TRUE, FALSE, TRUE))
    expect_equal(be$centroided, c(TRUE, FALSE, TRUE))

    expect_error(centroided(be) <- c("a", "b", "c"), "logical")
    expect_error(centroided(be) <- c(TRUE, FALSE), "logical")

    centroided(be) <- FALSE
    expect_equal(centroided(be), rep(FALSE, 3))
})

test_that("ionCount,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_equal(ionCount(be), numeric())

    be <- backendInitialize(be, test_df)
    expect_equal(ionCount(be),
                 vapply(test_df$intensity, sum, numeric(1)))
})

test_that("isEmpty,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_identical(isEmpty(be), logical(0))

    be <- backendInitialize(be, test_df)
    expect_identical(isEmpty(be), c(FALSE, FALSE, FALSE))
})

test_that("isolationWindowLowerMz,isolationWindowLowerMz<-,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_identical(isolationWindowLowerMz(be), numeric())

    be <- backendInitialize(be, test_df)
    expect_identical(isolationWindowLowerMz(be), rep(NA_real_, 3))
    isolationWindowLowerMz(be) <- c(1.3, 1.4, 1.5)
    expect_identical(isolationWindowLowerMz(be), c(1.3, 1.4, 1.5))

    expect_error(isolationWindowLowerMz(be) <- 1.3, "length 3")
    expect_error(isolationWindowLowerMz(be) <- letters[1:3], "numeric")
})

test_that("isolationWindowTargetMz,isolationWindowTargetMz<-,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_identical(isolationWindowTargetMz(be), numeric())

    be <- backendInitialize(be, test_df)
    expect_identical(isolationWindowTargetMz(be), rep(NA_real_, 3))
    isolationWindowTargetMz(be) <- c(1.3, 1.4, 1.5)
    expect_identical(isolationWindowTargetMz(be), c(1.3, 1.4, 1.5))

    expect_error(isolationWindowTargetMz(be) <- 1.3, "length 3")
    expect_error(isolationWindowTargetMz(be) <- letters[1:3], "numeric")
})

test_that("isolationWindowUpperMz,isolationWindowUpperMz<-,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_identical(isolationWindowUpperMz(be), numeric())

    be <- backendInitialize(be, test_df)
    expect_identical(isolationWindowUpperMz(be), rep(NA_real_, 3))
    isolationWindowUpperMz(be) <- c(1.3, 1.4, 1.5)
    expect_identical(isolationWindowUpperMz(be), c(1.3, 1.4, 1.5))

    expect_error(isolationWindowUpperMz(be) <- 1.3, "length 3")
    expect_error(isolationWindowUpperMz(be) <- letters[1:3], "numeric")
})

test_that("msLevel,msLevel<-,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_identical(msLevel(be), integer())

    be <- backendInitialize(be, test_df)
    expect_identical(msLevel(be), test_df$msLevel)
    expect_identical(be$msLevel, test_df$msLevel)
    msLevel(be) <- 3:1
    expect_identical(msLevel(be), 3:1)
    expect_identical(be$msLevel, 3:1)

    be$msLevel <- 1:3
    expect_identical(msLevel(be), 1:3)
    expect_identical(be$msLevel, 1:3)

    expect_error(msLevel(be) <- 1L, "length 3")
    expect_error(msLevel(be) <- letters[1:3], "integer")
})

test_that("polarity,polarity<-,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_identical(polarity(be), integer())

    be <- backendInitialize(be, test_df)
    expect_identical(polarity(be), rep(NA_integer_, 3))
    expect_identical(be$polarity, rep(NA_integer_, 3))
    polarity(be) <- 3:1
    expect_identical(polarity(be), 3:1)
    expect_identical(be$polarity, 3:1)

    be$polarity <- 1:3
    expect_identical(polarity(be), 1:3)
    expect_identical(be$polarity, 1:3)

    polarity(be) <- 1L
    expect_identical(polarity(be), rep(1L, 3))
    expect_error(polarity(be) <- letters[1:3], "integer")
})

test_that("precScanNum,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_identical(precScanNum(be), integer())

    be <- backendInitialize(be, test_df)
    expect_identical(precScanNum(be), rep(NA_integer_, 3))
    expect_identical(be$precScanNum, rep(NA_integer_, 3))
    be$precScanNum <- 1:3
    expect_identical(precScanNum(be), 1:3)
    expect_identical(be$precScanNum, 1:3)
})

test_that("precursorCharge,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_identical(precursorCharge(be), integer())

    be <- backendInitialize(be, test_df)
    expect_identical(precursorCharge(be), rep(NA_integer_, 3))
    expect_identical(be$precursorCharge, rep(NA_integer_, 3))
    be$precursorCharge <- 1:3
    expect_identical(precursorCharge(be), 1:3)
    expect_identical(be$precursorCharge, 1:3)
})

test_that("precursorIntensity,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_identical(precursorIntensity(be), numeric())

    be <- backendInitialize(be, test_df)
    expect_identical(precursorIntensity(be), rep(NA_real_, 3))
    expect_identical(be$precursorIntensity, rep(NA_real_, 3))
    be$precursorIntensity <- c(12.2, 12.5, 15.2)
    expect_identical(precursorIntensity(be), c(12.2, 12.5, 15.2))
    expect_identical(be$precursorIntensity, c(12.2, 12.5, 15.2))
})

test_that("precursorMz,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_identical(precursorMz(be), numeric())

    be <- backendInitialize(be, test_df)
    expect_identical(precursorMz(be), rep(NA_real_, 3))
    expect_identical(be$precursorMz, rep(NA_real_, 3))
    be$precursorMz <- c(1.2, 1.3, 1.4)
    expect_identical(precursorMz(be), c(1.2, 1.3, 1.4))
    expect_identical(be$precursorMz, c(1.2, 1.3, 1.4))
})

test_that("rtime,rtime<-,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_identical(rtime(be), numeric())

    be <- backendInitialize(be, test_df)
    expect_identical(rtime(be), rep(NA_real_, 3))
    expect_identical(be$rtime, rep(NA_real_, 3))
    rtime(be) <- c(1.2, 1.3, 1.4)
    expect_identical(rtime(be), c(1.2, 1.3, 1.4))
    expect_identical(be$rtime, c(1.2, 1.3, 1.4))

    be$rtime <- c(2.1, 2.2, 2.3)
    expect_identical(rtime(be), c(2.1, 2.2, 2.3))
    expect_identical(be$rtime, c(2.1, 2.2, 2.3))

    expect_error(rtime(be) <- 1.3, "length 3")
    expect_error(rtime(be) <- letters[1:3], "numeric")
})

test_that("scanIndex,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_identical(scanIndex(be), integer())

    be <- backendInitialize(be, test_df)
    expect_identical(scanIndex(be), test_df$scanIndex)
    expect_identical(be$scanIndex, test_df$scanIndex)
    be$scanIndex <- 1:3
    expect_identical(scanIndex(be), 1:3)
    expect_identical(be$scanIndex, 1:3)
})

test_that("smoothed,smoothed<-,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_identical(smoothed(be), logical())

    be <- backendInitialize(be, test_df)
    expect_identical(smoothed(be), rep(NA, 3))
    expect_identical(be$smoothed, rep(NA, 3))
    smoothed(be) <- c(TRUE, FALSE, TRUE)
    expect_identical(smoothed(be), c(TRUE, FALSE, TRUE))
    expect_identical(be$smoothed, c(TRUE, FALSE, TRUE))

    smoothed(be) <- FALSE
    expect_identical(smoothed(be), rep(FALSE, 3))
    expect_error(smoothed(be) <- letters[1:3], "logical")
})

test_that("spectraNames,spectraNames<-,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    expect_identical(spectraNames(be), character())

    be <- backendInitialize(be, test_df)
    expect_identical(spectraNames(be), as.character(1:3))
    spectraNames(be) <- c("a", "b", "c")
    expect_identical(spectraNames(be), c("a", "b", "c"))

    expect_error(spectraNames(be) <- 1, "length")
})

test_that("backendMerge,MsBackendMemory works", {
    be <- new("MsBackendMemory")
    tmp <- list(be, be, be)
    res <- backendMerge(tmp)
    expect_equal(be, res)

    ## With empty one in between.
    be2 <- backendInitialize(be, test_df)
    res <- backendMerge(list(be2, be, be2))
    expect_s4_class(res, "MsBackendMemory")
    expect_true(length(res) == 6)
    expect_equal(res[1:3], be2)
    expect_equal(res$msLevel[4:6], be2$msLevel)
    expect_equal(res$mz[4:6], be2$mz)
    expect_equal(res$intensity[4:6], be2$intensity)
})

test_that("selectSpectraVariables,MsBackendMemory works", {
    be <- MsBackendMemory()
    expect_error(selectSpectraVariables(be, c("msLevel", "other")),
                 "not available")

    be <- backendInitialize(be, test_df)
    expect_error(selectSpectraVariables(be, c("msLevel")), "dataStorage")
    res <- selectSpectraVariables(be, c("msLevel", "dataStorage"))
    expect_equal(colnames(res@spectraData), c("msLevel", "dataStorage"))
    expect_equal(peaksVariables(res), NULL)

    res <- selectSpectraVariables(be, c("mz", "dataStorage"))
    expect_equal(colnames(res@spectraData), c("dataStorage"))
    expect_equal(peaksVariables(res), "mz")

    be$peak_anno <- list(c("a", "", "b"), c("a", "b"), c("a", "b", "c", "d"))
    res <- selectSpectraVariables(be, c("msLevel", "dataStorage",
                                        "mz", "intensity"))
    expect_equal(res@peaksDataFrame, list())

    be$second_col <- list(1:3, 1:2, 1:4)
    expect_equal(peaksVariables(be), c("mz", "intensity",
                                       "peak_anno", "second_col"))
    res <- selectSpectraVariables(be, c("msLevel", "dataStorage",
                                        "mz", "intensity", "second_col"))
    expect_equal(sort(peaksVariables(res)),
                 sort(c("mz", "intensity", "second_col")))
})

test_that("tic,MsBackendMemory works", {
    be <- MsBackendMemory()
    expect_equal(tic(be), numeric())

    be <- backendInitialize(be, test_df)
    expect_equal(tic(be), rep(NA_real_, 3))
    expect_equal(tic(be, initial = FALSE)[1], sum(intensity(be)[[1L]]))
    expect_equal(tic(be, initial = FALSE)[2], sum(intensity(be)[[2L]]))
})
