test_df <- DataFrame(msLevel = c(1L, 2L, 2L), scanIndex = 4:6)
test_df$mz <- list(c(1.1, 1.3, 1.5), c(4.1, 5.1), c(1.6, 1.7, 1.8, 1.9))
test_df$intensity <- list(c(45.1, 34, 12), c(234.4, 1333), c(42.1, 34.2, 65, 6))

test_that("backendInitialize,MsBackendDF works", {
    be <- new("MsBackendDF")
    expect_true(validObject(be))

    be <- backendInitialize(be, data = DataFrame(msLevel = 2L))
    expect_true(validObject(be))
    expect_equal(dataStorage(be), "<memory>")
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

test_that("show,MsBackendDF works", {
    be <- new("MsBackendDF")
    expect_output(show(be), "MsBackendDF")
})

test_that("dataStorage,dataStorage<-MsBackendDF works", {
    be <- new("MsBackendDF")
    expect_equal(dataStorage(be), character())

    be <- backendInitialize(be, test_df)
    expect_equal(dataStorage(be), rep("<memory>", 3))

    dataStorage(be) <- c("other", "storage",  "mode")
    expect_equal(dataStorage(be), c("other", "storage", "mode"))

    expect_error(dataStorage(be) <- 3, "length 3")
    expect_error(dataStorage(be) <- 1:3, "'character'")
})

test_that("length,MsBackendDF works", {
    be <- new("MsBackendDF")
    expect_true(length(be) == 0)
    be <- backendInitialize(be, test_df)
    expect_true(length(be) == 3)
})

test_that("spectraVariables,MsBackendDF works", {
    be <- new("MsBackendDF")
    res <- spectraVariables(be)
    expect_equal(res, names(coreSpectraVariables()))

    df <- data.frame(new_col = "a")
    df$mz <- list(1:3)
    be <- backendInitialize(be, df)
    res <- spectraVariables(be)
    expect_equal(res, c(names(coreSpectraVariables()), "new_col"))
})

test_that("peaksVariables,MsBackendDF works", {
    be <- new("MsBackendDF")
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

test_that("lengths,MsBackendDF works", {
    be <- new("MsBackendDF")
    expect_equal(lengths(be), integer())

    be <- backendInitialize(be, test_df)
    expect_equal(lengths(be), c(3L, 2L, 4L))
})

test_that("mz,mz<-,MsBackendDF works", {
    be <- new("MsBackendDF")
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

test_that("intensity,intensity<-,MsBackendDF works", {
    be <- new("MsBackendDF")
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

test_that("spectraData,MsBackendDF works", {
    be <- new("MsBackendDF")

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

test_that("spectraData<-,MsBackendDF works", {
    be <- backendInitialize(new("MsBackendDF"), test_df)

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

test_that("peaksData,MsBackendDF works", {
    be <- new("MsBackendDF")
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

test_that("peaksData<-,MsBackendDF works", {
    be <- new("MsBackendDF")
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

test_that("$,MsBackendDF works", {
    be <- new("MsBackendDF")
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

test_that("$<-,MsBackendDF works", {
    be <- new("MsBackendDF")

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

test_that("[,MsBackendDF works", {
    be <- new("MsBackendDF")
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

test_that("split,MsBackendDF works", {
    be <- new("MsBackendDF")
    be <- backendInitialize(be, test_df)
    f <- factor(c("b", "a", "a"))
    res <- split(be, f)
    expect_true(is.list(res))
    expect_true(length(res) == 2)
    expect_s4_class(res[[1L]], "MsBackendDF")
    expect_s4_class(res[[2L]], "MsBackendDF")
    expect_equal(res[[1L]]$scanIndex, c(5, 6))
    expect_equal(res[[2L]]$scanIndex, c(4))

    res <- split(be, factor(c("b", "a", "a"), levels = c("b", "a")))
    expect_true(is.list(res))
    expect_true(length(res) == 2)
    expect_s4_class(res[[1L]], "MsBackendDF")
    expect_s4_class(res[[2L]], "MsBackendDF")
    expect_equal(res[[2L]]$scanIndex, c(5, 6))
    expect_equal(res[[1L]]$scanIndex, c(4))
})
