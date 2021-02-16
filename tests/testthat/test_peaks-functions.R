test_that(".peaks_replace_intensity works", {
    int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 1, 2, 0,
             0, 1, 5, 10, 5, 1)
    x <- cbind(mz = 1:length(int), intensity = int)
    res <- .peaks_replace_intensity(x, 1L, centroided = FALSE)
    expect_equal(res, x)
    res <- .peaks_replace_intensity(x, 1L, centroided = FALSE, threshold = 3)
    int_exp <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 0, 0, 0, 0,
                 1, 5, 10, 5, 1)
    expect_equal(res[, "intensity"], int_exp)
    expect_warning(res <- .peaks_replace_intensity(x, 1L, threshold = max, value = 3))
    expect_equal(res, x)
    res <- .peaks_replace_intensity(x, 1L, threshold = max,
                                    value = 3, centroided = TRUE)
    expect_true(all(res[, "intensity"] == 3))

    res <- .peaks_replace_intensity(x, 1L, function(z, ...) min(z[z > 0]),
                                    centroided = TRUE)
    int_exp <- c(0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 3, 10, 6, 2, 0, 0, 0, 2, 0, 0,
                 0, 5, 10, 5, 0)
    expect_equal(res[, "intensity"], int_exp)

    res <- .peaks_replace_intensity(x, 1L, centroided = TRUE, msLevel = 2L)
    expect_equal(res[, "intensity"], int)

    ## Working with NAs.
    x[4, 2] <- NA
    res <- .peaks_replace_intensity(x, 1L, centroided = TRUE, threshold = 4)
    expect_equal(res[, 2], c(0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 10, 6,
                             0, 0, 0, 0, 0, 0, 0, 0, 5, 10, 5, 0))
})

test_that(".peaks_filter_intensity works", {
    ints <- c(5, 3, 12, 14.4, 13.3, 9, 3, 0, NA, 21, 89, 55, 33, 5, 2)
    x <- cbind(mz = seq_along(ints), intensity = ints)
    res <- .peaks_filter_intensity(x, spectrumMsLevel = 1L)
    expect_equal(res[, 1], c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15))

    res <- .peaks_filter_intensity(x, spectrumMsLevel = 1L, msLevel = 2L)
    expect_equal(res, x)

    res <- .peaks_filter_intensity(x, spectrumMsLevel = 1L, intensity = c(5, 15))
    expect_equal(res[, 1], c(1, 3, 4, 5, 6, 14))
})

test_that(".peaks_filter_intensity_function works", {
    ints <- c(5, 3, 12, 14.4, 13.3, 9, 3, 0, NA, 21, 89, 55, 33, 5, 2)
    x <- cbind(mz = seq_along(ints), intensity = ints)
    res <- .peaks_filter_intensity_function(
        x, spectrumMsLevel = 1L, intfun = function(x) x > 0)
    expect_equal(res[, 1], c(1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15))

    res <- .peaks_filter_intensity_function(
        x, spectrumMsLevel = 1L,
        intfun = function(z) z > max(z, na.rm = TRUE) / 2)
    expect_true(all(res[, "intensity"] > 89/2))

    res <- .peaks_filter_intensity_function(x, spectrumMsLevel = 1L,
                                            msLevel = 2L)
    expect_equal(res, x)

    expect_error(.peaks_filter_intensity_function(
        x, spectrumMsLevel = 1L, function(x) FALSE),
        "does not return a")
    expect_error(.peaks_filter_intensity_function(
        x, spectrumMsLevel = 1L, function(x) which(is.na(x))),
        "does not return a")
})

test_that(".peaks_bin works", {
    int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 1, 2, 0,
             0, 1, 5, 10, 5, 1)
    x <- cbind(mz = 1:length(int), intensity = int)

    res <- .peaks_bin(x, spectrumMsLevel = 1L)
    expect_identical(res[, 2], x[, 2])
    expect_identical(res[, 1], x[, 1] + 0.5)

    res <- .peaks_bin(x, spectrumMsLevel = 1L, binSize = 2L)
    expect_equal(res[, 2], c(1, 5, 1, 0, 1, 13, 8, 1, 3, 0, 6, 15, 1))
    res <- .peaks_bin(x, spectrumMsLevel = 1L, binSize = 2L, FUN = max)
    expect_equal(res[, 2], c(1, 3, 1, 0, 1, 10, 6, 1, 2, 0, 5, 10, 1))
    res <- .peaks_bin(x, spectrumMsLevel = 1L, msLevel = 2L, binSize = 2L)
    expect_identical(res, x)
})

test_that("joinPeaks works", {
    x <- cbind(c(31.34, 50.14, 60.3, 120.9, 230, 514.13, 874.1),
               1:7)
    y <- cbind(c(12, 31.35, 70.3, 120.9 + ppm(120.9, 5),
                 230 + ppm(230, 10), 315, 514.14, 901, 1202),
               1:9)

    res <- joinPeaks(x, y, ppm = 0)
    expect_true(nrow(res$x) == nrow(res$y))
    expect_true(nrow(res$x) == nrow(x) + nrow(y))
    res <- joinPeaks(x, y, ppm = 0, type = "inner")
    expect_true(nrow(res$x) == nrow(res$y))
    expect_true(nrow(res$x) == 0)

    ## ppm 5
    res <- joinPeaks(x, y, ppm = 5)
    expect_true(nrow(res$x) == nrow(res$y))
    expect_true(nrow(res$x) == nrow(x) + nrow(y) - 1)
    expect_true(!is.na(res$x[7, 1]))
    expect_true(!is.na(res$y[7, 1]))
    res <- joinPeaks(x, y, ppm = 5, type = "inner")
    expect_true(nrow(res$x) == nrow(res$y))

    ## ppm 10
    res <- joinPeaks(x, y, ppm = 10)
    expect_true(nrow(res$x) == nrow(res$y))
    expect_true(nrow(res$x) == nrow(x) + nrow(y) - 2)
    res <- joinPeaks(x, y, ppm = 10, type = "inner")
    expect_true(nrow(res$x) == nrow(res$y))
    expect_true(nrow(res$x) == 2)

    ## tolerance 0.01
    res <- joinPeaks(x, y, tolerance = 0.01)
    expect_true(nrow(res$x) == nrow(res$y))
    expect_true(nrow(res$x) == nrow(x) + nrow(y) - 4)
    res <- joinPeaks(x, y, tolerance = 0.01, type = "left")
    expect_true(nrow(res$x) == nrow(res$y))
    expect_equal(res$x, x)
    res <- joinPeaks(x, y, tolerance = 0.01, type = "right")
    expect_true(nrow(res$x) == nrow(res$y))
    expect_equal(res$y, y)
})

test_that(".peaks_pick works", {
    e <- matrix(NA_real_, nrow = 0, ncol = 2,
                dimnames = list(c(), c("mz", "intensity")))
    expect_warning(.peaks_pick(e, spectrumMsLevel = 1), "empty")
    expect_equal(suppressWarnings(.peaks_pick(e, spectrumMsLevel = 1)), e)

    int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 1, 2, 0,
             0, 1, 5, 10, 5, 1)
    x <- cbind(mz = seq_along(int), intensity = int)
    expect_equal(.peaks_pick(x, spectrumMsLevel = 2, msLevel = 1), x)
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 2,
                             snr = 1),
                 cbind(mz = c(4, 12, 18, 23), intensity = c(3, 10, 2, 10)))
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 2,
                             snr = 5),
                 cbind(mz = c(12, 23), intensity = c(10, 10)))
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 2,
                             snr = 10), e)
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 10,
                             snr = 1),
                 cbind(mz = c(12, 23), intensity = c(10, 10)))
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 2,
                             method = "SuperSmoother", snr = 2),
                 cbind(mz = c(4, 12, 23), intensity = c(3, 10, 10)))

    ## use refineCentroids
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 2,
                             method = "MAD", snr = 2, k = 1L, ),
                 cbind(mz = mapply(weighted.mean,
                                   x = list(3:5, 11:13, 22:24),
                                   w = list(int[3:5], int[11:13], int[22:24])),
                       intensity = c(3, 10, 10)))
})

test_that(".peaks_smooth works", {
    e <- matrix(NA_real_, nrow = 0, ncol = 2,
                dimnames = list(c(), c("mz", "intensity")))
    expect_warning(.peaks_smooth(e, spectrumMsLevel = 1), "empty")
    expect_equal(suppressWarnings(.peaks_smooth(e, spectrumMsLevel = 1)), e)

    cf <- matrix(0.2, nrow = 5, ncol = 5)
    int <- 1:10
    x <- cbind(mz = seq_along(int), intensity = int)
    expect_equal(.peaks_smooth(x, spectrumMsLevel = 1, coef = cf),
                 cbind(mz = x[, 1L], intensity = rep(3:8, c(3, 1, 1, 1, 1, 3))))
})

test_that(".peaks_filter_mz_range works", {
    p <- cbind(mz = c(2, 5.6, 123.2, 422.8, 599.3, 599.4, 599.5, 743.1),
               intensity = 1:8)
    expect_warning(res <- .peaks_filter_mz_range(
                       p, 1L, mz = range(numeric())), "Inf")
    expect_equal(res, p)

    res <- .peaks_filter_mz_range(p, 1L, msLevel = 2)
    expect_equal(res, p)

    res <- .peaks_filter_mz_range(p, 1L, mz = c(200, 600))
    expect_equal(res[, "intensity"], c(4, 5, 6, 7))
})

test_that(".peaks_match_mz_value works", {
    p <- cbind(mz = c(2, 5.6, 123.2, 422.8, 599.3, 599.4, 599.5, 743.1),
               intensity = 1:8)
    res <- .peaks_filter_mz_value(p, 1L, mz = numeric())
    expect_true(is.matrix(res))
    expect_true(nrow(res) == 0)
    expect_equal(colnames(res), c("mz", "intensity"))
    res <- .peaks_filter_mz_value(p, 1L, mz = NA)
    expect_true(is.matrix(res))
    expect_equal(colnames(res), c("mz", "intensity"))
    expect_true(nrow(res) == 0)

    res <- .peaks_filter_mz_value(p, 1L, mz = 5, tolerance = 1)
    expect_equal(unname(res[, "intensity"]), 2)
    res <- .peaks_filter_mz_value(p, 1L, mz = c(5.5, 599.41), tolerance = 0.1)
    expect_equal(unname(res[, "intensity"]), c(2, 6))

    res <- .peaks_filter_mz_value(p, 1L, mz = c(123, 742.2),
                                  tolerance = c(0.2, 1))
    expect_equal(unname(res[, "intensity"]), c(3, 8))
})

test_that("joinPeaksGnps works", {
    a <- cbind(mz = c(10, 36, 63, 91, 93),
               intensity = c(14, 15, 999, 650, 1))
    a_pmz <- 91

    b <- cbind(mz = c(10, 12, 50, 63, 105),
               intensity = c(35, 5, 16, 999, 450))
    b_pmz <- 105

    expect_equal(joinPeaksGnps(a, b), joinPeaks(a, b))
    expect_equal(joinPeaksGnps(a, b, xPrecursorMz = 91),
                 joinPeaks(a, b))
    expect_equal(joinPeaksGnps(a, b, xPrecursorMz = 3, yPrecursorMz = 3),
                 joinPeaks(a, b))

    res <- joinPeaksGnps(a, b, xPrecursorMz = a_pmz, yPrecursorMz = b_pmz,
                         type = "left")
    expect_true(nrow(res$x) > nrow(a))  # peaks can match multiple
    ## peak 36 matches no peak in b directly, but matches peak 50 in b
    ## considering the precursor difference
    expect_true(sum(res$x[, 1] == 36) == 2)

    res_2 <- joinPeaksGnps(b, a, xPrecursorMz = b_pmz, yPrecursorMz = a_pmz,
                           type = "left")
    expect_equal(res_2$x[, 1], c(10, 12, 50, 50, 63, 105, 105))

    ## example with multi matches.
    a <- cbind(mz = c(10, 26, 36, 63, 91, 93),
               intensity = c(14, 10, 15, 999, 650, 1))
    a_pmz <- 91

    b <- cbind(mz = c(10, 12, 24, 50, 63, 105),
               intensity = c(35, 10, 5, 16, 999, 450))
    b_pmz <- 105

    res <- joinPeaksGnps(a, b, xPrecursorMz = a_pmz, yPrecursorMz = b_pmz,
                         type = "left")
    expect_equal(res$x[, 1], c(10, 10, 26, 36, 36, 63, 91, 91, 93))
    expect_equal(res$y[, 1], c(10, 24, NA, NA, 50, 63, NA, 105, NA))

    ## Example with peaks from a matching multiple peaks in b and vice versa
    a <- cbind(mz = c(10, 12, 14, 16, 19),
               intensity = 1:5)
    a_pmz <- 4
    b <- cbind(mz = c(10, 14, 16, 17, 19, 22),
               intensity = 1:6)
    b_pmz <- 8

    exp_a <- c(10, 10, 12, 12, 14, 16, 19, NA, NA)
    exp_b <- c(10, 14, NA, 16, 14, 16, 19, 17, 22)

    res <- joinPeaksGnps(a, b, a_pmz, b_pmz, type = "outer")
    expect_equal(res$x[, 1], exp_a)
    expect_equal(res$y[, 1], exp_b)
})

test_that("gnps works", {
    x <- cbind(mz = c(10, 12, 14, 16, 19),
               intensity = 1:5)
    xpmz <- 4
    y <- cbind(mz = c(10, 14, 16, 17, 19, 22),
               intensity = 1:6)
    ypmz <- 8

    expect_error(gnps(x, y), "aligned")
    a <- cbind(c(1, 2), c(0, 0))
    expect_equal(gnps(a, a), 0)
    a <- rbind(c(NA, NA), a)
    b <- a
    b[2:3, 1] <- NA
    expect_equal(gnps(a, b), 0)

    map <- joinPeaksGnps(x, y, xpmz, ypmz)
    res <- gnps(map$x, map$y)
    expect_true(res > 0.5)
    map <- joinPeaksGnps(y, x, ypmz, xpmz)
    res_2 <- gnps(map$x, map$y)
    expect_equal(res, res_2)

    ## Unit tests with reference values from GNPS online kindly provided by
    ## Liesa Salzer
    apmz <- 184.097
    a <- cbind(
        mz = c(59.049, 82.065002, 100.075996, 110.059998, 124.075996,
               125.059998, 142.085999),
        intensity = c(346.036011, 705.127014, 752.387024, 2193.239014,
                      2194.920898, 2455.741943, 3176.063965))
    bpmz <- 202.108
    b <- cbind(
        mz = c(55.018002, 59.049, 82.065002, 84.044998, 100.075996, 110.059998,
               125.059998, 142.085999, 152.070999, 156.065002),
        intensity = c(163.115997, 788.68103, 345.580994, 298.221008, 556.021973,
                      1938.749023, 1202.94104, 1801.137939, 1728.411011,
                      77.516998))
    map <- joinPeaksGnps(a, b, apmz, bpmz, tolerance = 0.1, ppm = 0)
    expect_equal(round(gnps(map$x, map$y), 2), 0.76)
    mapr <- joinPeaksGnps(b, a, bpmz, apmz, tolerance = 0.1, ppm = 0)
    expect_equal(gnps(map$x, map$y), gnps(mapr$x, mapr$y))

    dpmz <- 220.12
    d <- cbind(
        mz = c(58.064999, 59.049, 60.060001, 82.065002, 100.075996, 110.059998,
               125.059998, 142.085999, 152.070999, 161.059998),
        intensity = c(49.41, 319.764008, 76.835999, 180.686005, 278.931,
                      969.950012, 653.91803, 1037.552979, 1123.036987,
                      245.953003))
    map <- joinPeaksGnps(a, d, apmz, dpmz, tolerance = 0.1, ppm = 0)
    expect_equal(round(gnps(map$x, map$y), 2), 0.74)
    mapr <- joinPeaksGnps(d, a, dpmz, apmz, tolerance = 0.1, ppm = 0)
    expect_equal(gnps(map$x, map$y), gnps(mapr$x, mapr$y))

    map <- joinPeaksGnps(b, d, bpmz, dpmz, tolerance = 0.1, ppm = 0)
    expect_equal(round(gnps(map$x, map$y), 2), 0.93)
    mapr <- joinPeaksGnps(d, b, dpmz, bpmz, tolerance = 0.1, ppm = 0)
    expect_equal(gnps(map$x, map$y), gnps(mapr$x, mapr$y))

    ## Second set.
    apmz <- 371.564
    a <- cbind(
        mz = c(73.028999, 80.055, 87.043999, 89.059998, 103.039001, 133.085999,
               147.065002, 177.112, 195.121994, 207.123001, 221.138, 239.149002,
               279.144012, 283.174988, 309.156006, 327.201996, 354.225006,
               415.257996, 425.276001, 443.290009, 459.277008, 469.303986,
               487.311005, 503.30899, 531.336975, 547.328979, 575.362976,
               576.367004),
        intensity = c(324.171997, 29.08, 179.298996, 1954.572021, 194.485992,
                      1176.598999, 31.034, 316.217987, 19.666, 21.794001,
                      55.692001, 40.776001, 27.982, 42.688999, 48.007999,
                      145.162994, 34.383999, 0.624, 2.032, 1.357, 3.382, 1.338,
                      5.238, 4.136, 8.294, 3.376, 6.845, 0.766)
    )
    bpmz <- 489.659
    b <- cbind(
        mz = c(73.028, 87.043999, 89.059998, 99.07, 103.039001, 133.085999,
               143.106995, 147.065002, 177.112, 199.125, 221.138, 231.158997,
               235.117996, 283.175995, 293.196014, 327.201996, 337.223999,
               353.21701, 371.227997, 397.242004, 425.171997, 427.253998,
               454.276001, 455.248993, 472.286011, 529.312012, 547.328979,
               548.344971, 591.356995, 608.393982, 635.382996, 652.403015,
               679.411987, 696.434021, 723.435974, 740.45697, 767.45697,
               784.479004, 811.471985),
        intensity = c(121.129997, 394.515991, 1239.718994, 206.699997,
                      427.983002, 830.174988, 34.665001, 20.660999, 277.040985,
                      16.813, 85.750999, 15.945, 11.498, 7.878, 4.811, 6.001,
                      2.187, 3.93, 1.947, 2.531, 16.639, 3.748, 7.272, 2.639,
                      7.758, 0.797, 1.22, 0.7, 2.127, 0.797, 2.697, 1.089,
                      3.461, 2.321, 2.603, 2.696, 2.262, 2.222, 0.866)
    )
    dpmz <- 459.772
    d <- cbind(
        mz = c(73.028, 87.043999, 89.059998, 99.073997, 103.039001, 133.085999,
               143.106995, 177.112, 187.132996, 195.121994, 221.139008,
               231.158997, 239.149002, 283.174988, 293.196014, 306.065002,
               327.200989, 371.22699, 378.049011, 395.075012, 413.084991,
               415.253998, 442.277008, 503.307007, 513.320984, 531.333008,
               557.346008, 575.366028, 601.372009, 619.390991, 645.396973,
               663.41803, 664.414001, 707.440002, 708.440979, 751.46698),
        intensity = c(278.23999, 305.428009, 1729.343994, 254.378998,
                      123.156998, 1135.667969, 103.471001, 305.052002,
                      48.874001, 38.331001, 46.912998, 30.42, 33.466999,
                      25.868999, 19.299999, 15.968, 35.019001, 17.221001,
                      18.408001, 22.473, 13.028, 40.317001, 19.250999, 1.299,
                      1.099, 0.404, 0.361, 1.844, 0.476, 2.804, 0.359, 4.112,
                      0.398, 3.889, 0.447, 2.446)
    )
    epmz <- 432.28
    e <- cbind(
        mz = c(73.065002, 74.096001, 87.043999, 89.059998, 130.158997,
               131.070007, 133.085999, 175.097, 177.112, 195.123001, 221.138,
               239.149002, 265.165009, 283.174988, 309.190002, 327.201996,
               337.165009, 353.213013, 371.227997, 379.234009),
        intensity = c(18.827, 35.034, 110.706001, 1887.61499, 62.372002, 16.854,
                      1189.276978, 9.96, 385.527008, 3.877, 73.189003,
                      21.197001, 13.633, 28.17, 11.297, 42.102001, 2.14, 4.764,
                      38.998001, 1.858)
    )
    fpmz <- 344.229
    f <- cbind(
        mz = c(59.049, 73.042999, 74.096001, 87.043999, 89.059998, 107.07,
               130.158997, 133.085999, 177.112, 195.123001, 221.138, 239.149002,
               265.165009, 280.321014, 283.174988, 309.19101, 326.197998),
        intensity = c(0.702, 12.998, 14.456, 35.036999, 988.768005, 13.287,
                      6.312, 635.674988, 178.367004, 11.298, 28.547001, 29.573,
                      1.358, 1.233, 30.393999, 2.984, 2.935)
    )
    map <- joinPeaksGnps(a, b, apmz, bpmz, tolerance = 0.1, ppm = 0)
    expect_equal(round(gnps(map$x, map$y), 2), 0.90)
    mapr <- joinPeaksGnps(b, a, bpmz, apmz, tolerance = 0.1, ppm = 0)
    expect_equal(gnps(map$x, map$y), gnps(mapr$x, mapr$y))
    map <- joinPeaksGnps(a, d, apmz, dpmz, tolerance = 0.1, ppm = 0)
    expect_equal(round(gnps(map$x, map$y), 2), 0.90)
    mapr <- joinPeaksGnps(d, a, dpmz, apmz, tolerance = 0.1, ppm = 0)
    expect_equal(gnps(map$x, map$y), gnps(mapr$x, mapr$y))
    map <- joinPeaksGnps(a, e, apmz, epmz, tolerance = 0.1, ppm = 0)
    expect_equal(round(gnps(map$x, map$y), 2), 0.91)
    mapr <- joinPeaksGnps(e, a, epmz, apmz, tolerance = 0.1, ppm = 0)
    expect_equal(gnps(map$x, map$y), gnps(mapr$x, mapr$y))
    map <- joinPeaksGnps(a, f, apmz, fpmz, tolerance = 0.1, ppm = 0)
    expect_equal(round(gnps(map$x, map$y), 2), 0.91)
    mapr <- joinPeaksGnps(f, a, fpmz, apmz, tolerance = 0.1, ppm = 0)
    expect_equal(gnps(map$x, map$y), gnps(mapr$x, mapr$y))

    ## Third set
    apmz <- 488.358
    a <- cbind(
        mz = c(57.068001, 74.096001, 87.043999, 89.059998, 121.065002,
               133.085999, 147.080002, 165.091003, 177.112, 209.117004,
               221.139008, 233.190002, 253.143005, 271.153992, 277.216003,
               315.181, 321.243011, 359.207001, 365.21701, 383.179993,
               401.166992, 419.177002, 427.303986, 442.194, 471.330994),
        intensity = c(7114.173828, 1001.093994, 770.18103, 13194.55957,
                      6228.618164, 8398.287109, 1492.104004, 5111.643066,
                      2547.698975, 1840.756958, 245.207001, 3336.004883,
                      313.119995, 626.439026, 4527.278809, 461.109985,
                      682.629028, 8239.036133, 43.747002, 16.875999, 240.197998,
                      276.912994, 65.282997, 31.399, 14.15)
    )
    bpmz <- 356.28
    b <- cbind(
        mz = c(57.07, 74.096001, 89.059998, 121.065002, 133.085999, 139.074997,
               147.080002, 165.091003, 209.117996, 227.128006, 233.190002,
               250.901001, 268.912994, 277.216003, 303.231995, 309.901001,
               315.194, 321.239014),
        intensity = c(7689.937012, 1007.257996, 4579.356934, 8001.039062,
                      2327.36792, 881.854004, 446.681, 6326.22998, 497.213013,
                      17089.107422, 153.132996, 17.469999, 53.928001, 93.917,
                      5.337, 11.661, 3.787, 73.805)
    )
    map <- joinPeaksGnps(a, b, apmz, bpmz, tolerance = 0.1, ppm = 0)
    expect_equal(round(gnps(map$x, map$y), 2), 0.86)
    mapr <- joinPeaksGnps(b, a, bpmz, apmz, tolerance = 0.1, ppm = 0)
    expect_equal(gnps(map$x, map$y), gnps(mapr$x, mapr$y))
})
