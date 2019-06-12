#' @include hidden_aliases.R
NULL

#' @title The Spectra class to manage and access MS data
#'
#' @aliases Spectra-class [,Spectra-method
#'
#' @name Spectra
#'
#' @description
#'
#' The `Spectra` class encapsules spectral mass spectrometry data and
#' related metadata.
#'
#' It supports multiple data backends, e.g. in-memory ([MsBackendDataFrame()]),
#' on-disk as mzML ([MsBackendMzR()]).
#'
#' @details
#'
#' The `Spectra` class uses by default a lazy data manipulation strategy,
#' i.e. data manipulations such as performed with `removePeaks` are not applied
#' immediately to the data, but applied on-the-fly to the spectrum data once it
#' is retrieved.
#'
#' @section Creation of objects, conversion and changing the backend:
#'
#' `Spectra` classes can be created with the `Spectra` constructor function
#' which supports the following formats:
#'
#' - parameter `object` is a `DataFrame` containing the spectrum data. The
#'   provided `backend` (by default a [MsBackendDataFrame-class]) will be
#'   initialized with that data.
#'
#' - parameter `object` is a [MsBackend-class] (assumed to be already
#'   initialized).
#'
#' - parameter `object` is missing, in which case it is supposed that the data
#'   is provided by the [MsBackend-class] class passed along with the `backend`
#'   argument.
#'
#' `Spectra` classes are usually created with the `readSpectra`
#' function that reads general spectrum metadata information from the mass
#' spectrometry data files.
#'
#' The backend of a `Spectra` object can be changed with the `setBackend`
#' method that takes an instance of the new backend as second parameter
#' `backend`. A call to `setBackend(sps, backend = MsBackendDataFrame())` would
#' for example change the backend to the *in-memory* `MsBackendDataFrame`.
#' Note that it might not be possible to change from any backend to any other
#' backend. Changing from a `MsBackendDataFrame` to a `MsBackendMzR` would only
#' be possible if the original `MsBackendDataFrame` was generated on data from
#' e.g. mzML files and if these original file names are set in the backend
#' (slot `@files`). In contrast, it should be possible to convert almost every
#' backend into a `MsBackendDataFrame` (given sufficient memory is available).
#'
#' The definition of the function is:
#' `setBackend(object, backend, ..., f = fromFile(object), BPPARAM = bpparam())`
#' and its parameters are:
#'
#' - parameter `object`: the `Spectra` object.
#'
#' - parameter `backend`: an instance of the new backend, e.g.
#'   `MsBackendDataFrame()`.
#'
#' - parameter `f`: factor allowing to parallelize the change of the backends.
#'   By default the process of copying the spectra data from the original to the
#'   new backend is performed separately (and in parallel) for each file.
#'
#' - parameter `...`: optional additional arguments passed to the
#'   [backendInitialize()] method of the new `backend`.
#'
#' - parameter `BPPARAM`: setup for the parallel processing. See [bpparam()] for
#'   details.
#'
#' @section Accessing spectra data:
#'
#' - `$`, `$<-`: get (or set) a spectra variable for all spectra in `object`.
#'   See examples for details.
#'
#' - `acquisitionNum`: returns the acquisition number of each
#'   spectrum. Returns an `integer` of length equal to the number of
#'   spectra (with `NA_integer_` if not available).
#'
#' - `centroided`, `centroided<-`: gets or sets the centroiding
#'   information of the spectra. `centroided` returns a `logical`
#'   vector of length equal to the number of spectra with `TRUE` if a
#'   spectrum is centroided, `FALSE` if it is in profile mode and `NA`
#'   if it is undefined. See also `isCentroided` for estimating from
#'   the spectrum data whether the spectrum is centroided.  `value`
#'   for `centroided<-` is either a single `logical` or a `logical` of
#'   length equal to the number of spectra in `object`.
#'
#' - `collisionEnergy`, `collisionEnergy<-`: gets or sets the
#'   collision energy for all spectra in `object`. `collisionEnergy`
#'   returns a `numeric` with length equal to the number of spectra
#'   (`NA_real_` if not present/defined), `collisionEnergy<-` takes a
#'   `numeric` of length equal to the number of spectra in `object`.
#'
#' - `fileNames`: returns a `character` with the file names, or
#'   `NA_character_` if not relevant.
#'
#' - `fromFile`: get the file/sample assignment of each spectrum. Returns an
#'   integer vector of length equal to the number of spectra.
#'
#' - `intensity`: gets the intensity values from the spectra. Returns
#'   a [NumericList()] of `numeric` vectors (intensity values for each
#'   spectrum). The length of the `list` is equal to the number of
#'   `spectra` in `object`.
#'
#' - `ionCount`: returns a `numeric` with the sum of intensities for
#'   each spectrum. If the spectrum is empty (see `isEmpty`),
#'   `NA_real_` is returned.
#'
#' - `isCentroided`: a heuristic approach assessing if the spectra in
#'   `object` are in profile or centroided mode. The function takes
#'   the `qtl`th quantile top peaks, then calculates the difference
#'   between adjacent m/z value and returns `TRUE` if the first
#'   quartile is greater than `k`. (See `Spectra:::.isCentroided` for
#'   the code.)
#'
#' - `isEmpty`: checks whether a spectrum in `object` is empty
#'   (i.e. does not contain any peaks). Returns a `logical` vector of
#'   length equal number of spectra.
#'
#' - `isolationWindowLowerMz`, `isolationWindowLowerMz<-`: get or set the lower
#'   m/z boundary of the isolation window.
#'
#' - `isolationWindowTargetMz`, `isolationWindowTargetMz<-`: get or set the
#'   target m/z of the isolation window.
#'
#' - `isolationWindowUpperMz`, `isolationWindowUpperMz<-`: get or set the upper
#'   m/z boundary of the isolation window.
#'
#' - `length`: get the number of spectra in the object.
#'
#' - `msLevel`: get the spectra's MS level. Returns an integer vector (names
#'   being spectrum names, length equal to the number of spectra) with the MS
#'   level for each spectrum.
#'
#' - `mz`: gets the mass-to-charge ratios (m/z) from the
#'   spectra. Returns a [NumericList()] or length equal to the number of
#'   spectra, each element a `numeric` vector with the m/z values of
#'   one spectrum.
#'
#' - `peaks`: get the *peaks* matrices for all spectra in `object`. The function
#'   returns a [SimpleList()] of matrices, each `matrix` with columns `mz` and
#'   `intensity` with the m/z and intensity values for all peaks of a spectrum.
#'
#' - `peaksCount`: gets the number of peaks (m/z-intensity values) per
#'   spectrum. Returns an `integer` vector (length equal to the
#'   number of spectra). For empty spectra, `NA_integer_` is returned.
#'
#' - `polarity`, `polarity<-`: gets or sets the polarity for each
#'   spectrum.  `polarity` returns an `integer` vector (length equal
#'   to the number of spectra), with `0` and `1` representing negative
#'   and positive polarities, respectively. `polarity<-` expects an
#'   `integer` vector of length 1 or equal to the number of spectra.
#'
#' - `precursorCharge`, `precursorIntensity`, `precursorMz`,
#'   `precScanNum`, `precAcquisitionNum`: get the charge (`integer`),
#'   intensity (`numeric`), m/z (`numeric`), scan index (`integer`)
#'   and acquisition number (`interger`) of the precursor for MS level
#'   > 2 spectra from the object. Returns a vector of length equal to
#'   the number of spectra in `object`. `NA` are reported for MS1
#'   spectra of if no precursor information is available.
#'
#' - `rtime`, `rtime<-`: gets or sets the retention times for each
#'   spectrum.  `rtime` returns a `numeric` vector (length equal to
#'   the number of spectra) with the retention time for each spectrum.
#'   `rtime<-` expects a numeric vector with length equal to the
#'   number of spectra.
#'
#' - `scanIndex`: returns an `integer` vector with the *scan index*
#'   for each spectrum. This represents the relative index of the
#'   spectrum within each file. Note that this can be different to the
#'   `acquisitionNum` of the spectrum which represents the index of the
#'   spectrum during acquisition/measurement (as reported in the mzML file).
#'
#' - `smoothed`,`smoothed<-`: gets or sets whether a spectrum is
#'   *smoothed*. `smoothed` returns a `logical` vector of length equal
#'   to the number of spectra. `smoothed<-` takes a `logical` vector
#'   of length 1 or equal to the number of spectra in `object`.
#'
#' - `spectraData`, `spectraData<-`: get or sets general spectrum
#'   metadata (annotation, also called header). `spectraData` returns
#'   a `DataFrame`, `spectraData<-` expects a `DataFrame`. Note that not all
#'   backends support replacing all spectra variables (the [MsBackendMzR()]
#'   does for example not allow to replace `mz` and `intensity` values with the
#'   `spectraData<-` method.
#'
#' - `spectraNames`, `spectraNames<-`: get or set the spectra names.
#'
#' - `spectraVariables`: returns a `character` vector with the
#'   available spectra variables (columns, fields or attributes)
#'   available in `object`.
#'
#' - `tic`: gets the total ion current/count (sum of signal of a
#'   spectrum) for all spectra in `object`. By default, the value
#'   reported in the original raw data file is returned. For an empty
#'   spectrum, `0` is returned.
#'
#' @section Data subsetting, filtering and merging:
#'
#' Subsetting and filtering of `Spectra` objects can be performed with the below
#' listed methods.
#'
#' - `[`: subset the spectra keeping only selected elements (`i`). The method
#'   **always** returns a `Spectra` object.
#'
#' - `filterAcquisitionNum`: filter the object keeping only spectra matching the
#'   provided acquisition numbers (argument `n`). If `file` is also provided,
#'   `object` is subsetted to the spectra with an acquisition number equal to
#'   `n` **in this/these file(s)** and all spectra for the remaining files (not
#'   specified with `file`). Returns the filtered `Spectra`.
#'
#' - `filterEmptySpectra`: remove empty spectra (i.e. spectra without peaks).
#'
#' - `filterFile`: retain data of files matching the file index or file name
#'    provided with parameter `file`. Returns the filtered `Spectra`.
#'
#' - `filterIsolationWindow`: retain spectra that contain `mz` in their
#'   isolation window m/z range (i.e. with an `isolationWindowLowerMz` <= `mz`
#'   and `isolationWindowUpperMz` >= `mz`.
#'
#' - `filterMsLevel`: filter object by MS level keeping only spectra matching
#'   the MS level specified with argument `msLevel`. Returns the filtered
#'   `Spectra`.
#'
#' - `filterPolarity`: filter the object keeping only spectra matching the
#'   provided polarity. Returns the subsetted `Spectra`.
#'
#' - `filterPrecursorMz`: retain spectra with an m/z matching the provided `mz`
#'   accepting also a small difference in m/z which can be defined by parameter
#'   `ppm` (parts per million). With the default (`ppm = 0`) only spectra with
#'   m/z identical to `mz` are retained.
#'
#' - `filterPrecursorScan`: retain parent (e.g. MS1) and children scans (e.g.
#'    MS2) of acquisition number `acquisitionNum`. Returns the filtered
#'    `Spectra`.
#'
#' - `filterRt`: retain spectra of MS level `msLevel` with retention times
#'    within (`>=`) `rt[1]` and (`<=`) `rt[2]`.
#'
#' - `selectSpectraVariables`: reduce the information within the object to
#'   the selected spectra variables: all data for variables not specified will
#'   be dropped. For mandatory columns (such as *msLevel*, *rtime* ...) only
#'   the values will be dropped, while additional (user defined) spectra
#'   variables will be completely removed. Returns the filtered `Spectra`.
#'
#' Several `Spectra` objects can be merged into a single object with the
#' `merge` function. Merging will fail if the processing queue of any of the
#' `Spectra` objects is not empty, or if different backends are used in the
#' various `Spectra` objects. The spectra variables of the resulting `Spectra`
#' object is the union of the spectra variables of the individual `Spectra`
#' objects.
#'
#' @section Data manipulation and analysis methods:
#'
#' Many data manipulation operations, such as those listed in this section, are
#' not applied immediately to the spectra, but added to a
#' *lazy processinq queue*. Operations stored in this queue are applied
#' on-the-fly to spectra data each time it is accessed. This lazy
#' execution guarantees the same functionality for `Spectra` objects with
#' any backend, i.e. backends supporting to save changes to spectrum data
#' ([MsBackendDataFrame()] as well as read-only backends (such
#' as the [MsBackendMzR()]).
#'
#' - `addProcessing`: add an arbitrary function that should be applied to the
#'   peaks matrix of every spectrum in `object`. The function (can be passed
#'   with parameter `FUN`) is expected to take a peaks matrix as input and to
#'   return a peaks matrix. A peaks matrix is a numeric matrix with two columns,
#'   the first containing the m/z values of the peaks and the second the
#'   corresponding intensities. The function has to have `...` in its
#'   definition. Additional arguments can be passed with `...`. Examples are
#'   provided in the package vignette.
#'
#' - `clean`: remove 0-intensity data points. For `all = FALSE` (the default)
#'   0-intensity peaks next to non-zero intensity peaks are retained while with
#'   `all = TRUE` all 0-intensity peaks are removed.
#'
#' - `removePeaks`: *remove* peaks lower or equal to a threshold intensity
#'   value `t` by setting their intensity to `0`. With the default `t = "min"`
#'   all peaks with an intensity smaller or equal to the minimal non-zero
#'   intensity is set to `0`. If the spectrum is in profile mode, ranges of
#'   successive non-0 peaks <= `t` are set to 0. If the spectrum is centroided,
#'   then individual peaks <= `t` are set to 0. Note that the number of peaks
#'   is not changed unless `clean` is called after `removePeaks`.
#'
#' @return See individual method description for the return value.
#'
#' @param all for `clean`: `logical(1)` whether all 0 intensity peaks should be
#'     removed (`TRUE`) or whether 0-intensity peaks directly adjacent to a
#'     non-zero intensity peak should be kept (`FALSE`).
#'
#' @param acquisitionNum for `filterPrecursorScan`: `integer` with the
#'     acquisition number of the spectra to which the object should be
#'     subsetted.
#'
#' @param backend For `Spectra`: [MsBackend-class] to be used as backend. See
#'     section on creation of `Spectra` objects for details. For `setBackend`:
#'     instance of [MsBackend-class]. See section on creation of `Spectra`
#'     objects for details.
#'
#' @param BPPARAM Parallel setup configuration. See [bpparam()] for more
#'     information. This is passed directly to the [backendInitialize()] method
#'     of the [MsBackend-class].
#'
#' @param columns For `spectraData` accessor: optional `character` with column
#'     names (spectra variables) that should be included in the
#'     returned `DataFrame`. By default, all columns are returned.
#'
#' @param drop For `[`: not considered.
#'
#' @param f For `setBackend`: factor defining how to split the data should
#'     for the parallelized copying of the spectra data to the new backend.
#'
#' @param file For `filterFile`: index or name of the file(s) to which the data
#'     should be subsetted.
#'
#' @param FUN For `addProcessing`: function to be applied to the peak matrix
#'     of each spectrum in `object`. See section *Data manipulations* below
#'     for more details.
#'
#' @param i For `[`: `integer`, `logical` or `character` to subset the object.
#'
#' @param j For `[`: not supported.
#'
#' @param initial For `tic`: `logical(1)` whether the initially
#'     reported total ion current should be reported, or whether the
#'     total ion current should be (re)calculated on the actual data
#'     (`initial = FALSE`, same as `ionCount`).
#'
#' @param n for `filterAcquisitionNum`: `integer` with the acquisition numbers
#'     to filter for.
#'
#' @param name For `$` and `$<-`: the name of the spectra variable to return
#'     or set.
#'
#' @param metadata For `Spectra`: optional `list` with metadata information.
#'
#' @param msLevel For `filterMsLevel`: the MS level to which `object` should be
#'     subsetted.
#'
#' @param msLevel. `integer` defining the MS level(s) of the spectra to which
#'     the function should be applied. For `filterMsLevel`: the MS level to
#'     which `object` should be subsetted.
#'
#' @param mz For `filterIsolationWindow` and `filterPrecursorMz`: `numeric(1)`
#'     with the m/z value to filter the object.
#'
#' @param object For `Spectra`: either a `DataFrame` or `missing`. See section
#'     on creation of `Spectra` objects for details. For all other methods a
#'     `Spectra` object.
#'
#' @param polarity for `filterPolarity`: `integer` specifying the polarity to
#'     to subset `object`.
#'
#' @param ppm For `filterPrecursorMz`: `numeric(1)` defining the accepted
#'     difference between the provided m/z and the spectrum's m/z in parts per
#'     million.
#'
#' @param processingQueue For `Spectra`: optional `list` of
#'     [ProcessingStep-class] objects.
#'
#' @param spectraVariables For `selectSpectraVariables`: `character` with the
#'     names of the spectra variables to which the backend should be subsetted.
#'
#' @param rt for `filterRt`: `numeric(2)` defining the retention time range to
#'     be used to subset/filter `object`.
#'
#' @param t for `removePeaks`: a `numeric(1)` defining the threshold or `"min"`.
#'
#' @param x A `Spectra` object.
#'
#' @param y A `Spectra` object.
#'
#' @param value replacement value for `<-` methods. See individual
#'     method description or expected data type.
#'
#' @param ... Additional arguments.
#'
#' @author Sebastian Gibb, Johannes Rainer
#'
#' @md
#'
#' @exportClass Spectra
#'
#' @exportMethod Spectra
#'
#' @examples
#'
#' ## Create a Spectra providing a `DataFrame` containing the spectrum data.
#'
#' spd <- DataFrame(msLevel = c(1L, 2L), rtime = c(1.1, 1.2))
#' spd$mz <- list(c(100, 103.2, 104.3, 106.5), c(45.6, 120.4, 190.2))
#' spd$intensity <- list(c(200, 400, 34.2, 17), c(12.3, 15.2, 6.8))
#'
#' data <- Spectra(spd)
#' data
#'
#' ## Create a Spectra from a mzML file.
#' sciex_file <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' sciex_mzr <- backendInitialize(MsBackendMzR(), files = sciex_file)
#' sciex <- Spectra(sciex_mzr)
#' sciex
#'
#' ## The MS data is on-disk and will be read into memory on-demand. We can
#' ## however change the backend to a MsBackendDataFrame backend which will
#' ## keep all of the data in memory.
#' sciex_im <- setBackend(sciex, MsBackendDataFrame())
#' sciex_im
#'
#' ## The size of the objects will obviously be different
#' object.size(sciex)
#' object.size(sciex_im)
#'
#' ## ---- ACCESSING AND ADDING DATA ----
#'
#' ## Get the MS level for each spectrum.
#' msLevel(data)
#'
#' ## Alternatively, we could also use $ to access a specific spectra variable.
#' ## This could also be used to add additional spectra variables to the
#' ## object (see further below).
#' data$msLevel
#'
#' ## Get the intensity and m/z values.
#' intensity(data)
#' mz(data)
#'
#' ## Get the m/z for the first spectrum.
#' mz(data)[[1]]
#'
#' ## Get the peak data (m/z and intensity values).
#' pks <- peaks(data)
#' pks
#' pks[[1]]
#' pks[[2]]
#'
#' ## List all available spectra variables (i.e. spectrum data and metadata).
#' spectraVariables(data)
#'
#' ## For all *core* spectrum variables accessor functions are available. These
#' ## return NA if the variable was not set.
#' centroided(data)
#' fromFile(data)
#' rtime(data)
#' precursorMz(data)
#'
#' ## Add an additional metadata column.
#' data$spectrum_id <- c("sp_1", "sp_2")
#'
#' ## List spectra variables, "spectrum_id" is now also listed
#' spectraVariables(data)
#'
#' ## Get the values for the new spectra variable
#' data$spectrum_id
#'
#' ## Extract specific spectra variables.
#' spectraData(data, columns = c("spectrum_id", "msLevel"))
#'
#' ## Drop spectra variable data and/or columns.
#' res <- selectSpectraVariables(data, c("mz", "intensity"))
#'
#' ## This removed the additional columns "spectrum_id" and deleted all values
#' ## for all spectra variables, except "mz" and "intensity".
#' spectraData(res)
#'
#' ## Compared to the data before selectSpectraVariables.
#' spectraData(data)
#'
#'
#' ## ---- SUBSETTING, FILTERING AND COMBINING
#'
#' ## Subset to all MS2 spectra.
#' data[msLevel(data) == 2]
#'
#' ## Same with the filterMsLevel function
#' filterMsLevel(data, 2)
#'
#' ## Below we combine the `data` and `sciex_im` objects into a single one.
#' data_comb <- merge(data, sciex_im)
#'
#' ## The combined Spectra contains a union of all spectra variables:
#' head(data_comb$spectrum_id)
#' head(data_comb$rtime)
#' head(data_comb$fromFile)
#'
#' ## ---- DATA MANIPULATIONS ----
#'
#' ## Set the data to be centroided
#' centroided(data) <- TRUE
#'
#' ## Remove peaks with an intensity below 40.
#' res <- removePeaks(data, t = 40)
#' res
#'
#' ## Get the intensities of the first and second spectrum.
#' intensity(res)[[1]]
#' intensity(res)[[2]]
#'
#' ## Clean all spectra removing all 0-intensity peaks.
#' res <- clean(res, all = TRUE)
#'
#' ## Get the intensities of the first and second spectrum.
#' intensity(res)[[1]]
#' intensity(res)[[2]]
#'
#' ## Second spectrum is now empty:
#' isEmpty(res)
NULL

#' The Spectra class
#'
#' The [Spectra-class] encapsulates data and meta-data for mass
#' spectrometry experiments.
#'
#'
#' @slot backend A derivate of [MsBackend-class] holding/controlling the spectra
#' data.
#' @slot processingQueue `list` of `ProcessingStep` objects.
#' @slot processing A `character` storing logging information.
#' @slot metadata A `list` storing experiment metadata.
#' @slot version A `characher(1)` containing the class version.
#'
#' @name Spectra-class
#' @docType class
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#'
#' @importClassesFrom S4Vectors DataFrame
#'
#' @importMethodsFrom S4Vectors lapply
#'
#' @importFrom S4Vectors DataFrame
#'
#' @noRd
setClass(
    "Spectra",
    slots = c(
        backend = "MsBackend",
        processingQueue = "list",
        ## logging
        processing = "character",
        ## metadata
        metadata = "list",
        version = "character"
    ),
    prototype = prototype(version = "0.1")
)

setValidity("Spectra", function(object) {
    msg <- .valid_processing_queue(object@processingQueue)
    if (length(msg)) msg
    else TRUE
})

#' @rdname hidden_aliases
#'
#' @importMethodsFrom methods show
#'
#' @importFrom utils capture.output
#'
#' @exportMethod show
setMethod("show", "Spectra",
    function(object) {
        cat("MSn data (", class(object)[1L], ") with ",
            length(object@backend), " spectra in a ", class(object@backend),
            " backend:\n", sep = "")
        if (length(object@backend)) {
            txt <- capture.output(show(object@backend))
            cat(txt[-1], sep = "\n")
        }
        if (length(object@processingQueue))
            cat("Lazy evaluation queue:", length(object@processingQueue),
                "processing step(s)\n")
        cat("Processing:\n", paste(object@processing, collapse="\n "), "\n")
    })

#' @rdname Spectra
setMethod("Spectra", "DataFrame", function(object, processingQueue = list(),
                                           metadata = list(), ...,
                                           backend = MsBackendDataFrame(),
                                           BPPARAM = bpparam()) {
    object$fromFile <- rep(1L, nrow(object))
    new("Spectra", metadata = metadata, processingQueue = processingQueue,
        backend = backendInitialize(
            backend, files = if (nrow(object)) NA_character_ else character(),
            object, BPPARAM = BPPARAM))
})

#' @rdname Spectra
setMethod("Spectra", "missing", function(object, processingQueue = list(),
                                         metadata = list(), ...,
                                         backend = MsBackendDataFrame(),
                                         BPPARAM = bpparam()) {
    new("Spectra", metadata = metadata, processingQueue = processingQueue,
        backend = backend)
})

#' @rdname Spectra
setMethod("Spectra", "MsBackend", function(object, processingQueue = list(),
                                           metadata = list(), ...,
                                           BPPARAM = bpparam()) {
    new("Spectra", metadata = metadata, processingQueue = processingQueue,
        backend = object)
})

#' @rdname Spectra
#'
#' @exportMethod setBackend
setMethod("setBackend", c("Spectra", "MsBackend"),
          function(object, backend, f = fromFile(object), ...,
                   BPPARAM = bpparam()) {
              backend_class <- class(object@backend)
              f <- factor(f, levels = unique(f))
              if (length(f) != length(object))
                  stop("length of 'f' has to match the length of 'object'")
              bknds <- bplapply(split(object@backend, f = f), function(z, ...) {
                  if (isReadOnly(backend) && any(z@modCount))
                      stop(class(backend), " backends are read-only but the ",
                           "original data appears to be modified")
                  res <- backendInitialize(backend, files = z@files,
                                           spectraData = spectraData(z),
                                           ...)
                  res@modCount <- z@modCount
                  res
              }, ..., BPPARAM = BPPARAM)
              bknds <- backendMerge(bknds)
              if (is.unsorted(f))
                  bknds <- bknds[order(unlist(split(seq_along(bknds), f),
                                              use.names = FALSE))]
              object@backend <- bknds
              object@processing <- .logging(object@processing,
                                            "Switch backend from ",
                                            backend_class, " to ",
                                            class(object@backend))
              object
          })

#' @rdname Spectra
#'
#' @importMethodsFrom S4Vectors merge
#'
#' @exportMethod merge
setMethod("merge", c("Spectra", "Spectra"), function(x, y, ...) {
    objs <- unname(c(x, y, ...))
    pqs <- lapply(objs, function(z) z@processingQueue)
    ## For now we stop if there is any of the processingQueues not empty. Later
    ## we could even test if they are similar, and if so, merge.
    if (any(lengths(pqs)))
        stop("Can not merge Spectra objects with non-empty processing queue")
    metad <- do.call(c, lapply(objs, function(z) z@metadata))
    procs <- unique(unlist(lapply(objs, function(z) z@processing)))
    object <- new(
        "Spectra", metadata = metad,
        backend = backendMerge(lapply(objs, function(z) z@backend)),
        processing = c(procs, paste0("Merge ", length(objs),
                                     " Spectra into one [", date(), "]"))
    )
    validObject(object)
    object
})

#### ---------------------------------------------------------------------------
##
##                          ACCESSOR METHODS
##
#### ---------------------------------------------------------------------------

#' @rdname Spectra
setMethod("acquisitionNum", "Spectra", function(object)
    acquisitionNum(object@backend))

#' @rdname Spectra
setMethod("centroided", "Spectra", function(object) {
    centroided(object@backend)
})

#' @rdname Spectra
setReplaceMethod("centroided", "Spectra", function(object, value) {
    centroided(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("collisionEnergy", "Spectra", function(object) {
    collisionEnergy(object@backend)
})

#' @rdname Spectra
setReplaceMethod("collisionEnergy", "Spectra", function(object, value) {
    collisionEnergy(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("fileNames", "Spectra", function(object) {
    fileNames(object@backend)
})

#' @rdname Spectra
setMethod("fromFile", "Spectra", function(object) fromFile(object@backend))

#' @rdname Spectra
setMethod("intensity", "Spectra", function(object, ...) {
    NumericList(lapply(.peaksapply(object, ...), function(z) z[, 2]),
                compress = FALSE)
})

#' @rdname Spectra
setMethod("ionCount", "Spectra", function(object) {
    if (length(object))
        unlist(.peaksapply(object, FUN = function(pks, ...)
            sum(pks[, 2], na.rm = TRUE)), use.names = FALSE)
    else numeric()
})

#' @rdname Spectra
setMethod("isCentroided", "Spectra", function(object, ...) {
    if (length(object))
        unlist(.peaksapply(object, FUN = .is_centroided_peaks),
               use.names = FALSE)
    else logical()
})

#' @rdname Spectra
setMethod("isEmpty", "Spectra", function(x) {
    if (length(x))
        unlist(.peaksapply(x, FUN = function(pks, ...) nrow(pks) == 0),
               use.names = FALSE)
    else logical()
})

#' @rdname Spectra
setMethod("isolationWindowLowerMz", "Spectra", function(object) {
    isolationWindowLowerMz(object@backend)
})

#' @rdname Spectra
setReplaceMethod("isolationWindowLowerMz", "Spectra", function(object, value) {
    isolationWindowLowerMz(object) <- value
    object
})

#' @rdname Spectra
setMethod("isolationWindowTargetMz", "Spectra", function(object) {
    isolationWindowTargetMz(object@backend)
})

#' @rdname Spectra
setReplaceMethod("isolationWindowTargetMz", "Spectra", function(object, value) {
    isolationWindowTargetMz(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("isolationWindowUpperMz", "Spectra", function(object) {
    isolationWindowUpperMz(object@backend)
})

#' @rdname Spectra
setReplaceMethod("isolationWindowUpperMz", "Spectra", function(object, value) {
    isolationWindowUpperMz(object@backend) <- value
    object
})

#' @rdname Spectra
#'
#' @exportMethod length
setMethod("length", "Spectra", function(x) length(x@backend))

#' @rdname Spectra
setMethod("msLevel", "Spectra", function(object) msLevel(object@backend))

#' @rdname Spectra
setMethod("mz", "Spectra", function(object, ...) {
    NumericList(lapply(.peaksapply(object, ...), function(z) z[, 1]),
                compress = FALSE)
})

#' @rdname Spectra
#'
#' @exportMethod peaks
setMethod("peaks", "Spectra", function(object, ...) {
    SimpleList(.peaksapply(object, ...))
})

#' @rdname Spectra
setMethod("peaksCount", "Spectra", function(object) {
    if (length(object))
        unlist(.peaksapply(object, FUN = function(pks, ...) nrow(pks)),
               use.names = FALSE)
    else integer()
})

#' @rdname Spectra
setMethod("polarity", "Spectra", function(object) {
    polarity(object@backend)
})

#' @rdname Spectra
setReplaceMethod("polarity", "Spectra", function(object, value) {
    polarity(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("precScanNum", "Spectra", function(object) {
    precScanNum(object@backend)
})

#' @rdname Spectra
setMethod("precursorCharge", "Spectra", function(object) {
    precursorCharge(object@backend)
})

#' @rdname Spectra
setMethod("precursorIntensity", "Spectra", function(object) {
    precursorIntensity(object@backend)
})

#' @rdname Spectra
setMethod("precursorMz", "Spectra", function(object) {
    precursorMz(object@backend)
})

#' @rdname Spectra
setMethod("rtime", "Spectra", function(object) {
    rtime(object@backend)
})

#' @rdname Spectra
setReplaceMethod("rtime", "Spectra", function(object, value) {
    rtime(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("scanIndex", "Spectra", function(object) {
    scanIndex(object@backend)
})

#' @rdname Spectra
setMethod("selectSpectraVariables", "Spectra",
          function(object, spectraVariables = spectraVariables(object)) {
              spectraVariables <- union(spectraVariables, "fromFile")
              object@backend <- selectSpectraVariables(
                  object@backend, spectraVariables = spectraVariables)
              object
})

#' @rdname Spectra
setMethod("smoothed", "Spectra", function(object) {
    smoothed(object@backend)
})

#' @rdname Spectra
setReplaceMethod("smoothed", "Spectra", function(object, value) {
    smoothed(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("spectraData", "Spectra", function(object,
                                             columns = spectraVariables(object))
{
    skip_cols <- c("mz", "intensity") # Eventually also tic, peak count etc
    if (all(columns %in% skip_cols))
        res <- DataFrame(matrix(ncol = 0, nrow = length(object)))
    else res <- spectraData(object@backend,
                            columns = columns[!(columns %in% skip_cols)])
    if (any(columns %in% skip_cols)) {
        pks <- .peaksapply(object)
        if (any(columns == "mz"))
            res$mz <- NumericList(lapply(pks, function(z) z[, 1]),
                                  compress = FALSE)
        if (any(columns == "intensity"))
            res$intensity <- NumericList(lapply(pks, function(z) z[, 2]),
                                         compress = FALSE)
    }
    res[, columns, drop = FALSE]
})

#' @rdname Spectra
setReplaceMethod("spectraData", "Spectra", function(object, value) {
    spectraData(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("spectraNames", "Spectra", function(object) {
    spectraNames(object@backend)
})

#' @rdname Spectra
setReplaceMethod("spectraNames", "Spectra", function(object, value) {
    spectraNames(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("spectraVariables", "Spectra", function(object) {
    spectraVariables(object@backend)
})

#' @rdname Spectra
setMethod("tic", "Spectra", function(object, initial = TRUE) {
    if (!length(object))
        return(numeric())
    if (initial)
        tic(object@backend, initial = initial)
    else ionCount(object)
})

#' @rdname Spectra
#'
#' @importMethodsFrom S4Vectors $
setMethod("$", "Spectra", function(x, name) {
    if (!(name %in% spectraVariables(x)))
        stop("No spectra variable '", name, "' available")
    ## Use spectraData instead of x@backend$name to support the processing
    ## queue.
    spectraData(x, column = name)[, 1]
})

#' @rdname Spectra
setReplaceMethod("$", "Spectra", function(x, name, value) {
    ## if (any(name %in% c("mz", "intensity")))
    ##     stop("Replacing mz or intensity values is currently not supported")
    x@backend <- do.call("$<-", list(x@backend, name, value))
    x
})

#### ---------------------------------------------------------------------------
##
##                      FILTERING AND SUBSETTING
##
#### ---------------------------------------------------------------------------

#' @rdname Spectra
setMethod("[", "Spectra", function(x, i, j, ..., drop = FALSE) {
    if (!missing(j))
        stop("Subsetting 'Spectra' by columns is not (yet) supported")
    if (missing(i))
        i <- seq_len(length(x))
    x@backend <- x@backend[i = i, ..., drop = drop]
    x
})

#' @rdname Spectra
setMethod("filterAcquisitionNum", "Spectra", function(object, n = integer(),
                                                      file = integer()) {
    object@backend <- filterAcquisitionNum(object@backend, n, file)
    object@processing <- .logging(object@processing,
                                  "Filter: select by: ", length(n),
                                  " acquisition number(s) in ", length(file),
                                  " file(s)")
    object
})

#' @rdname Spectra
setMethod("filterEmptySpectra", "Spectra", function(object) {
    object@backend <- object@backend[as.logical(peaksCount(object))]
    object@processing <- .logging(object@processing,
                                  "Filter: removed empty spectra.")
    object
})

#' @rdname Spectra
setMethod("filterFile", "Spectra", function(object, file = integer()) {
    object@backend <- filterFile(object@backend, file = file)
    object@processing <- .logging(object@processing,
                                  "Filter: select file(s) ",
                                  paste0(file, collapse = ", "))
    object
})

#' @rdname Spectra
setMethod("filterIsolationWindow", "Spectra", function(object, mz = numeric()) {
    object@backend <- filterIsolationWindow(object@backend, mz = mz)
    object@processing <- .logging(object@processing,
                                  "Filter: select spectra containing m/z ",
                                  mz, " in their isolation window")
    object
})

#' @rdname Spectra
setMethod("filterMsLevel", "Spectra", function(object, msLevel. = integer()) {
    object@backend <- filterMsLevel(object@backend, msLevel = msLevel.)
    object@processing <- .logging(object@processing,
                                  "Filter: select MS level(s) ",
                                  paste0(unique(msLevel.), collapse = " "))
    object
})

#' @rdname Spectra
setMethod("filterPolarity", "Spectra", function(object, polarity = integer()) {
    object@backend <- filterPolarity(object@backend, polarity = polarity)
    object@processing <- .logging(object@processing,
                                  "Filter: select spectra with polarity ",
                                  paste0(polarity, collapse = " "))
    object
})

#' @rdname Spectra
setMethod("filterPrecursorMz", "Spectra",
          function(object, mz = numeric(), ppm = 0) {
              object@backend <- filterPrecursorMz(object@backend, mz, ppm)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select spectra with a precursor m/z of ",
                  mz, " accepting ", ppm, " ppm difference")
              object
          })

#' @rdname Spectra
setMethod("filterPrecursorScan", "Spectra",
          function(object, acquisitionNum= integer()) {
              object@backend <- filterPrecursorScan(object@backend,
                                                    acquisitionNum)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select parent/children scans for ",
                  paste0(acquisitionNum, collapse = " "))
              object
          })

#' @rdname Spectra
setMethod("filterRt", "Spectra",
          function(object, rt = numeric(), msLevel. = unique(msLevel(object))) {
              suppressWarnings(rt <- range(rt))
              object@backend <- filterRt(object@backend, rt, msLevel.)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select retention time [", rt[1], "..", rt[2],
                  "] on MS level(s) ", paste0(msLevel., collapse = " "))
              object
          })

#### ---------------------------------------------------------------------------
##
##                      DATA MANIPULATION METHODS
##
#### ---------------------------------------------------------------------------

#' @rdname Spectra
#'
#' @exportMethod removePeaks
setMethod("removePeaks", "Spectra",
          function(object, t = "min", msLevel. = unique(msLevel(object))) {
              if (!is.numeric(t) & t != "min")
                  stop("Argument 't' has to be either numeric of 'min'.")
              if (!is.numeric(msLevel.))
                  stop("'msLevel.' must be numeric.")
              object <- addProcessing(object, .remove_peaks, t = t,
                                      msLevel = msLevel.)
              object@processing <- .logging(object@processing,
                                            "Signal <= ", t, " in MS level(s) ",
                                            paste0(msLevel., collapse = ", "),
                                            " set to 0")
              object
          })

#' @rdname Spectra
#'
#' @exportMethod clean
setMethod("clean", "Spectra",
          function(object, all = FALSE, msLevel. = unique(msLevel(object))) {
              if (!is.logical(all) || length(all) != 1)
                  stop("Argument 'all' must be a logical of length 1")
              if (!is.numeric(msLevel.))
                  stop("'msLevel' must be numeric.")
              object <- addProcessing(object, .clean_peaks, all = all,
                                      msLevel = msLevel.)
              object@processing <- .logging(object@processing,
                                            "Spectra of MS level(s) ",
                                            paste0(msLevel., collapse = ", "),
                                            " cleaned ")
              object
          })

applyProcessing <- function(object, f = fromFile(object), BPPARAM = bpparam(),
                            ...) {
    if (!length(object@processingQueue))
        return(object)
}

## applyProcessing:
## bknds <- bplapply(split(object@backend, f = f), function(z, ...) {
##     if (isReadOnly(z))
##         stop("Can not replace peaks data because ", class(z), " backends ",
##         "are read-only")
##     peaks(z) <- .apply_processing_queue(peaks(z), msLevel(z), centroided(z), queue)
##     z@modCount <- 0
##     z
## }, ..., BPPARAM = BPPARAM)
## bknds <- backendMerge(bknds)
## if (is.unsorted(f))
##     bknds <- bknds[order(unlist(split(seq_along(bknds), f),
##                                 use.names = FALSE))]
## object@backend <- bknds
