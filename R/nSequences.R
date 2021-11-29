#' @title Count the number of identifications per scan
#'
#' @description
#'
#' The function takes a proteomics `Spectra` object containing
#' identification results as input. It then counts the number of
#' identifications each scan (or their descendants) has lead to - this
#' is either 0 or 1 for MS2 scans, or, for MS1 scans, the number of
#' MS2 scans originating from any MS1 peak that lead to an
#' identification.
#'
#' This function can be used to generate id-annotated total ion
#' chromatograms, as can illustrated
#' [here](https://rformassspectrometry.github.io/docs/sec-id.html#an-identification-annotated-chromatogram).
#'
#' @details
#'
#' The computed number of identifications is stored in a new spectra
#' variables named `"nSequences"`. If it already exists, the function
#' throws a message and returns the object unchanged. To force the
#' recomputation of the `"nSequences"` variables, users should either
#' delete or rename name.
#'
#' @param object An instance of class `Spectra()` that contains
#'     identification data, as defined by the `sequence` argument.
#'
#' @param sequence `character(1)` with the name of the spectra
#'     variable that defines wheather a scan lead to an
#'     identification, typically containing the idenfified peptides
#'     sequence. The absence of sequence or identification is encode
#'     by an `NA`. Default is `"sequence"`.
#'
#' @param f A `factor` defining how to split `object` for parallelized
#'     processing. Default is `dataOrigin(x)`, i.e. each raw data
#'     files is processed in parallel.
#'
#' @param BPPARAM Parallel setup configuration. See
#'     [BiocParallel::bpparam()] for details.
#'
#' @return An updated [Spectra()] object that now contains an integer
#'     spectra variable `nSequence` with the number of identification
#'     for each scan.
#'
#' @author Laurent Gatto
#'
#' @export
#'
#' @examples
#' spdf <- new("DFrame", rownames = NULL, nrows = 86L,
#'    listData = list(
#'        msLevel = c(1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
#'                    2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L,
#'                    2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
#'                    2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
#'                    2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
#'                    2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L,
#'                    2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
#'                    2L, 2L),
#'        acquisitionNum = 8975:9060,
#'        precScanNum = c(NA, 8956L, 8956L, 8956L, 8956L, 8956L, 8956L,
#'                        8956L, 8956L, 8956L, 8956L, 8956L, 8956L,
#'                        8956L, 8956L, 8956L, 8956L, 8956L, 8956L, NA,
#'                        8975L, 8975L, 8975L, 8975L, 8975L, 8975L,
#'                        8975L, 8975L, 8975L, 8975L, 8975L, 8975L,
#'                        8975L, 8975L, 8975L, 8975L, 8975L, NA, 8994L,
#'                        8994L, 8994L, 8994L, 8994L, 8994L, 8994L,
#'                        8994L, 8994L, 8994L, 8994L, 8994L, 8994L, NA,
#'                        9012L, 9012L, 9012L, 9012L, 9012L, 9012L,
#'                        9012L, 9012L, 9012L, 9012L, 9012L, 9012L,
#'                        9012L, 9012L, 9012L, 9012L, 9012L, 9012L, NA,
#'                        9026L, 9026L, 9026L, 9026L, 9026L, 9026L,
#'                        9026L, 9026L, 9026L, 9026L, 9026L, 9026L,
#'                        9026L, 9026L, 9026L),
#'        sequence = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
#'                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
#'                     "LSEHATAPTR", NA, NA, NA, NA, NA, NA, NA,
#'                     "EGSDATGDGTK", NA, NA, "NEDEDSPNK", NA, NA, NA,
#'                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
#'                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
#'                     NA, NA, NA, NA, NA, NA, NA, NA, NA, "GLTLAQGGVK",
#'                     NA, NA, NA, NA, "STLPDADRER", NA, NA, NA, NA, NA,
#'                     NA, NA, NA)),
#'    elementType = "ANY", elementMetadata = NULL, metadata = list())
#'
#' sp <- Spectra(spdf)
#'
#' ## We have in this data 5 MS1 and 81 MS2 scans
#' table(msLevel(sp))
#'
#' ## The acquisition number of the MS1 scans
#' acquisitionNum(filterMsLevel(sp, 1))
#'
#' ## And the number of MS2 scans with precursor ions selected
#' ## from MS1 scans (those in the data and others)
#' table(precScanNum(sp))
#'
#' ## Count number of sequences/identifications per scan
#' sp <- nSequences(sp)
#'
#' ## MS2 scans either lead to an identification (5 instances) or none
#' ## (76). Among the five MS1 scans in the experiment, 3 lead to MS2
#' ## scans being matched to no peptides and two MS1 scans produced two
#' ## and three PSMs respectively.
#' table(sp$nSequences, sp$msLevel)
nSequences <- function(object,
                       sequence = "sequence",
                       f = dataStorage(object),
                       BPPARAM = bpparam()) {
    stopifnot(is(object, "Spectra"))
    if ("nSequences" %in% spectraVariables(object)) {
        message("Spectra variable 'nSequences' already present.")
        return(object)
    }
    sequence <- sequence[1]
    res <- .lapply(object,
                   FUN = .nSequences,
                   sequence = "sequence",
                   f = f,
                   BPPARAM = BPPARAM)
    object$nSequences <- unsplit(res, f)
    object
}

.nSequences <- function(object, sequence) {
    if (length(unique(dataOrigin(object))) != 1)
        stop(".nSequences must be called on data from a single dataOrigin.")
    ## Count 0 is sequence is NA, 1 if otherwise
    res <- as.integer(!is.na(object[[sequence]]))
    names(res) <- acquisitionNum(object)
    ## Iterate over MS1 acquition numbers
    N <- acquisitionNum(filterMsLevel(object, 1))
    for (i in seq_along(N)) {
        ## curent MS1 acquisition number
        an <- N[i]
        an_char <- as.character(an)
        ## get logical with all descendant acquisitions
        j <- .filterSpectraHierarchy(object$acquisitionNum,
                                     object$precScanNum, an)
        res[an_char] <- sum(!is.na(object[[sequence]][j]))
    }
    unname(res)
}
