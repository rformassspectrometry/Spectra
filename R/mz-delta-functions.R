#' @title MZ delta Quality Control
#'
#' @aliases computeMzDeltas plotMzDelta
#'
#' @description
#'
#' The M/Z delta plot illustrates the suitability of MS2 spectra for
#' identification by plotting the M/Z differences of the most intense
#' peaks. The resulting histogram should optimally show modes at
#' amino acid residu masses. The plots have been described in Foster
#' et al. 2011.
#'
#' Only a certain percentage of most intense MS2 peaks are taken into
#' account to use the most significant signal. Default value is 20%
#' (see `percentage` argument). The difference between peaks is then
#' computed for all individual spectra and their distribution is
#' plotted as a histogram. Delta M/Z between 40 and 200 are plotted
#' by default, to encompass the residue masses of all amino acids and
#' several common contaminants, although this can be changes with the
#' `xlim` argument.
#'
#' In addition to the processing described above, isobaric reporter
#' tag peaks and the precursor peak can also be removed from the MS2
#' spectrum, to avoid interence with the fragment peaks.
#'
#' Note that figures in Foster et al. 2011 have been produced and
#' optimised for centroided data. While running the function on
#' profile mode is likely fine, it is recommended to use centroided
#' data.
#'
#' @references
#'
#' Foster JM, Degroeve S, Gatto L, Visser M, Wang R, Griss J, et al. A
#' posteriori quality control for the curation and reuse of public
#' proteomics data. Proteomics. 2011;11:
#' 2182-2194. http://dx.doi.org/10.1002/pmic.201000602
#'
#' @param object An instance of class `Spectra()`.
#'
#' @param percentage `numeric(1)` between 0 and 1 indicating the
#'     percentage of the most intense peaks in each MS2 spectrum to
#'     include in the calculation. Default is 0.2.
#'
#' @param xlim `numeric(2)` with the upper and lower M/Z to be used to
#'     the MZ deltas. Default is `c(40, 200)`.
#'
#' @param BPPARAM An optional `BiocParallelParam` instance determining
#'     the parallel back-end to be used during evaluation. Default is
#'     to use `BiocParallel::bpparam()`. See `?BiocParallel::bpparam`
#'     for details.
#'
#' @param x A list of M/Z delta values, as returned by
#'     `computeMzDeltas()`.
#'
#' @param aaLabels `logical(1)` defining whether the amino acids
#'     should be labelled on the histogram. Default is `TRUE`.
#'
#' @return `computeMzDeltas()` returns a `list` of numeric
#'     vectors. `plotMzDelta()` is used to visualise of M/Z delta
#'     distributions.
#'
#' @author Laurent Gatto with contributions (to MSnbase) of
#'     Guangchuang Yu.
#'
#' @name plotMzDelta
#'
#' @examples
#'
#' library(msdata)
#' f <- proteomics(pattern = "TMT.+20141210.mzML.gz", full.names = TRUE)
#' sp <- Spectra(f)
#'
#' d <- computeMzDeltas(sp)
#' plotMzDelta(d)
NULL


#' @importFrom stats quantile
#'
#' @noRd
compute_mz_deltas <- function(pks,
                              percentage,
                              xlim) {
    ## only keep top intensity peaks
    sel <- pks[, "intensity"] >= quantile(pks[, "intensity"], 1 - percentage)
    ## keep mz values of these top peaks
    mzs <- pks[sel, "mz"]
    ## prepare list for all delta mzs
    delta <- vector("list", length = length(mzs))
    i <- 1
    ## compute all delta mzs for 1st, 2nd, ... to last peak
    while (length(mzs) > 1) {
        m <- mzs[1]
        mzs <- mzs[-1]
        delta[[i]] <- abs(mzs - m)
        i <- i + 1
    }
    delta <- unlist(delta)
    ## only keep deltas that are relevant for aa masses
    delta[delta > xlim[1] & delta < xlim[2]]
}

#' @rdname plotMzDelta
#'
#' @import BiocParallel
#'
#' @export
computeMzDeltas <- function(object,
                            percentage = 0.2,
                            xlim = c(40, 200),
                            BPPARAM = BiocParallel::bpparam()) {
    BiocParallel::bplapply(peaksData(filterMsLevel(object, 2)),
                           compute_mz_deltas,
                           percentage = percentage,
                           xlim = xlim,
                           BPPARAM = BPPARAM)
}


#' @rdname plotMzDelta
#'
#' @importFrom graphics hist segments
#'
#' @export
plotMzDelta <- function(x, aaLabels = TRUE) {
    ## from PSM::getAminoAcids()
    amino_acids <-
        structure(list(AA = c("peg", "A", "R", "N", "D", "C", "E", "Q",
                              "G", "H", "I", "L", "K", "M", "F", "P", "S",
                              "T", "W", "Y", "V"),
                       ResidueMass = c(44, 71.03711, 156.10111, 114.04293,
                                       115.02694, 103.00919, 129.04259,
                                       128.05858, 57.02146, 137.05891,
                                       113.08406, 113.08406, 128.09496,
                                       131.04049, 147.06841, 97.05276,
                                       87.03203, 101.04768, 186.07931,
                                       163.06333, 99.06841),
                       Abbrev3 = c(NA, "Ala", "Arg", "Asn", "Asp", "Cys",
                                   "Glu", "Gln", "Gly", "His", "Ile",
                                   "Leu", "Lys", "Met", "Phe", "Pro",
                                   "Ser", "Thr", "Trp", "Tyr", "Val"),
                       ImmoniumIonMass = c(NA, 44.05003, 129.114,
                                           87.05584, 88.03986, 76.0221,
                                           102.0555, 101.0715, 30.03438,
                                           110.0718, 86.09698, 86.09698,
                                           101.1079, 104.0534, 120.0813,
                                           70.06568, 60.04494, 74.06059,
                                           159.0922, 136.0762, 72.08133),
                       Name = c("Polyethylene glycol", "Alanine",
                                "Arginine", "Asparagine", "Aspartic acid",
                                "Cysteine", "Glutamic acid", "Glutamine",
                                "Glycine", "Histidine", "Isoleucine",
                                "Leucine", "Lysine", "Methionine",
                                "Phenylalanine", "Proline", "Serine",
                                "Threonine", "Tryptophan", "Tyrosine",
                                "Valine"),
                       Hydrophobicity = c(NA, 0.62, -2.53, -0.78, -0.9,
                                          0.29, -0.74, -0.85, 0.48, -0.4,
                                          1.38, 1.06, -1.5, 0.64, 1.19,
                                          0.12, -0.18, -0.05, 0.81, 0.26,
                                          1.08),
                       Hydrophilicity = c(NA, -0.5, 3, 0.2, 3, -1, 3, 0.2,
                                          0, -0.5, -1.8, -1.8, 3, -1.3,
                                          -2.5, 0, 0.3, -0.4, -3.4, -2.3,
                                          -1.5),
                       SideChainMass = c(NA, 15, 101, 58, 59, 47, 73, 72,
                                         1, 82, 57, 57, 73, 75, 91, 42,
                                         31, 45, 130, 107, 43),
                       pK1 = c(NA, 2.35, 2.18, 2.18, 1.88, 1.71, 2.19,
                               2.17, 2.34, 1.78, 2.32, 2.36, 2.2, 2.28,
                               2.58, 1.99, 2.21, 2.15, 2.38, 2.2, 2.29),
                       pK2 = c(NA, 9.87, 9.09, 9.09, 9.6, 10.78, 9.67,
                               9.13, 9.6, 8.97, 9.76, 9.6, 8.9, 9.21,
                               9.24, 10.6, 9.15, 9.12, 9.39, 9.11, 9.74),
                       pI = c(NA, 6.11, 10.76, 10.76, 2.98, 5.02, 3.08,
                              5.65, 6.06, 7.64, 6.04, 6.04, 9.47, 5.74,
                              5.91, 6.3, 5.68, 5.6, 5.88, 5.63, 6.02)),
                  class = "data.frame",
                  row.names = c(NA, -21L))
    col <- "grey30"
    hist(unlist(x), breaks = 200, col = col, border = col,
         ylab = "Frequency", xlab = "M/Z delta",
         main = "Histrogram of Mass Delta Distributions")
    if (aaLabels) {
        usr <- par("usr")
        aa_labels <- amino_acids$AA
        aa_labels[aa_labels == "I"] <- "I/L"
        aa_labels[aa_labels == "L"] <- ""
        aa_labels[aa_labels == "Q"] <- "Q/K"
        aa_labels[aa_labels == "K"] <- ""
        amino_acids$label <- aa_labels
        amino_acids$x_adj <- 1
        amino_acids$y_adj <- 0.97
        amino_acids$y_adj[amino_acids$label == "I/L"] <- 0.95
        amino_acids$x_adj[amino_acids$label == "I/L"] <- 0.985
        amino_acids$y_adj[amino_acids$label == "D"] <- 0.95
        amino_acids$x_adj[amino_acids$label == "D"] <- 1.008
        amino_acids$y_adj[amino_acids$label == "Q/K"] <- 0.95
        amino_acids$x_adj[amino_acids$label == "Q/K"] <- 0.985
        amino_acids$y_adj[amino_acids$label == "M"] <- 0.95
        amino_acids$x_adj[amino_acids$label == "M"] <- 1.0075
        segments(amino_acids$ResidueMass, 0,
                 amino_acids$ResidueMass, usr[4] * 0.95,
                 col = "red", lwd = 1, lty =  "dashed")
        text(amino_acids$ResidueMass  * amino_acids$x_adj,
             usr[4] * amino_acids$y_adj,
             amino_acids$label)
    }
}
