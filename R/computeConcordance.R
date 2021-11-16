#' Compute consistency between correlations of cell line sensitivity measures
#'
#' Computes Harrell's Concordance Index between drug sensitivity measures
#' across PharmacoSets after removing inconsistent cell lines.
#'
#' @param pSets A character vector outlining the two PharmacoSets from which
#'     the correlations were derived from.
#' @param correlations A list of dataframes that contain the correlation
#'     coefficients of the sensitivity measures of interest, as returned
#'     by computeCorrelation().
#' @param cellLines A list of cell lines deemed consistent, as returned
#'     by getConsistentCellLines().
#' @param sensMeasure A character vector specifying the sensitivity
#'     measure of interest. Must be one of the measures included in
#'     correlations parameter.
#' @param coefName A character vector specifying the correlation coefficient
#'     of interest. Must be one of the coefficients included in the
#'     correlations parameter.
#'
#' @return Returns a dataframe outlining concordance indices, where rows
#'     are ??? and columns are ???
#'
#' @example
#' # Intersect PharmacoSets of interest based on common cell lines
#' CTRP <- PharmacoGx::downloadPSet("CTRPv2_2015")
#' GRAY <- PharmacoGx::downloadPSet("GRAY_2013")
#' intersected <- PharmacoGx::intersectPSet(c(CTRP, GRAY), intersectOn = c("drugs", "cell.lines"))
#' correlations_CTRPvGRAY <- computeCorrelation(intersected, coefs = c("pearson", "spearman"), TRUE)
#' consistentLines_CTRPvGRAY <- getConsistentCellLines(correlations,
#'     "aac_recomputed_corrs", "pearson")
#' concordance_CTRPvGRAY <- computeConcordance(c("CTRP", "GRAY"), correlations,
#'     rownames(consistentLines), "aac_recomputed", "pearson")
#'
#' @author {Casey Hon, \email{casey.hon@mail.utoronto.ca}}
#'
#' @import survival
#'



computeConcordance <- function(correlations, cellLines, sensMeasure, coef) {

    # Performing Checks
  if (is.list(correlations) == TRUE) {
    if (all(lapply(correlations, class) != "data.frame")) {
      stop("All items in correlations should be dataframes containing
           computed correlation coefficients for your sensitivity
           measure of interest, as returned by computeCorrelation().")
    }
  } else if (is.list(correlations) == FALSE) {
    stop("Correlations should be of class list containing dataframes
         for each sensitivity measure of interest, as returned by
         computeCorrelation().")
  }

  if (is.character(sensMeasure) == TRUE) {
    sensUsed <- names(correlations)
    if ((sensMeasure %in% sensUsed) == FALSE) {
      stop("sensMeasure must be one of the measures included in correlations
           as returned by names(correlations).")
    }
  } else if (is.character(sensMeasure) == FALSE) {
    stop("sensMeasure must be of class character specifying one of the
         measures included in correlations as returned by
         names(correlations).")
  }

  if (is.character(cellLines) == TRUE) {
    linesUsed <- rownames(correlations[[sensMeasure]])
    if (all((cellLines == linesUsed) == FALSE)) {
      stop("cellLines must be a list of class character specifying at least
           one cell line as returned by
           rownames(correlations[[sensMeasure]]).")
    }
  } else if (is.character(cellLines) == FALSE) {
    stop("cellLines must be a list of class character.")
  }

  if (is.character(coefName) == TRUE) {
    coefUsed <- names(correlations[[sensUsed]])
    if (all(coefName == coefUsed) == FALSE) {
      stop("coefName must be one of the correlation coefficients included
           in correlations, as returned by
           names(correlations[[\"sensMeasure\"]]).")
    }
  } else if (is.character(coefName) == FALSE) {
    stop("coefName must be of class character specifying one of the correlation
         coefficients included in correlations, as returned by
         names(correlations[[\"sensMeasure\"]]).")
  }

  # subset
  origCor <- correlations[[sensMeasure]]
  whichCor <- origCor[cellLines, ] # only keep consistent cell lines
  toConcord <- merge(origCor, whichCor, by = 0, all = TRUE)



  # get concordance

  concord <- survival::concordance(object = pearson.x ~ pearson.y,
                                   data = toConcord)

}






