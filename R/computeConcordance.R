#' Compute consistency between correlations of cell line sensitivity measures
#'
#' Computes Harrell's Concordance Index between drug sensitivity measures
#' across PharmacoSets after removing inconsistent cell lines.
#'
#' @param allCorrelations A list of dataframes that contain the correlation
#'     coefficients across drugs for all
#'     cell lines, as returned by computeDrugCorrelation().
#' @param subsettedCorrelations A list of dataframes that contain the
#'     correlation coefficients across drugs
#'     for only subsetted cell lines, as returned by computeDrugCorrelation().
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
#' intersected <- PharmacoGx::intersectPSet(c(CTRP, GRAY),
#'     intersectOn = c("drugs", "cell.lines"))
#' correlations <- computeCellLineCorrelation(intersected,
#'     coefs = c("pearson", "spearman"), TRUE)
#' consistentLines <- getConsistentCellLines(correlations,
#'     sensMeasure = "aac_recomputed_corrs", coefName = "pearson")
#' drugAllCorrelations <- computeDrugCorrelation(pSet = intersected,
#'     coefs = "pearson",
#'     sensMeasures = "aac_recomputed")
#' drugConsistentCorrelations <- computeDrugCorrelation(pSet = intersected,
#'     cellLines = consistentLines,
#'     coefs = "pearson",
#'     sensMeasures = "aac_recomputed")
#' concordance <- computeConcordance(allCorrelations = drugAllCorrelations,
#'     subsettedCorrelations = drugConsistentCorrelations,
#'     sensMeasure = "aac_recomputed",
#'     coefName = "pearson")
#'
#' @author {Casey Hon, \email{casey.hon@mail.utoronto.ca}}
#'
#' @import survival
#'

computeConcordance <- function(allCorrelations,
                               subsettedCorrelations,
                               sensMeasure,
                               coefName) {

    # Performing Checks
  if (is.list(allCorrelations) == TRUE) {
    if (all(lapply(allCorrelations, class) != "data.frame")) {
      stop("All items in allCorrelations should be dataframes containing
           computed correlation coefficients for your sensitivity
           measure of interest, as returned by computeCorrelation().")
    }
  } else if (is.list(allCorrelations) == FALSE) {
    stop("allCorrelations should be of class list containing dataframes
         for each sensitivity measure of interest, as returned by
         computeCorrelation().")
  }

  if (is.list(subsettedCorrelations) == TRUE) {
    if (all(lapply(subsettedCorrelations, class) != "data.frame")) {
      stop("All items in subsettedCorrelation should be dataframes containing
           computed correlation coefficients for your sensitivity
           measure of interest, as returned by computeCorrelation().")
    }
  } else if (is.list(subsettedCorrelations) == FALSE) {
    stop("subsettedCorrelation should be of class list containing dataframes
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
  # TO-DO: check if sensMeasure exists in both correlation dataframes

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


  toSurv <- transform(merge(allCorrelations[sensMeasure],
                            subsettedCorrelations[sensMeasure],
                            by = 0,
                            all = TRUE), row.names=Row.names, Row.names=NULL)
  concorded <- survival::concordance(object = paste0(sensMeasure, "pearson.x") ~ paste0(sensMeasure, "pearson.y"),
                                   data = toSurv)

}






