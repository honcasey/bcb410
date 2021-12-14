#' Get list of consistent cell lines
#'
#' Identifies cell lines that are deemed consistent based on the
#'     user-specified minimum correlation of sensitivity measures
#'     across PharmacoSets.
#'
#' @param correlations A list of dataframes that contain the correlation
#'     coefficients of the sensitivity measures of interest, as returned
#'     by computeCorrelation().
#' @param sensMeasure A character vector defining the sensitivity measure
#'     of interest. Must be one of the measures included in the correlations
#'     dataframe.
#' @param coefName A character vector defining the correlation coefficient
#'     of interest. Must be one of the coefficients included in correlations
#'     dataframe.
#' @param min A numeric value representing the minimum correlation coefficient
#'     desired. All correlations greater than or equal to min are deemed
#'     consistent. Default is 0.5.
#'
#' @return A dataframe containing only correlations of consistent cell lines.
#'
#' @examples
#' # Intersect PharmacoSets of interest based on common cell lines
#' CTRP <- PharmacoGx::downloadPSet("CTRPv2_2015")
#' GRAY <- PharmacoGx::downloadPSet("GRAY_2013")
#' intersected <- PharmacoGx::intersectPSet(c(CTRP, GRAY), intersectOn = c("drugs", "cell.lines"))
#' correlations <- computeCellLineCorrelation(pSet = intersected,
#'     coefs = "pearson",
#'     sensMeasures = "aac_recomputed",
#'     pval = TRUE)
#' consistentLines <- getConsistentCellLines(correlations,
#'     sensMeasure = "aac_recomputed_corrs",
#'     coefName = "pearson")
#'
#' @author {Casey Hon, \email{casey.hon@mail.utoronto.ca}}
#'
#' @export
#'

getConsistentCellLines <- function(correlations,
                                   sensMeasure,
                                   coefName,
                                   min = 0.5) {

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

  if (is.character(coefName) == TRUE) {
    coefUsed <- names(correlations[[sensUsed]])
    if (all(coefName %in% coefUsed) == FALSE) {
      stop("coefName must be one of the correlation coefficients included
               in correlations.")
    }
  } else if (is.character(coefName) == FALSE) {
    stop("coefName must be of class character specifying one of the correlation
         coefficients included in correlations.")
  }

  if (missing(min) == TRUE) {
    min = 0.5 # DEFAULT
  } else if (!missing(min)) {
    if (!is.numeric(min)) {
      stop("Min should be a numeric value.")
    } else if (is.numeric(min)) {
      if (min < 0 || min > 1) {
        stop("Min should be a numeric value from 0 to 1.")
      }
    }
  }

  sens_cor <- correlations[[sensMeasure]]   # subset correlations
  # only keep correlations geq to min
  cons <- sens_cor[which(sens_cor[coefName] >= min), ]

  return(cons)
}
#[END]
