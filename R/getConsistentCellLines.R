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
#' @param min A numeric value representing the minimum consistency desired.
#'     Default is 0.5.
#'
#'
#' @return A list of consistent cell lines.
#'
#' @author {Casey Hon, \email{casey.hon@mail.utoronto.ca}}
#'

getConsistentCellLines <- function(correlations, sensMeasure, min) {

  # Performing Checks
  if (is.list(correlations)) {

  }
}
