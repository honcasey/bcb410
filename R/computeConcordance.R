#' Compute consistency between correlations of cell line sensitivity measures
#'
#' Computes Harrell's Concordance Index between drug sensitivity measures
#' across PharmacoSets after removing inconsistent cell lines.
#'
#' @param correlations A list of dataframes that contain the correlation
#'     coefficients of the sensitivity measures of interest, as returned
#'     by computeCorrelation().
#' @param cellLines A list of cell lines deemed consistent, as returned
#'     by getConsistentCellLines().
#' @param sensMeasure A character vector specifying the sensitivity
#'     measure of interest. Must be one of the measures included in
#'     correlations parameter.
#' @param coef A character vector specifying the correlation coefficient
#'     of interest. Must be one of the coefficients included in the
#'     correlations parameter.
#'
#' @return Returns a dataframe outlining concordance indices, where rows
#'     are ??? and columns are ???
#'
#' @author {Casey Hon, \email{casey.hon@mail.utoronto.ca}}
#'



computeConcordance <- function(correlations, cellLines, sensMeasure, coef) {

}






