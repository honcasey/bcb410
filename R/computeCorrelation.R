#' Compute correlations between drug sensitivity measures of cell lines
#'
#' Computes Pearson and Spearman Correlation Coefficients for specified
#'     drug sensitivity measures (such as IC50 of AUC) and specified
#'     subset of cell lines
#' @param pSet Intersected PharmacoSet object containing common drugs and/or
#'     cell lines, as returned by PharmacoGx::intersectPSet(). PharmacoSets
#'     should have sensitivity information which can be checked using
#'     PharmacoGx::availablePSets() in the 'type' column.
#' @param coefs List of correlation coefficients of interest. Current options
#'     are: "pearson", "spearman", "kendall".
#' @param pval Logical of whether or not to include p-values of each
#'     correlation coefficient. Default is TRUE.
#'
#' @return Returns a dataframe where rows correspond to drugs and columns
#'     correspond to each type of correlation coefficient from coefs.
#'     Also includes p-values in columns is pval was set to TRUE.
#'
#' @example
#' # Intersect PharmacoSets of interest based on common cell lines
#' CTRP <- PharmacoGx::downloadPSet("CTRPv2_2015")
#' GRAY <- PharmacoGx::downloadPSet("GRAY_2013")
#' intersected <- PharmacoGx::intersectPSet(c(CTRP, GRAY), intersectOn = c("drugs", "cell.lines"))
#' rm(CTRP)
#' rm(GRAY)
#'
#' set.seed(1001)
#'
#'
#' @author {Casey Hon, \email{casey.hon@mail.utoronto.ca}}
#'
#' @references
#'
#'

install.packages(c("devtools", "tidyverse", "fs", "PharmacoGx"))
library(devtools)
library(tidyverse)
library(fs)
library(PharmacoGx)

computeCorrelation <- function(pSet,
                               coefs,
                               pval) {

  # Performing checks


  if (is.character(coefs) == TRUE) {
    coefsUsed <- c("pearson", "spearman", "kendall")
    if(all((coefs == coefsUsed) == FALSE)) {
      stop("coefs should be of class character, specifying either:
           pearson, spearman, or kendall.")
    }
  } else if (is.character(coefs) != TRUE) {
      stop("coefs should be of class character, specifying either:
           pearson, spearman, or kendall.")
  }

  if (is.logical(pval) != TRUE) {
    stop("pval should be of class logical, specifying either TRUE or FALSE.")
  }

  # separate each pset after intersection
  pSetList <- list()
  for (set in pSet) {
    x <- set@annotation[["name"]]
    pSetList <- append(pSetList, x)
    assign(x, set, envir = globalenv())
  }

  # get number of pairs of PSets to be compared, each column is a pair
  pSetPairs <- combn(pSetList, 2)

  # include p-values of each correlation coefficient
  if (isTRUE(pval)) {
    coefs <- c(coefs, paste(coefs, "p-value", sep = ' '))
  }


  # initialize data frame to be outputted
  cors = data.frame(matrix(ncol = length(coefs) * ncol(pSetPairs),
                           nrow = length(colnames(intersected)),
                           row.names = colnames(intersected)))
  colnames(cors) <- coefs
  for (cell.line in rownames(cors)) {
    tryCatch(
      expr = {
      }
    )
  }





}
