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
#' @param sensMeasures List of drug sensitivity measures of interest. Options
#'     depend on the PharmacoSets chosen, and can be identified by looking at
#'     names(pSetName@sensitivity[["profiles"]])
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
  ## TO-DO: check if pSet is valid (contains a list of at least two pSets)


  if (is.character(coefs) == TRUE) {
    coefsUsed <- c("pearson", "spearman", "kendall")
    if (all((coefs == coefsUsed) == FALSE)) {
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
  sensUsed <- list()
  for (set in pSet) {
    x <- set@annotation[["name"]]
    pSetList <- append(pSetList, x)
    sensUsed <- append(sensUsed, names(set@sensitivity[["profiles"]]))
    assign(x, set, envir = globalenv())
  }

  sensUsed <- unique(sensUsed) # keep only unique values of sensitivity measures

  if (is.character(sensMeasures) == TRUE) {
    if (all((sensMeasures == sensUsed) == FALSE)) {
      stop("sensMeasures should be of class character, specifying the drug
         sensitivity measures as outputted by
         > names(PSET@sensitivity[[\"profiles\"]])")
    }
  } else if (is.character(sensMeasures) != TRUE) {
    stop("sensMeasures should be of class character, specifying the drug
         sensitivity measures as outputted by
         > names(PSET@sensitivity[[\"profiles\"]])")
  }

  # get number of pairs,combinations of PSets to be compared
  # (each column is a pair)
  pSetPairs <- as.data.frame(combn(pSetList, 2))

  # include p-values of each correlation coefficient
  if (isTRUE(pval)) {
    coefs <- c(coefs, paste(coefs, "p-value", sep = ' '))
  }


  # initialize data frame to be outputted
  cors = data.frame(matrix(ncol = length(coefs) * ncol(pSetPairs),
                           nrow = length(colnames(intersected)),
                           row.names = colnames(intersected)))
  colnames(cors) <- coefs

  # get sensitivity profiles for each sensitivity measure for each pset
  for (pSet in pSet) {
    for (sensName in sensMeasures) {
      var_name <- paste0(coefName, pSet@annotation[["name"]], sep = "_")
      assign(var_name,
             as.data.frame(PharmacoGx::summarizeSensitivityProfiles(pSet,
                                                                    sensitivity.measure = coefName)),
             envir = globalenv())
    }
  } ## TO-DO: not sure if this works yet

  # compute each correlation in coef between each pair of psets
  for (cell.line in rownames(cors)) {
    tryCatch(
      expr = {
        if ("pearson" %in% coefs) {
          for (pair in colnames(pSetPairs)) {
            pearson.cor <- cor.test(x = pSetPairs["1", pair], y = pSetPairs["2", pair]) # not finished
          }

        }
      }
    )
  }





}
