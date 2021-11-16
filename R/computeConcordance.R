#' Compute consistency between correlations of cell line sensitivity measures
#'
#' Computes Harrell's Concordance Index between drug sensitivity measures
#' across PharmacoSets after removing inconsistent cell lines. Harrell’s
#' C-index, also known as the concordance index, is a goodness of fit
#' measure commonly used to evaluate risk models in survival analysis.
#' To measure the concordance, the number of concordant and discordant pairs
#' is used.  Concordant and discordant pairs are used to describe the
#' relationship between pairs of observations, which in this case is the pair
#' of the sensitivity measures reported in each of the two studies. A value
#' between 0 and 1 is returned, where a value below 0.5
#' indicates a very poor model, a value of 0.5 indicates that the model
#' predicts at the same accuracy as random chance, and a value of 1
#' indicates that the model’s prediction accuracy is perfect. The concordance
#' serves as a means of quantifying the improvement in consistency of drug
#' sensitivity measures.
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
#' @return Returns an S3 concordance object outlining concordance between
#'     the coefName correlations of sensMeasure in allCorrelations and in
#'     subsettedCorrelations.
#'
#' @examples
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
#'     sensMeasure = "aac_recomputed_corrs",
#'     coefName = "pearson")
#'
#' @author {Casey Hon, \email{casey.hon@mail.utoronto.ca}}
#'
#' @references
#' Harrell FE, Califf RM, Pryor DB, Lee KL, Rosati RA, (1982). Evaluating the
#'     Yield of Medical Tests. \emph{JAMA}, 247(18):2543–2546.
#'
#' Statistical Odds and Ends (2019). What is Harrell's C-index?
#'     https://statisticaloddsandends.wordpress.com/2019/10/26/what-is-harrells-c-index/
#'
#' Therneau T (2021). A Package for Survival Analysis in R. R package version
#'     3.2-13, https://CRAN.R-project.org/package=survival.
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
    sensUsed <- names(allCorrelations)
    if ((sensMeasure %in% sensUsed) == FALSE) {
      stop("sensMeasure must be one of the measures included in correlations
           as returned by names(allCorrelations).")
    }
  } else if (is.character(sensMeasure) == FALSE) {
    stop("sensMeasure must be of class character specifying one of the
         measures included in correlations as returned by
         names(allCorrelations).")
  }
  # TO-DO: check if sensMeasure exists in both correlation dataframes

  if (is.character(coefName) == TRUE) {
    coefUsed <- names(allCorrelations[[sensMeasure]])
    if (all((coefName == coefUsed) == FALSE)) {
      stop("coefName must be one of the correlation coefficients included
           in correlations, as returned by
           names(allCorrelations[[\"sensMeasure\"]]).")
    }
  } else if (is.character(coefName) == FALSE) {
    stop("coefName must be of class character specifying one of the correlation
         coefficients included in correlations, as returned by
         names(allCorrelations[[\"sensMeasure\"]]).")
  }

  toSurv <- transform(merge(allCorrelations[sensMeasure],
                            subsettedCorrelations[sensMeasure],
                            by = 0,
                            all = TRUE), row.names=Row.names, Row.names=NULL)

  concorded <- survival::concordance(object = get(paste0(sensMeasure, ".pearson.x"))
                                     ~ get(paste0(sensMeasure, ".pearson.y")),
                                   data = toSurv)

  return(concorded)
}

