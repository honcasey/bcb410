library(testthat)
library(cellLineConsistency)

# Generate example data - takes a while so putting out here
intersected <- PharmacoGx::intersectPSet(c(CTRPv2, GRAY),
                                         intersectOn = c("drugs", "cell.lines"))
correlations <- computeCellLineCorrelation(pSet = intersected,
                                           coefs = "pearson",
                                           sensMeasures = "aac_recomputed",
                                           pval = TRUE)
consistentLines <- getConsistentCellLines(correlations,
                                          sensMeasure = "aac_recomputed_corrs",
                                          coefName = "pearson")
drugAllCorrelations <- computeDrugCorrelation(pSet = intersected,
                                              coefs = "pearson",
                                              sensMeasures = "aac_recomputed")
drugConsistentCorrelations <- computeDrugCorrelation(pSet = intersected,
                                                     cellLines = rownames(consistentLines),
                                                     coefs = "pearson",
                                                     sensMeasures = "aac_recomputed")

test_that("Check if function returns right concordance object", {

  allCorrelations = drugAllCorrelations
  subsettedCorrelations = drugConsistentCorrelations
  sensMeasure = "aac_recomputed_corrs"
  coefName = "pearson"

  concordance <- computeConcordance(allCorrelations = allCorrelations,
                                    subsettedCorrelations = subsettedCorrelations,
                                    sensMeasure = sensMeasure,
                                    coefName = coefName)

  expect_s3_class(concordance, "concordance")
  expect_length(concordance, 7)
})

test_that("check for invalid user inputs", {
  allCorrelations = drugAllCorrelations
  subsettedCorrelations = drugConsistentCorrelations
  sensMeasure = "aac_recomputed_corrs"
  coefName = "pearson"

  # correlation parameters is not a list
  random_df <- data.frame()
  expect_error(computeConcordance(allCorrelations = random_df,
                                  subsettedCorrelations = subsettedCorrelations,
                                  sensMeasure = sensMeasure,
                                  coefName = coefName),
               "All items in allCorrelations should be dataframes containing
           computed correlations for your sensitivity measure.")

  # correlation doesn't contain dataframes
  random_list <- list()
  expect_error(computeConcordance(allCorrelations = random_list,
                                  subsettedCorrelations = subsettedCorrelations,
                                  sensMeasure = sensMeasure,
                                  coefName = coefName),
               "All items in allCorrelations should be dataframes containing
           computed correlations for your sensitivity measure.")

  # sensMeasure not of class character
  expect_error(computeConcordance(allCorrelations = allCorrelations,
                                  subsettedCorrelations = subsettedCorrelations,
                                  sensMeasure = 123,
                                  coefName = coefName),
               "sensMeasure must be of class character specifying one of the
         measures included in correlations.")

  # sensMeasure is not in list
  expect_error(computeConcordance(allCorrelations = allCorrelations,
                                  subsettedCorrelations = subsettedCorrelations,
                                  sensMeasure = "auc_published",
                                  coefName = coefName),
               "sensMeasure must be one of the measures included in correlations.")

  # coefName is not in list
  expect_error(computeConcordance(allCorrelations = allCorrelations,
                                  subsettedCorrelations = subsettedCorrelations,
                                  sensMeasure = sensMeasure,
                                  coefName = "parson"),
               "coefName must be one of the correlation coefficients included
           in correlations.")

  # coefName is not of class character
  expect_error(computeConcordance(allCorrelations = allCorrelations,
                                  subsettedCorrelations = subsettedCorrelations,
                                  sensMeasure = sensMeasure,
                                  coefName = 123),
               "coefName must be of class character specifying one of the correlation
         coefficients included in correlations.")
})
