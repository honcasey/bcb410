library(testthat)
library(cellLineConsistency)

# Generate example data - takes a while so putting out here
intersected <- PharmacoGx::intersectPSet(c(CTRPv2, GRAY),
                                         intersectOn = c("drugs", "cell.lines"))
correlations <- computeCellLineCorrelation(pSet = intersected,
                                           coefs = "pearson",
                                           sensMeasures = "aac_recomputed",
                                           pval = TRUE)

test_that("Check if function returns right cell lines", {
  correlations = correlations
  sensMeasure = "aac_recomputed_corrs"
  coefName = "pearson"
  consistentLines <- getConsistentCellLines(correlations = correlations,
                                            sensMeasure = sensMeasure,
                                            coefName = coefName)

  expect_type(consistentLines, "list")
  expect_true(coefName %in% names(consistentLines))
  expect_length(rownames(consistentLines), 17)

})

test_that("Check for invalid user inputs", {
  correlations = correlations
  sensMeasure = "aac_recomputed_corrs"
  coefName = "pearson"

  # correlations is not a list
  expect_error(getConsistentCellLines(correlations = "correlations",
                                      sensMeasure = sensMeasure,
                                      coefName = coefName),
               "Correlations should be of class list containing dataframes
         for each sensitivity measure of interest, as returned by
         computeCorrelation().")

  # correlations doesn't contain dataframes
  random_list <- list()
  expect_error(getConsistentCellLines(correlations = random_list,
                                      sensMeasure = sensMeasure,
                                      coefName = coefName),
               "All items in correlations should be dataframes containing
           computed correlation coefficients for your sensitivity
           measure of interest, as returned by computeCorrelation().")


  # coefName not in correlations
  expect_error(getConsistentCellLines(correlations = correlations,
                                      sensMeasure = sensMeasure,
                                      coefName = "auc_recomputed"),
               "coefName must be one of the correlation coefficients included
               in correlations.")

  # coefName not a character
  expect_error(getConsistentCellLines(correlations = correlations,
                                      sensMeasure = sensMeasure,
                                      coefName = 123),
               "coefName must be of class character specifying one of the correlation
         coefficients included in correlations.")

  # min not a numeric value
  expect_error(getConsistentCellLines(correlations = correlations,
                                      sensMeasure = sensMeasure,
                                      coefName = coefName,
                                      min = "0.5"),
               "Min should be a numeric value from 0 to 1.")

  # min not in bounds
  expect_error(getConsistentCellLines(correlations = correlations,
                                      sensMeasure = sensMeasure,
                                      coefName = coefName,
                                      min = 1.5),
               "Min should be a numeric value from 0 to 1.")

})
