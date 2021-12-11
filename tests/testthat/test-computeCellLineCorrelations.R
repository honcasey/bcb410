library(testthat)
library(cellLineConsistency)

# Generate example data - takes a while so putting out here
intersected <- PharmacoGx::intersectPSet(c(CTRPv2, GRAY),
                                         intersectOn = c("drugs", "cell.lines"))

test_that("Checking function with one sensitivity measure and one coefficient parameters", {

  pSet = intersected
  coefs = c("pearson")
  sensMeasures = c("aac_recomputed")
  pval = TRUE
  correlations <- computeCellLineCorrelation(pSet = pSet,
                                     coefs = coefs,
                                     sensMeasures = sensMeasures,
                                     pval = pval)

  sensName <- paste0(sensMeasures, "_corrs")

  expect_type(correlations, "list")
  expect_type(correlations[[sensName]][[coefs]], "double")
  expect_length(correlations, length(coefs))
  expect_length(correlations[[sensName]][[coefs]], 34)
  expect_equal(names(correlations[1]), sensName)
  expect_type(correlations[[sensName]], "list")
  expect_length(correlations[[sensName]], 2)
})

test_that("Checking function with two sensitivity measures and two coefficient parameters", {

  pSet = intersected
  coefs = c("pearson", "spearman")
  sensMeasures = c("aac_recomputed", "ic50_recomputed")
  pval = TRUE
  correlations <- computeCellLineCorrelation(pSet = pSet,
                                             coefs = coefs,
                                             sensMeasures = sensMeasures,
                                             pval = pval)

  sensName <- paste0(sensMeasures, "_corrs")

  expect_type(correlations, "list")
  expect_length(correlations, length(coefs))
  expect_identical(names(correlations), sensName)
  expect_type(correlations[[sensName[1]]], "list")
  expect_type(correlations[[sensName[2]]], "list")
  expect_length(correlations[[sensName[1]]], 4)
  expect_length(correlations[[sensName[2]]], 4)

})

test_that("Checking function with FALSE pval ", {

  pSet = intersected
  coefs = c("pearson")
  sensMeasures = c("aac_recomputed")
  pval = FALSE
  correlations <- computeCellLineCorrelation(pSet = pSet,
                                             coefs = coefs,
                                             sensMeasures = sensMeasures,
                                             pval = pval)

  sensName <- paste0(sensMeasures, "_corrs")

  expect_type(correlations, "list")
  expect_type(correlations[[sensName]][[coefs]], "double")
  expect_length(correlations, length(coefs))
  expect_length(correlations[[sensName]][[coefs]], 34)
  expect_equal(names(correlations[1]), sensName)
  expect_type(correlations[[sensName]], "list")
  expect_length(correlations[[sensName]], 1)

})

test_that("Checking for invalid user inputs", {

  pSet = intersected
  coefs = c("pearson")
  sensMeasures = c("aac_recomputed")
  pval = TRUE

  # pSet not provided in list
  expect_error(computeCellLineCorrelation(pSet = CTRPv2,
                                          coefs = coefs,
                                          sensMeasures = sensMeasures,
                                          pval = pval),
               "pSet should be of class list, which contains a list of
         intersected PharmacoSets.")

  # only 1 pSet provided
  expect_error(computeCellLineCorrelation(pSet = c(CTRPv2),
                                         coefs = coefs,
                                         sensMeasures = sensMeasures,
                                         pval = pval),
               "pSet list should have at least two PharmacoSets to
           compare correlations.")

  # coef spelled wrong
  expect_error(computeCellLineCorrelation(pSet = intersected,
                                          coefs = c("parson"),
                                          sensMeasures = sensMeasures,
                                          pval = pval),
               "coefs should be of class character, specifying either:
           pearson, spearman, and/or kendall.")

  # sensMeasures spelled wrong
  expect_error(computeCellLineCorrelation(pSet = intersected,
                                          coefs = coefs,
                                          sensMeasures = "abc_recomputed",
                                          pval = pval),
               "sensMeasures should be of class character, specifying the drug
         sensitivity measures. The list must be a subset of the measures in
           common to all PSets intersected.")

  # pval provided not a logical
  expect_error(computeCellLineCorrelation(pSet = intersected,
                                          coefs = coefs,
                                          sensMeasures = sensMeasures,
                                          pval = "tru"),
               "pval should be of class logical, specifying either TRUE or FALSE.")

})

## TO-DO: tests for more than 2 psets
