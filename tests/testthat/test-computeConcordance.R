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

test_that("Check if function returns right concordance object" {

  allCorrelations = drugAllCorrelations
  subsettedCorrelations = drugConsistentCorrelations
  sensMeasure = "aac_recomputed_corrs"
  coefName = "pearson"

  concordance <- computeConcordance(allCorrelations = allCorrelations,
                                    subsettedCorrelations = subsettedCorrelations,
                                    sensMeasure = sensMeasure,
                                    coefName = coefName)





})
