#' Intersect PharmacoSets on common sensitivity measures
#'
#' Given an object of intersected PharmacoSets, the function finds sensitivity
#' measures in common.
#' @param pSet intersected PharmacoSets of interest
#'
#' @return Returns a list of sensitivity measures in common.
#'
#' @example
#' CTRP <- PharmacoGx::downloadPSet("CTRPv2_2015")
#' GRAY <- PharmacoGx::downloadPSet("GRAY_2013")
#' intersected <- PharmacoGx::intersectPSet(c(CTRP, GRAY), intersectOn = c("drugs", "cell.lines"))
#'
#' commonSensMeasures <- intersectSensMeasures(intersected)
#' commonSensMeasures
#' > [1] "aac_recomputed"  "ic50_recomputed" "HS"              "E_inf"           "EC50"
#'
#' @author {Casey Hon, \email{casey.hon@mail.utoronto.ca}}
#'


library(PharmacoGx)

intersectSensMeasures <- function(pSet) {
  measList <- list()
  for (set in pSet) {
    setList <- list(names(set@sensitivity[["profiles"]]))
    measList <- append(measList, setList)
  }
  return(Reduce(intersect, measList))
}
