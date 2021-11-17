#' Intersect PharmacoSets on common sensitivity measures
#'
#' Given an object of intersected PharmacoSets, the function finds sensitivity
#' measures in common.
#'
#' @param pSet intersected PharmacoSets of interest
#'
#' @return Returns a list of sensitivity measures in common.
#'
#' @examples
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
#' @references
#' Smirnov P, Safikhani Z, El-Hachem N, Wang D, She A, Olsen C, Freeman M,
#'     Selby H, Gendoo D, Grossman P, Beck A, Aerts H, Lupien M,
#'     Haibe-Kains AG, (2016). PharmacoGx: an R package for analysis of
#'     large pharmacogenomic datasets. \emph{Bioinformatics (Oxford, England)}.
#'
#' @import PharmacoGx
#'

intersectSensMeasures <- function(pSet) {
  # Performing checks
  if (is.list(pSet) == TRUE) {
    if (length(pSet) < 2) {
      stop("pSet list should have at least two PharmacoSets to
           compare correlations.")
    }
  } else if (is.list(pSet) == FALSE) {
    stop("pSet should be of class list, which contains a list of
         intersected PharmacoSets.")
  }

  measList <- list()
  for (set in pSet) {
    setList <- list(names(set@sensitivity[["profiles"]]))
    measList <- append(measList, setList)
  }
  return(Reduce(intersect, measList))
}

