#' Plot correlation coefficients in bar plots
#'
#' Plots correlation coefficients in bar plot from given
#'     dataframe, as returned by computeCellLineCorrelation
#'     or computeDrugCorrelation.
#'
#' @param correlations A dataframe that contains correlation
#'     coefficients.
#' @param sensMeasure A character vector specifying which sensitivity
#'     measure to plot. Must be one of the measures in correlations, as
#'     returned by names(correlations).
#' @param coefficient A character vector specifying which
#'     correlation coefficient to plot. Must be one of the coefficients
#'     in correlations, as returned by names(correlations$sensMeasure).
#' @param title A character vector specifying the desired plot title.
#'
#' @return A bar plot with coefficient values on the y-axis, and either
#'     drugs or cell lines on the x-axis depending on the correlation
#'     dataframe inputted.
#'
#' @examples
#' # Intersect PharmacoSets of interest based on common cell lines
#' CTRP <- PharmacoGx::downloadPSet("CTRPv2_2015")
#' GRAY <- PharmacoGx::downloadPSet("GRAY_2013")
#' intersected <- PharmacoGx::intersectPSet(c(CTRP, GRAY),
#'     intersectOn = c("drugs", "cell.lines"))
#' correlations <- computeCorrelation(pSet = intersected,
#'     coefs = "pearson",
#'     sensMeasures = "aac_recomputed",
#'     pval = TRUE)
#' plotCorrelations(correlations = correlations,
#'     sensMeasure = "aac_recomputed_corrs",
#'     coefficient = "pearson",
#'     title = "Pearson Correlations of Recomputed AAC Values")
#'
#' @importFrom grDevices rainbow
#' @importFrom graphics par text
#' @importFrom utils combn
#' @export
#'

plotCorrelations <- function(correlations,
                             sensMeasure,
                             coefficient,
                             title) {
  # Performing Checks
  if (is.character(coefficient) == FALSE) {
    stop("coefficient must of class character, and must be one of the
         coefficients included in the correlations dataframe.")
  } else if (is.character(coefficient) == TRUE) {
    if (is.character(sensMeasure) == TRUE) {
      if (is.list(correlations) == TRUE) {
        if (!sensMeasure %in% names(correlations)) {
          stop("sensMeasure must be one of the objects in correlations.")
        } else if (is.list(correlations[[sensMeasure]]) == FALSE) {
          stop("correlations[sensMeasure] must be of class data.frame.")
        }
        if (!coefficient %in% names(correlations[[sensMeasure]])) {
          stop("coefficient must be one of the columns in
               correlations$sensMeasure.")
        }
      } else if (is.list(correlations) == FALSE) {
        stop("correlations must be of class list.")
      }
    } else if (is.character(sensMeasure) == FALSE) {
      stop("sensMeasure must be of class character.")
    }
  }

  rr <- as.matrix(t(cbind(correlations[[sensMeasure]][[coefficient]])))
  rr[!is.na(rr) & rr < 0] <- 0 # set negative values to 0
  names(rr) <- rownames(correlations[[sensMeasure]])
  par(mar = c(8, 5, 5, 5))
  rb <- graphics::barplot(rr,
                          beside = TRUE,
                          space = c(0.1, 2),
                          col=rep(rainbow(length(rr), v=0.9),
                                  each=2),
                          ylab= paste(coefficient, "correlation coefficient", sep = ' '),
                          ylim = c(0, 1),
                          density=c(100) ,
                          angle=c(0),
                          main = title,
                          xaxt = "n")
  text(x=apply(rb, 2, mean) + 3,
       y=par("usr")[3] - (par("usr")[4] * 0.05) + 0.03, pos=2,
       labels=toupper(names(rr)),
       srt=45,
       xpd=NA,
       font=1,
       cex = 0.7) # x-axis label size
}
#[END]
