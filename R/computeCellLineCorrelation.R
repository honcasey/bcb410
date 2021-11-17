#' Compute correlations between drug sensitivity measures of cell lines
#'
#' Computes Pearson, Spearman, and/or Kendall Correlation Coefficients and
#'     associated p-values for specified drug sensitivity measures
#'     (such as IC50 of AUC) between common cell lines of two or more
#'     PharmacoSets.
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
#' @return Returns a list of dataframes, one dataframe for each sensitivity
#'     measure of interest, where rows correspond to cell lines and columns
#'     correspond to each type of correlation coefficient from coefs.
#'     Also includes p-values in columns is pval was set to TRUE.
#'
#' @examples
#' # Intersect PharmacoSets of interest based on common cell lines
#' CTRP <- PharmacoGx::downloadPSet("CTRPv2_2015")
#' GRAY <- PharmacoGx::downloadPSet("GRAY_2013")
#' intersected <- PharmacoGx::intersectPSet(c(CTRP, GRAY),
#'     intersectOn = c("drugs", "cell.lines"))
#' correlations <- computeCellLineCorrelation(pSet = intersected,
#'     coefs = "pearson",
#'     sensMeasures = "aac_recomputed",
#'     pval = TRUE)
#'
#'
#'
#' @author {Casey Hon, \email{casey.hon@mail.utoronto.ca}}
#'
#' @references
#' Smirnov P, Safikhani Z, El-Hachem N, Wang D, She A, Olsen C, Freeman M,
#'     Selby H, Gendoo D, Grossman P, Beck A, Aerts H, Lupien M, Haibe-Kains
#'     AG, (2016). “PharmacoGx: an R package for analysis of large
#'     pharmacogenomic datasets.” \emph{Bioinformatics (Oxford, England)}.
#'
#'
#' @import PharmacoGx
#' @import stats
#'

computeCellLineCorrelation <- function(pSet,
                                       coefs,
                                       sensMeasures,
                                       pval = TRUE) {

  # Performing checks
  if (is.list(pSet) == TRUE) {
    if (length(pSet) < 2) {
      stop("pSet list should have at least two PharmacoSets to
           compare correlations.")
    }
  } else if (is.list(pSet) != TRUE) {
    stop("pSet should be of class list, which contains a list of
         intersected PharmacoSets.")
  }

  if (is.character(coefs) == TRUE) {
    coefsUsed <- c("pearson", "spearman", "kendall")
    if (all((coefs %in% coefsUsed) == FALSE)) {
      stop("coefs should be of class character, specifying either:
           pearson, spearman, and/or kendall.")
    }
  } else if (is.character(coefs) != TRUE) {
      stop("coefs should be of class character, specifying either:
           pearson, spearman, and/or kendall.")
  }
  if (is.logical(pval) != TRUE) {
      stop("pval should be of class logical, specifying either TRUE or FALSE.")
  }

  # separate each pset after intersection
  pSetList <- list()
  # sensUsed <- list()
  for (set in pSet) {
    x <- set@annotation[["name"]]
    pSetList <- append(pSetList, x)
    assign(x, set, envir = globalenv())
  }

  sensUsed <- intersectSensMeasures(pSet)
  # put in package name here? cellLineConsistency::intersectSensMeasures()?

  if (is.character(sensMeasures) == TRUE) {
    # check if all measures inputted are in common in all psets of choice
    if (all((sensMeasures %in% sensUsed) == FALSE)) {
      stop("sensMeasures should be of class character, specifying the drug
         sensitivity measures. The list must be a subset of the measures in
           common to all PSets intersected.")
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

  # if there is more than two psets then add pset names for each coef and pair
  if (length(pSet) > 2) {
    # go through each pair in pSetList and paste new coef to list
    for (pair in colnames(pSetPairs)) {
      tempString <- paste(pSetPairs[1, pair], pSetPairs[2, pair], sep = '&')
      coefs <- c(coefs, paste(coefs, tempString, sep = '_'))
    }
  }

  # get sensitivity profiles for each sensitivity measure for each pset
  profList <- list() # list of each psets sensitivity profiles
  for (tempPSet in pSet) {
    setProfList <- list() # list sensitivity profiles of current pset
    for (sensName in sensMeasures) {
      pset_name <- tempPSet@annotation[["name"]]
      var_name <- paste(sensName, pset_name, sep = "_")
      prof <- PharmacoGx::summarizeSensitivityProfiles(tempPSet,
                                                       sensitivity.measure = sensName)
      assign(var_name, as.data.frame(prof))
      setProfList[sensName] <- list(get(var_name))
    }
    profList[pset_name] <- list(setProfList)
  }
  drugs <- rownames(pSet[[1]]@drug)
  cells <- rownames(pSet[[1]]@cell)

  # initialize data frames to be outputted
  cor_list <- list()
  for (sens in sensMeasures) {
    cors <- data.frame(matrix(ncol = length(coefs) * ncol(pSetPairs),
                             # columns are each type of correlation coefficient
                             nrow = length(cells)),
                             row.names = cells)
    colnames(cors) <- coefs
    var_name <- paste(sens, "corrs", sep = '_')
    assign(var_name, cors, envir = globalenv())
    cor_list[[var_name]] <- get(var_name) # TO-DO: fix warning here (cus they're dataframes)
  }

  # compute each correlation in coef between each pair of psets
  for (pair in colnames(pSetPairs)) { # each PSet
    set1 <- profList[[pSetPairs[1, pair][[1]]]]
    set2 <- profList[[pSetPairs[2, pair][[1]]]]
    for (sens in sensMeasures) { # each sensitivity measure (aac, ic50)
      for (cell.line in rownames(cors)) { # each cell line
        tryCatch(
          expr = {
      tofill <- get(paste(sens, "corrs", sep = '_'))
      set1sens <- as.numeric(unlist(set1[[sens]][cell.line]))
      set2sens <- as.numeric(unlist(set2[[sens]][cell.line]))
      # PEARSON CORRELATION
      if ("pearson" %in% coefs) {
        pearson.cor <- stats::cor.test(x = set1sens,
                                       y = set2sens,
                                       method = 'pearson',
                                       use = 'pairw',
                                       exact = FALSE)
        tofill[cell.line, "pearson"] <- pearson.cor$estimate
        if ("pearson p-value" %in% coefs) {
          tofill[cell.line, "pearson p-value"] <- pearson.cor$p.value
        }
      }
      # SPEARMAN CORRELATION
      if ("spearman" %in% coefs) {
        spearman.cor <- stats::cor.test(x = set1sens,
                                        y = set2sens,
                                        method = 'spearman',
                                        use = 'pairw',
                                        exact = FALSE)
        tofill[cell.line, "spearman"] <- spearman.cor$estimate
        if ("spearman p-value" %in% coefs) {
          tofill[cell.line, "spearman p-value"] <- spearman.cor$p.value
        }
      }
      # KENDALL CORRELATION
      if ("kendall" %in% coefs) {
        kendall.cor <- stats::cor.test(x = set1sens,
                                       y= set2sens,
                                        method = 'kendall',
                                        use = 'pairw',
                                       exact = FALSE)
        tofill[cell.line, "kendall"] <- kendall.cor$estimate
        if ("kendall p-value" %in% coefs) {
          tofill[cell.line, "kendall p-value"] <- kendall.cor$p.value
        }
      }
      curr_sens <- paste(sens, "corrs", sep = "_")
      assign(curr_sens, tofill)
      cor_list[[curr_sens]] <- get(curr_sens)
    },
      error = {function(e) {
        message(cell.line, " does not have enough finite observations to compute
                a correlation, so is left as NA.")
        }
      }
    )
      }
    }
  }
  return(cor_list)
}
