#' Nurses hom data
#'
#' This simulated data set consists of 1000 data points representing nurses from
#' 25 hospitals. The outcome of interest is job-related stress among nurses
#' (standardized) and the covariates are age (standardized), experience
#' (standardized), sex (0 = male, 1 = female), the type of ward in which the
#' nurse works (0 = general care, 1 = special care), and hospital (1, 2, ...,
#' 25).
#'
#' @docType data
#'
#' @name nurses_hom
#'
#' @format ## `nurses_hom`
#' \describe{
#'   \item{data_ipd}{Individual participant data. A data frame with 1000 rows and 6 variables:
#'     \describe{
#'       \item{gender}{factor: 0 = male, 1 = female}
#'       \item{age}{dbl: Standardized age}
#'       \item{experience}{dbl: Standardized experience}
#'       \item{wardtype}{factor: 0 = general care, 1 = special care}
#'       \item{stress}{dbl: Standardized stress levels}
#'       \item{hospital}{factor: Hospital ID, 1 to 25}
#'     }
#'   }
#'   \item{summary_stats_hom}{List of summary statistics for each hospital: (X'X, X'y, y'y, n).
#'              Use this for models with a global intercept.
#'   }
#'   \item{summary_stats_het}{List of summary statistics for each hospital: (X'X, X'y, y'y, n).
#'              Use this for models with a site-specific intercepts.
#'   }
#' }
"nurses_hom"