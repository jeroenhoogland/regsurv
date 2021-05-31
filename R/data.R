#' Simulated data for illustrative purposes
#'
#' Data were simulated using the simsurv package. The data generating mechanism is Weibull with
#' administrative censoring at time=10. The first 4 covariates have non-zero effects and the
#' remainder are noise variables. The first 2 covariate effects are time-varying. Details on the
#' data generating mechanism can be found in the source files (data-raw folder).
#'
#' @format A data frame with 250 rows and 12 columns
#' \describe{
#'   \item{id}{identifier}
#'   \item{eventtime}{eventtime}
#'   \item{status}{status indicator (1=event, 0=censored)}
#'   \item{x1-x9}{covariates 1-9}
#' }
"simdata"
