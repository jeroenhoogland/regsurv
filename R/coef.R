#' Return model coefficients for a regsurv model
#'
#' @param object regsurv object
#' @param s lambda value
#' @param ... other arguments
#'
#' @return coefficients for a fitted regsurv model at all lambda's (s=NULL) or a specific value of lambda (s)
#' @export
#' @method coef regsurv
coef.regsurv <- function(object, s=NULL, ...){

  if(!"regsurv" %in% class(object)){
    stop("coef.regsurv only takes objects of class regsurv as a first argument")
  }

  if(!is.null(s)){
    return(object$betahat[ ,which(object$lambda.grid == s)])
  } else {
    return(object$betahat)
  }
}
