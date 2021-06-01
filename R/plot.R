#' Plot coefficients stored in a regsurv object
#'
#' @param x regsurv object
#' @param incl.baseline whether to include the parameters that relate to the baseline hazard or baseline cumulative hazard
#' @param ... additional arguments that are passed along to matplot()
#'
#' @return Nothing. Side-effect: plot.
#' @export
plot.regsurv <- function(x, incl.baseline=FALSE, ...){

  if(class(x) != "regsurv"){
    stop("plot.regsurv only takes objects of class regsurv as a first argument")
  }

  if(incl.baseline){
    cols <- unlist(x$which.param)
    lty <- rep(c(3,1,2), times=sapply(x$which.param, length))
  } else {
    cols <- unlist(x$which.param[2:3])
    lty <- rep(1:2, times=sapply(x$which.param[2:3], length))
  }

  defaults <- list(x=log(x$lambda.grid),
                   y=t(x$betahat[cols, ]),
                   type="l", lty=lty,
                   ylab="Coefficientcs",
                   xlab="log(lambda)",
                   main="Regularization path")

  arguments <- Reduce(utils::modifyList, list(defaults, list(...)))

  do.call(graphics::matplot, arguments)
}




