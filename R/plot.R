#' Plot coefficients stored in a regsurv object
#'
#' @param x regsurv object
#' @param incl.baseline whether to include the parameters that relate to the baseline hazard or baseline cumulative hazard
#' @param scaled.betas plots coefficients for scaled columns of the design matrix if TRUE
#' @param ... additional arguments that are passed along to matplot()
#'
#' @return Nothing. Side-effect: plot.
#' @export
#' @method plot regsurv
plot.regsurv <- function(x, incl.baseline=FALSE, scaled.betas=FALSE, ...){

  if(!"regsurv" %in% class(x)){
    stop("plot.regsurv only takes objects of class regsurv as a first argument")
  }

  if(incl.baseline){
    cols <- unlist(x$which.param)
    lty <- rep(c(3,1,2), times=sapply(x$which.param, length))
  } else {
    cols <- unlist(x$which.param[2:3])
    lty <- rep(1:2, times=sapply(x$which.param[2:3], length))
  }

  if(!scaled.betas){
    y <- t(x$betahat[cols, ])
  } else {
    y <- t(x$betahat.scaled[cols, ])
  }
  defaults <- list(x=log(x$lambda.grid),
                   y=y,
                   type="l", lty=lty,
                   ylab="Coefficientcs",
                   xlab="log(lambda)",
                   main="Regularization path")

  arguments <- Reduce(utils::modifyList, list(defaults, list(...)))

  do.call(graphics::matplot, arguments)
}

#' Plot cv.regsurv object
#'
#' @param x cv.regsurv object
#' @param ... .. additional arguments that are passed along to plot()
#'
#' @return Nothing. Side-effect: plot.
#' @export
#' @method plot cv.regsurv
plot.cv.regsurv <- function(x, ...){

  if(!"cv.regsurv" %in% class(x)){
    stop("plot.cv.regsurv only takes objects of class cv.regsurv as a first argument")
  }

  y <- apply(x$oosll, 1, mean)
  cvlo <- -2*(y + x$cvse)
  cvup <- -2*(y - x$cvse)
  y <- -2*y

  defaults <- list(x=log(x$lambda.grid),
                   y=y,
                   pch=19,
                   col="red",
                   ylim=range(c(cvup, cvlo), na.rm=TRUE),
                   ylab="Deviance",
                   xlab="log(lambda)",
                   main="")

  arguments <- Reduce(utils::modifyList, list(defaults, list(...)))

  do.call(graphics::plot, arguments)
  graphics::arrows(x0=log(x$lambda.grid), y0=cvlo, y1=cvup, col="grey", code=3, angle=90, length = .04)
  lambda.min.index <- which(x$cvm == max(x$cvm))
  graphics::abline(v=log(x$lambda.grid[lambda.min.index]), lty=3, lwd=1.2)
}




