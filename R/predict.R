#' Predictions for a regsurv model
#'
#' @param object regsurv object
#' @param prep survprep object used to fit the regsurv model
#' @param lambda.index lambda index for the required predictions (e.g. optimal value based on cross-validation)
#' @param newdata if provided, this should be a matrix with untransformed time-to-event in the first column and the model.matrix for
#'   for main main predictor effects as would be provided to survprep()
#' @param type type of prediction ("cumhazard" or "surv")
#' @param ... as for predict()
#'
#' @return prediction of the requested type
#' @export
#' @method predict regsurv
predict.regsurv <- function(object, prep, lambda.index, newdata=NULL, type=c("cumhazard", "surv"), ...){

  if(class(object) != "regsurv"){
    stop("predict.regsurv only takes objects of class regsurv as a first argument")
  }

  if(class(prep) != "survprep"){
    stop("regsurv only takes objects of class survprep as a second argument")
  }

  if(prep$survprep.id != object$survprep.id){
    stop("regsurv object and survprep object do not match")
  }

  mod <- object
  spline.type <- prep$spline.type
  time.scale <- prep$time.scale
  knots <- prep$knots
  iknots <- prep$iknots
  time.type <- prep$time.type
  itime.type <- prep$itime.type
  tv <- prep$tv
  scales <- prep$scales
  shifts <- prep$shifts
  if(prep$model.scale == "loghazard"){
    qpoints <- prep$qpoints
    rule <- legendre.quadrature.rules(qpoints)[[qpoints]]
  }
  betahat <- as.matrix(mod$betahat.scaled)[ ,lambda.index]

  if(is.null(newdata)){
    X <- prep$mm$d[ ,prep$which.param[[2]], drop=FALSE]
    tte <- prep$tte
  } else {
    X <- newdata[ ,-1, drop=FALSE]
    if(prep$time.scale == "logtime"){
      tte <- log(newdata[ ,1])
    } else {
      tte <- newdata[ ,1]
    }
  }

  if(prep$model.scale == "loghazard"){
    if(time.scale == "time"){
      glsbi <- lapply(1:nrow(X), function(i)
      {
        lower = 0 # time zero; might alter in case of left-truncation
        lambda = (tte[i] - lower)/2 # the upper limit of integration is x
        mu = (lower + tte[i])/2 # equivalent to lambda for current application (with lower = 0)
        y = lambda * rule$x + mu
        Xi <- as.matrix(X[rep(i, length(y)), ])
        colnames(Xi) <- colnames(X)

        z <- sbi(t=y, X=Xi, time.type=time.type, itime.type=itime.type, tv=tv,
                 knots=knots, iknots=iknots, spline.type=spline.type)$d
        z[ ,-1] <- t((t(z[ ,-1]) - shifts) / scales)
        list(w=rule$w, lambda=lambda, z=z)
      })
    } else { # so for logtime
      glsbi <- lapply(1:nrow(X), function(i)
      {
        lower = 0 # time zero; might alter in case of left-truncation
        lambda = (exp(tte[i]) - lower)/2 # the upper limit of integration is x
        mu = (lower + exp(tte[i]))/2 # equivalent to lambda for current application (with lower = 0)
        y = lambda * rule$x + mu
        Xi <- as.matrix(X[rep(i, length(y)), ])
        colnames(Xi) <- colnames(X)

        z <- sbi(t=log(y), X=Xi, time.type=time.type, itime.type=itime.type, tv=tv,
                 knots=knots, iknots=iknots, spline.type=spline.type)$d
        z[ ,-1] <- t((t(z[ ,-1]) - shifts) / scales)
        list(w=rule$w, lambda=lambda, z=z)
      })
    }
    if(type=="cumhazard"){return(sapply(glsbi, integrate, param=betahat))}
    if(type=="surv"){return(exp(-sapply(glsbi, integrate, param=betahat)))}
  }

  if(prep$model.scale == "logHazard"){
    sbt <- sbi(t=tte, X=X, time.type=time.type, itime.type=itime.type, tv=tv,
               knots=prep$knots, iknots=prep$iknots, spline.type=spline.type)

    sbt$d[ ,-1] <- t((t(sbt$d[ ,-1]) - shifts) / scales)

    # if(prep$spline.type == "rcs"){
    #   if(type=="cumhazard"){return(exp(sbt$d %*% betahat))}
    #   if(type=="surv"){return(exp(-exp(sbt$d %*% betahat)))}
    # } else {
      if(type=="cumhazard"){return(exp(sbt$d %*% betahat + tte))}
      if(type=="surv"){return(exp(-exp(sbt$d %*% betahat + tte)))}
    # }
  }
}


integrate <- function(x, param){
  x$lambda * sum(x$w * exp(x$z %*% param))
}

