#' Title
#'
#' @param object an object of class regsurv
#' @param prep an object of class survprep
#' @param nfolds number of cross-validation folds
#' @param plot plots deviance versus log(lambda) if TRUE
#' @param force.nnhazards the force.nnhazards status as used for regsurv() is always maintained during cross-validation, but the
#'   force.nnhazards parameter in cv.regsurv() allows to force these constraints on the whole sample (TRUE) (as opposed to just
#'   the in-fold cases (when FALSE)).
#' @param print prints progress if TRUE
#' @param maxit (integer) maxit the maximum number of iterations for the (ecos) solver, default 100L
#' @param feastol the tolerance on the primal and dual residual, default 1e-8
#' @param ... other parameters (only used by internal functions)
#'
#' @return
#'  \item{oosll}{out-of-sample log-likelihood}
#'  \item{cvm}{mean oosll per lambda}
#'  \item{cvse}{se of oosll per lambda}
#'  \item{msdr}{mean squared deviance residuals per lambda}
#'  \item{lambda.grid}{sed grid of lambda values (taken to be the same as those used to fit to regsurv model)}
#'  \item{lambda.min}{lambda value minimizing the (out-of-sample) deviance}
#'  \item{lambda.min.index}{index for the value of lambda minimizing the (out-of-sample) deviance}

#' @export
cv.regsurv <- function(object, prep, nfolds=10, plot=FALSE, force.nnhazards=TRUE, print=FALSE, maxit=100L, feastol=1e-08, ...){

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

  nrows <- nrow(prep$mm.scaled$d)
  sampled.rows <- sample(1:nrows, nrows, replace = FALSE)
  cv.index <- oosll <- deviance.res <- list()
  for(i in 1:nfolds){
    cv.index[[i]] <- sampled.rows[1:nrows %% nfolds + 1 == i]
  }

  # if(min(sapply(cv.index, function(x) sum(prep$delta[x]))) < 5){
  #   warning("some folds contain <5 events")
  # }

  if(prep$time.scale == "logtime"){
    tte <- exp(prep$tte)
  } else {
    tte <- prep$tte
  }

  if(prep$model.scale == "loghazard"){
    qpoints <- prep$qpoints
  } else {
    qpoints <- NULL
  }

  ntimebasis <- length(prep$which.param[[1]]) - 1

  if(is.null(prep$which.param[[3]])){
    nitimebasis <- NULL
  } else {
    nitimebasis <- length(prep$which.param[[3]]) / length(prep$tv)
  }

  check.negative.hazards.global <- list()
  i=1
  for(i in 1:nfolds){

    if(ncol(as.matrix(prep$X)) == 1){
      Xin <- as.matrix(prep$X[-cv.index[[i]]])
      colnames(Xin) <- colnames(prep$X)
    } else {
      Xin <- as.matrix(prep$X[-cv.index[[i]], ])
      colnames(Xin) <- colnames(prep$X)
    }

    # prep for in-sample cases
    prep.cv <- survprep(tte=tte[-cv.index[[i]]],
                        delta=prep$delta[-cv.index[[i]]],
                        X=Xin,
                        model.scale=prep$model.scale,
                        time.scale=prep$time.scale,
                        spline.type=prep$spline.type,
                        ntimebasis=ntimebasis,
                        tv=prep$tv,
                        nitimebasis=nitimebasis,
                        qpoints=qpoints)

    if(prep$model.scale == "logHazard" & force.nnhazards){
      prep.constrainst <- survprep(tte=tte,
                                   delta=prep$delta,
                                   X=as.matrix(prep$X),
                                   model.scale=prep$model.scale,
                                   time.scale=prep$time.scale,
                                   spline.type=prep$spline.type,
                                   ntimebasis=ntimebasis,
                                   time.knots=prep.cv$knots,
                                   tv=prep$tv,
                                   nitimebasis=nitimebasis,
                                   itime.knots=prep.cv$iknots,
                                   qpoints=qpoints,
                                   scales=prep.cv$scales,
                                   shifts=prep.cv$shifts)
      broad.Xd <- as.matrix(prep.constrainst$dmm.scaled$d)
      cv <- regsurv(prep=prep.cv, penpars=mod$penpars, l1l2=mod$l1l2, groups=mod$groups, lambda.grid=mod$lambda.grid,
                    force.nnhazards=mod$force.nnhazards, maxit=maxit, cv.constraint=broad.Xd, feastol=feastol)
    } else {
      cv <- regsurv(prep=prep.cv, penpars=mod$penpars, l1l2=mod$l1l2, groups=mod$groups, lambda.grid=mod$lambda.grid,
                    force.nnhazards=mod$force.nnhazards, maxit=maxit, feastol=feastol)
    }

    # prep for out-of-sample cases
    # takes care to use the estimated knots, iknots, scales, and shifts as in prep.cv
    if(ncol(as.matrix(prep$X)) == 1){
      Xout <- as.matrix(prep$X[cv.index[[i]]])
      colnames(Xout) <- colnames(prep$X)
    } else {
      Xout <- as.matrix(prep$X[cv.index[[i]], ])
      colnames(Xout) <- colnames(prep$X)
    }
    check <- survprep(tte=tte[cv.index[[i]]],
                      delta=prep$delta[cv.index[[i]]],
                      X=Xout,
                      model.scale=prep$model.scale,
                      time.scale=prep$time.scale,
                      spline.type=prep$spline.type,
                      ntimebasis=ntimebasis,
                      time.knots=prep.cv$knots,
                      tv=prep$tv,
                      nitimebasis=nitimebasis,
                      itime.knots=prep.cv$iknots,
                      qpoints=qpoints,
                      scales=prep.cv$scales,
                      shifts=prep.cv$shifts)

    if(prep$model.scale == "loghazard"){
      oosll[[i]] <- sapply(1:length(cv$lambda.grid), function(x){
                      with(check, loss.hazard(beta=cv$betahat.scaled[ ,x], X=mm.scaled$d, delta=delta, z=z, wl=wl))})
    }

    if(prep$model.scale == "logHazard"){
      ag.index <- c(prep.cv$which.param[[1]], prep.cv$which.param[[3]])
      if(prep$time.scale == "time"){
        check.negative.hazards.global[[i]] <- any(with(check, dmm.scaled$d[delta == 1, ] %*% cv$betahat.scaled[ag.index,] + 1) < 0)
        oosll[[i]] <- sapply(1:length(cv$lambda.grid), function(x){
          check.negative.hazards <- any(with(check, dmm.scaled$d[delta == 1, ] %*% cv$betahat.scaled[ ,x][ag.index] + 1) < 0)
          if(!is.na(check.negative.hazards) & check.negative.hazards){
            suppressWarnings(deltaloghazard <- with(check, mm.scaled$d[delta == 1, ] %*% cv$betahat.scaled[ ,x] + tte[delta == 1] +
                                       log(dmm.scaled$d[delta == 1, ] %*% cv$betahat.scaled[ ,x][ag.index] + 1)))
            deltaloghazard <- ifelse(is.na(deltaloghazard), min(deltaloghazard, na.rm=TRUE), deltaloghazard)
          } else {
            deltaloghazard <- with(check, mm.scaled$d[delta == 1, ] %*% cv$betahat.scaled[ ,x] + tte[delta == 1] +  # NB the offset doesn't cancel in case of linear time
                                     log(dmm.scaled$d[delta == 1, ] %*% cv$betahat.scaled[ ,x][ag.index] + 1))
          }
          cumhazard <- with(check, exp(mm.scaled$d %*% cv$betahat.scaled[ ,x] + tte))
          sum(deltaloghazard) - sum(cumhazard)})
      } else {
        check.negative.hazards.global[[i]] <- any(with(check, dmm.scaled$d[delta == 1, ] %*% cv$betahat.scaled[ag.index,] + 1) < 0)
        oosll[[i]] <- sapply(1:length(cv$lambda.grid), function(x){
          check.negative.hazards <- any(with(check, dmm.scaled$d[delta == 1, ] %*% cv$betahat.scaled[ ,x][ag.index] + 1) < 0)
          if(!is.na(check.negative.hazards) & check.negative.hazards){
            suppressWarnings(deltaloghazard <- with(check, mm.scaled$d[delta == 1, ] %*% cv$betahat.scaled[ ,x] +
                                   log(dmm.scaled$d[delta == 1, ] %*% cv$betahat.scaled[ ,x][ag.index] + 1)))
            deltaloghazard <- ifelse(is.na(deltaloghazard), min(deltaloghazard, na.rm=TRUE), deltaloghazard)
          } else {
            deltaloghazard <- with(check, mm.scaled$d[delta == 1, ] %*% cv$betahat.scaled[ ,x] +
                                     log(dmm.scaled$d[delta == 1, ] %*% cv$betahat.scaled[ ,x][ag.index] + 1))
          }
          cumhazard <- with(check, exp(mm.scaled$d %*% cv$betahat.scaled[ ,x] + tte))
          sum(deltaloghazard) - sum(cumhazard)})
      }
    }

    # deviance residuals
    deviance.res.j <- list()
    j=1
    for(j in 1:length(cv$lambda.grid)){
      if(prep$time.scale == "logtime"){
        checktte <- exp(check$tte)
      } else {
        checktte <- check$tte
      }
      pred <- predict.regsurv(cv, prep.cv, lambda.index=j, newdata=as.matrix(cbind(checktte, check$X)), type=c("cumhazard"))

      # martingale residuals m
      m <- check$delta - pred

      # deviance residuals d
      d <- sign(m) * sqrt(-2 * (m + check$delta * log(check$delta - m)))

      # squared deviance residuals
      deviance.res.j[[j]] <- mean(d^2)

    }

    deviance.res[[i]] <- unlist(deviance.res.j)

    if(print){print(paste("fold",i,"finished"))}
  }

  oosll <- matrix(unlist(oosll), nrow=length(mod$lambda.grid), ncol=nfolds)
  cvm <- apply(oosll, 1, mean)
  cvse <- apply(oosll, 1, stats::sd) / sqrt(nfolds)
  cvm <- ifelse(abs(cvm) == Inf, NA, cvm)
  cvse <- ifelse(abs(cvse) == Inf, NA, cvse)
  cvup <- cvm + cvse
  cvlo <- cvm - cvse

  msdr <- matrix(unlist(deviance.res), nrow=length(mod$lambda.grid), ncol=nfolds)

  if(plot){
    y <- cvm
    ylo <- -2*(y + cvse)
    yup <- -2*(y - cvse)
    y <- -2*y

    plot(y ~ log(mod$lambda.grid), pch=19, col="red", ylim=range(c(yup, ylo), na.rm=TRUE),
         main="", xlab="log(lambda)", ylab="Deviance")
    suppressWarnings(graphics::arrows(x0=log(mod$lambda.grid)[!is.na(ylo) & !is.na(yup)],
                     y0=ylo[!is.na(ylo) & !is.na(yup)],
                     y1=yup[!is.na(ylo) & !is.na(yup)],
                     col="grey", code=3, angle=90, length = .04))
    lambda.min.index <- which(cvm == max(cvm))
    graphics::abline(v=log(mod$lambda.grid[lambda.min.index]), lty=3)
  }

  if(any(unlist(check.negative.hazards.global)[!is.na(unlist(check.negative.hazards.global))])){
    warning("Model resulted in negative out-of-sample hazards which were set to the lowest valid value in the sample.
            Use force.nnhazards=TRUE to add the necesssary constraints")
  }

  return(
    structure(
      list(oosll=oosll,
         cvm=cvm,
         cvse=cvse,
         msdr=msdr,
         lambda.grid=mod$lambda.grid,
         lambda.min=mod$lambda.grid[which(cvm == max(cvm, na.rm = TRUE))],
         lambda.min.index=which(cvm == max(cvm, na.rm = TRUE)),
         cv.index=cv.index),
      class=c("cv.regsurv")))
}

