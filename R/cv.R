cv.regsurv <- function(object, prep, nfolds=10, plot=FALSE, print=FALSE){

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

  nrows <- nrow(prep$sbt$d)
  sampled.rows <- sample(1:nrows, nrows, replace = FALSE)
  cv.index <- oosll <- mr <- list()
  for(i in 1:nfolds){
    cv.index[[i]] <- sampled.rows[1:nrows %% nfolds + 1 == i]
  }

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

  i=1
  for(i in 1:nfolds){

    prep.cv <- survprep(tte=tte[-cv.index[[i]]],
                        delta=prep$delta[-cv.index[[i]]],
                        X=prep$X[-cv.index[[i]], ],
                        model.scale=prep$model.scale,
                        time.scale=prep$time.scale,
                        spline.type=prep$spline.type,
                        ntimebasis=ntimebasis,
                        tv=prep$tv,
                        nitimebasis=nitimebasis,
                        qpoints=qpoints)

    cv <- regsurv(prep=prep.cv, penpars=mod$penpars, l1l2=mod$l1l2, lambda.grid=mod$lambda.grid)

    check <- survprep(tte=tte[cv.index[[i]]],
                      delta=prep$delta[cv.index[[i]]],
                      X=prep$X[cv.index[[i]], ],
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
                      with(check, loss.hazard(beta=cv$betahat.scaled[ ,x], X=sbt$d, delta=delta, z=z, wl=wl))})
    }

    if(prep$model.scale == "logHazard"){
      if(prep$time.scale == "time"){
        oosll[[i]] <- sapply(1:length(cv$lambda.grid), function(x){
          deltaloghazard <- with(check, sbt$d[delta == 1, ] %*% cv$betahat.scaled[ ,x] + tte[delta == 1] +  # NB the offset doesn't cancel in case of linear time
                                   log(dsbt$d[delta == 1, ] %*% cv$betahat.scaled[ ,x][c(prep.cv$which.param[[1]], prep.cv$which.param[[3]])] + 1))
          cumhazard <- with(check, exp(sbt$d %*% cv$betahat.scaled[ ,x] + tte))
          sum(deltaloghazard) - sum(cumhazard)})
      } else {
        oosll[[i]] <- sapply(1:length(cv$lambda.grid), function(x){
          deltaloghazard <- with(check, sbt$d[delta == 1, ] %*% cv$betahat.scaled[ ,x] +
                                   log(dsbt$d[delta == 1, ] %*% cv$betahat.scaled[ ,x][c(prep.cv$which.param[[1]], prep.cv$which.param[[3]])] + 1))
          cumhazard <- with(check, exp(sbt$d %*% cv$betahat.scaled[ ,x] + tte))
          sum(deltaloghazard) - sum(cumhazard)})
      }
    }

    mrj <- list()
    for(j in 1:length(cv$lambda.grid)){
      if(prep$time.scale == "logtime"){
        preptte <- exp(check$tte)
      } else {
        preptte <- check$tte
      }
      pred <- predict.regsurv(cv, prep.cv, lambda.index=j, newdata=as.matrix(cbind(preptte, check$X)), type=c("cumhazard"))
      # mrj[[j]] <- mean((pred - check$delta)^2)
      mrj[[j]] <- mean(abs(pred - check$delta))
    }
    mr[[i]] <- unlist(mrj)

    if(print){print(paste("fold",i,"finished"))}
  }

  oosll <- matrix(unlist(oosll), nrow=length(mod$lambda.grid), ncol=nfolds)
  cvm <- apply(oosll, 1, mean)
  cvsd <- apply(oosll, 1, stats::sd)
  cvup <- cvm + cvsd
  cvlo <- cvm - cvsd

  mr <- matrix(unlist(mr), nrow=length(mod$lambda.grid), ncol=nfolds)

  if(plot){
    plot(apply(oosll, 1, mean) ~ log(mod$lambda.grid), pch=19, col="red", ylim=range(c(cvup, cvlo)),
         main="", xlab="log(lambda)", ylab="Deviance")
    graphics::arrows(x0=log(mod$lambda.grid), y0=cvlo, y1=cvup, col="grey", code=3, angle=90, length = .04)
    lambda.min.index <- which(cvm == max(cvm))
    graphics::abline(v=log(mod$lambda.grid[lambda.min.index]), lty=3)
  }

  return(
    structure(
      list(oosll=oosll,
         cvm=cvm,
         cvsd=cvsd,
         mr=mr,
         lambda.grid=mod$lambda.grid,
         lambda.min=mod$lambda.grid[which(cvm == max(cvm))],
         lambda.min.index=which(cvm == max(cvm))),
      class=c("cv.regsurv")))
}

