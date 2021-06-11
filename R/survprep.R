#' Prepare data for regsurv()
#'
#' @param tte numeric vector containing time-to-event data
#' @param delta numeric vector containing event indicator data (1=event, 0=censored)
#' @param X design matrix containing all necessary covariate data
#' @param model.scale choose between log hazard and log cumulative hazard
#' @param time.scale choose between log time and linear time for the time scale of the model
#' @param spline.type type of spline used to model time dependencies (rcs=restricted cubic splines,
#'     ns=natural splines (as implemented in the splines2 package))
#' @param ntimebasis integer number of basis columns for the spline representation of the log baseline hazard
#'    or log cumulative baseline hazard.
#' @param time.knots numeric vector containing baseline spline knot locations.
#'     Defaults to equally spaced knots (if NULL).
#' @param tv numeric vector of column indices for the columns in X with time-varying (i.e. non-proportional) effects.
#' @param nitimebasis number of basis columns for the spline representation of non-proportional
#'     covariate effects.
#' @param itime.knots numeric vector containing spline knot locations for time-varying covariate
#'     effects. Defaults to equally spaced knots (if NULL).
#' @param qpoints number of quadrature points for numerical integration over time (only applicable
#'     for models on the log hazard scale)
#' @param scales only used for internal functions
#' @param shifts only used for internal functions
#'
#' @return an object of class survprep
#' @export
#'
#' @examples
#' prep <- survprep(tte=simdata$eventtime,
#'                  delta=simdata$status,
#'                  X=as.matrix(simdata[ ,grep("x", names(simdata))]),
#'                  model.scale="loghazard",
#'                  time.scale="time",
#'                  spline.type="rcs",
#'                  ntimebasis=4,
#'                  qpoints=9)
survprep <- function(tte, delta, X,
                     model.scale=c("loghazard", "logHazard"),
                     time.scale=c("logtime", "time"),
                     spline.type=c("rcs","ns"),
                     ntimebasis, time.knots=NULL,
                     tv=NULL,
                     nitimebasis=NULL, itime.knots=NULL,
                     qpoints=9,
                     scales=NULL,
                     shifts=NULL){

  if(spline.type=="ns"){
    if (!requireNamespace("splines", quietly = TRUE)) {
      stop("Package \"splines2\" needed when requesting natural splines as the splines method. Please install it.",
           call. = FALSE)
    }
  }

  if(!is.numeric(tte)){
    stop("tte should be a numeric vector")
  }

  if(!is.numeric(delta)){
    stop("delta should be a numeric vector")
  }

  if(!is.matrix(X)){
    stop("X should be a matrix")
  }

  if(!is.numeric(ntimebasis)){
    stop("ntimebasis should be a numeric vector")
  }

  if(!is.null(tv) & !is.numeric(tv)){
    stop("tv should be a numeric vector")
  } else if (!all(tv %in% 1:ncol(X))){
    stop("not all of the provided tv indices are in the dimension of X")
  }

  if(!is.null(nitimebasis) & !is.numeric(nitimebasis)){
    stop("nitimebasis should be a numeric vector")
  } else if (!is.null(nitimebasis) & is.null(tv)){
    stop("nitimebasis specified in absence of any no time dependent variables (tv)")
  }

  if(time.scale == "logtime"){tte <- log(tte)}

  if(ntimebasis==0){time.type <- "constant"; knots <- NULL}
  if(ntimebasis==1){time.type <- "linear"; knots <- NULL}
  if(ntimebasis>1){time.type <- "spline"}

  if(is.null(nitimebasis)){
    itime.type <- "none"; iknots <- NULL
  } else {
    if(nitimebasis==0){itime.type <- "none"; iknots <- NULL}
    if(nitimebasis==1){itime.type <- "linear"; iknots <- NULL}
    if(nitimebasis>1){itime.type <- "spline"}
  }

  if(time.type != "spline" & itime.type != "spline"){spline.type <- NULL}

  # knots for the time representation of the baseline model
  if(time.type=="spline"){
    if(!is.null(time.knots)){
      if(length(time.knots) != (ntimebasis+1)){
        stop("Specified number of time.knots is not consistent with the
           requested number of basis columns (should follow length(time.knots) == (ntimebasis+1))")
      } else {
        knots <- time.knots
      }
    } else{
      knot.quantiles <- seq(0,1,1/(ntimebasis))
      knots <- stats::quantile(tte[delta==1], probs = knot.quantiles)
    }
  }

  ## Knots for the time-varying effects spline representation of time
  if(itime.type=="spline"){
    if(!is.null(itime.knots)){
      if(length(itime.knots) != (nitimebasis+1)){
        stop("Specified number of itime.knots is not consistent with the
           requested number of basis columns (should follow length(itime.knots) == (nitimebasis+1))")
      } else {
        iknots <- itime.knots
      }
    } else{
      iknot.quantiles <- seq(0,1,1/(ntimebasis))
      iknots <- stats::quantile(tte[delta==1], probs = iknot.quantiles)
    }
  }

  if(length(knots) != length(unique(knots))){
    knots <- unique(knots)
    warning("The number of time.knots was reduced to a set with unique elements")
  }

  if(length(iknots) != length(unique(iknots))){
    warning("The number of itime.knots was reduced to a set with unique elements")
  }

  survprep.id <- stats::runif(1)

  sbt <- sbi(t=tte, X=X, time.type=time.type, itime.type=itime.type, tv=tv,
             knots=knots, iknots=iknots, spline.type=spline.type)

  alpha <- rep(0, sbt$nbasis+1)
  beta <- rep(0, sbt$p)
  gamma <- rep(0, sbt$nibasis)

  if(length(gamma) == 0){
    which.param <- list(1:length(alpha),
                        (length(alpha)+1):(length(alpha)+length(beta)),
                        NULL)
  } else {
    which.param <- list(1:length(alpha),
                        (length(alpha)+1):(length(alpha)+length(beta)),
                        (length(alpha)+length(beta)+1):(length(alpha)+length(beta)+length(gamma)))
  }

  parameters <- c(alpha,beta,gamma)

  if(is.null(scales)){
    scales <- attr(scale(sbt$d[ ,-1]), "scaled:scale")
  }
  if(is.null(shifts)){
    shifts <- attr(scale(sbt$d[ ,-1]), "scaled:center")
  }

  if(model.scale == "loghazard"){
    # Gauss-Legendre weights and evaluation points for the log hazard
    rule <- legendre.quadrature.rules(qpoints)[[qpoints]]

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

    sbt.scaled <- sbt
    sbt.scaled$d[ ,-1] <- t((t(sbt.scaled$d[ ,-1]) - shifts) / scales)

    l <- sapply(glsbi, function(x) x$lambda)
    wl <- rep(rule$w, length(delta)) * rep(l, each=length(rule$w))
    z <- lapply(glsbi, function(x) x$z)
    z <- Reduce(rbind, z)

    return(
      structure(
        list("model.scale"=model.scale,
           "time.scale"=time.scale,
           "spline.type"=spline.type,
           "tte"=tte, # NB is log(tte) when time.scale equals "logtime"
           "delta"=delta,
           "sbt"=sbt.scaled,
           "scales"=scales,
           "shifts"=shifts,
           "w"=rule$w,
           "l"=l,
           "wl"=wl,
           "z"=z,
           "which.param"=which.param,
           "parameters"=parameters,
           "X"=X,
           "knots"=knots,
           "iknots"=iknots,
           "qpoints"=qpoints,
           "time.type"=time.type,
           "itime.type"=itime.type,
           "tv"=tv,
           "survprep.id"=survprep.id),
        class="survprep"))
  }

  if(model.scale == "logHazard"){

    sbt.scaled <- sbt
    sbt.scaled$d[ ,-1] <- t((t(sbt.scaled$d[ ,-1]) - shifts) / scales)

    dsbt <- dsbi(t=tte, X=X, time.type=time.type, itime.type=itime.type, tv=tv,
                 knots=knots, iknots=iknots, spline.type=spline.type)
    dsbt.scaled <- dsbt
    dsbt.scaled$d[ ,-1] <- t(t(dsbt.scaled$d[ ,-1]) / scales[grep("basis", names(scales))])

    return(
      structure(
        list("model.scale"=model.scale,
         "time.scale"=time.scale,
         "spline.type"=spline.type,
         "tte"=tte, # NB is log(tte) when time.scale equals "logtime"
         "delta"=delta,
         "sbt"=sbt.scaled,
         "dsbt"=dsbt.scaled,
         "scales"=scales,
         "shifts"=shifts,
         "which.param"=which.param,
         "parameters"=parameters,
         "X"=X,
         "knots"=knots,
         "iknots"=iknots,
         "time.type"=time.type,
         "itime.type"=itime.type,
         "tv"=tv,
         "survprep.id"=survprep.id),
        class="survprep"))
  }
}


sbi <- function(t, X, time.type, itime.type, tv=NULL, knots=NULL, iknots=NULL, spline.type=NULL){

  # basis matrix for baseline model
  if(time.type=="constant"){
    basis <- NULL
  }

  if(time.type=="linear"){
    basis <- t
  }

  if(time.type=="spline"){
    if(spline.type=="rcs"){
      basis <- rcs(t, knots)
    }
    if(spline.type=="ns"){
      outer.knots <- range(knots)
      knots <- knots[!knots %in% outer.knots]
      basis <- splines2::naturalSpline(t, knots=knots, Boundary.knots = outer.knots)
    }
  }

  # add basis matrix for time-varying effects
  if(itime.type=="none"){
    ibasis <- NULL
    d <- cbind(1, basis, X)
  }

  if(itime.type=="linear"){
    ibasis <- as.matrix(t)
    d <- cbind(1, basis, X, as.vector(ibasis) * X[ ,tv])
  }

  if(itime.type=="spline"){
    if(spline.type == "ns"){
      iknots <- iknots[!iknots %in% outer.knots]
      ibasis <- splines2::naturalSpline(t, knots=iknots, Boundary.knots = outer.knots)
      intbasis <- ibasis[ ,1] * X[ ,tv]
      for(i in 2:(ncol(ibasis))){
          intbasis <- cbind(intbasis, ibasis[ ,i] * X[ ,tv])
      }
      d <- cbind(1, basis, X, intbasis)
    }

    if(spline.type == "rcs"){
      ibasis <- rcs(t, iknots)
      intbasis <- ibasis[ ,1] * X[ ,tv]
      for(i in 2:(ncol(ibasis))){
          intbasis <- cbind(intbasis, ibasis[ ,i] * X[ ,tv])
      }
      d <- cbind(1, basis, X, intbasis)
    }
  }

  if(time.type=="constant"){
    if(itime.type=="none"){
      colnames(d) <- c("int", colnames(X))
    } else {
      colnames(d) <- c("int", colnames(X),
                       paste0("ibasis", 1:(ncol(ibasis)*ncol(as.matrix(X[ ,tv])))))
    }
  } else {
    if(itime.type=="none"){
      colnames(d) <- c("int", paste0("basis", 1:ncol(as.matrix(basis))), colnames(X))
    } else {
      colnames(d) <- c("int", paste0("basis", 1:ncol(as.matrix(basis))), colnames(X),
                       paste0("ibasis", 1:(ncol(ibasis)*ncol(as.matrix(X[ ,tv])))))
    }
  }

  return(list(d=d,
         nbasis=ifelse(!is.null(basis), ncol(as.matrix(basis)), 0),
         p=ncol(X),
         nibasis=ifelse(!is.null(ibasis), ncol(ibasis)*ncol(as.matrix(X[ ,tv])), 0)))
}

dsbi <- function(t, X, time.type, itime.type, tv, knots=NULL, iknots=NULL, spline.type=NULL){

  # derivatives of baseline model basis matrix
  if(time.type=="constant"){
    basis <- NULL
  }

  if (time.type=="linear"){
    basis <- matrix(1, nrow = length(t), ncol = 1)
  }

  if(time.type=="spline"){
    if(spline.type=="rcs"){
      basis <- drcs(t, knots)
    }
    if(spline.type=="ns"){
      outer.knots <- range(knots)
      knots <- knots[!knots %in% outer.knots]
      basis <- splines2::naturalSpline(t, knots=knots, Boundary.knots = outer.knots, derivs = 1)
    }
  }

  # derivatives of time-varying effects basis matrix
  if(itime.type=="none"){
    ibasis <- NULL
    d <- cbind(0, basis)
  }

  if(itime.type=="linear"){
    ibasis <- matrix(1, nrow=length(t), ncol=1)
    d <- cbind(0, basis, X[ ,tv])
  }

  if(itime.type=="spline"){
    if(spline.type == "ns"){
      outer.iknots <- range(iknots)
      iknots <- iknots[!iknots %in% outer.iknots]
      ibasis <- splines2::naturalSpline(t, knots=iknots, Boundary.knots = outer.iknots, derivs = 1)
      intbasis <- ibasis[ ,1] * X[ ,tv]
      for(i in 2:(ncol(ibasis))){
        intbasis <- cbind(intbasis, ibasis[ ,i] * X[ ,tv])
      }
      d <- cbind(0, basis, intbasis)
    }

    if(spline.type == "rcs"){
      ibasis <- drcs(t, iknots)
      intbasis <- ibasis[ ,1] * X[ ,tv]
      for(i in 2:(ncol(ibasis))){
        intbasis <- cbind(intbasis, ibasis[ ,i] * X[ ,tv])
      }
      d <- cbind(0, basis, intbasis)
    }
  }

  if(time.type=="constant"){
    if(itime.type=="none"){
      colnames(d) <- c("int")
    } else {
      colnames(d) <- c("int", paste0("ibasis", 1:(ncol(ibasis)*ncol(as.matrix(X[ ,tv])))))
    }
  } else {
    if(itime.type=="none"){
      colnames(d) <- c("int", paste0("basis", 1:ncol(as.matrix(basis))))
    } else {
      colnames(d) <- c("int", paste0("basis", 1:ncol(as.matrix(basis))),
                       paste0("ibasis", 1:(ncol(ibasis)*ncol(as.matrix(X[ ,tv])))))
    }
  }

  return(list(d=d,
              nbasis=ifelse(!is.null(basis), ncol(as.matrix(basis)), 0),
              nibasis=ifelse(!is.null(ibasis), ncol(ibasis)*ncol(as.matrix(X[ ,tv])), 0)))
}


