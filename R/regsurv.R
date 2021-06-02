#' Regularized parametric survival modeling
#'
#' @param prep an object of class survprep
#' @param penpars numeric vector indicating penalized (1) and unpenalized (0) parameters. The order should follow the order the survprep
#'   model matrix
#' @param l1l2 numeric vector indicating lasso (1 or TRUE) or ridge (0 or FALSE) penalty per parameter; order as for penpars
#' @param lambda.grid a grid of lambda values may be specified manually here. Default behavior is to fit all models from lambda
#'   equal to exp(-6) up to the moment where all lasso penalized parameters have disappeared and all the absolute
#'   value of all ridge penalized parameters is < 1e-2.
#' @param print if TRUE, prints progress (a line for each lambda for which the model was optimized)
#'
#' @return an object of class regsurv
#'  \item{optimal}{TRUE when the optimization converged to an optimal value for each lambda}
#'  \item{lambda.grid}{grid of lambda values}
#'  \item{obj.value}{objective function values per lambda}
#'  \item{betahat}{matrix of model coefficients with different parameters in rows and a column per
#'   lambda value}
#'  \item{num.iters}{number of iterations needed to reach the optimal solution (per lambda)}
#'  \item{solve.times}{solve times per lambda}
#'  \item{which.param}{includes 3 lists: the first with baseline model parameter indices, the
#'  second with main effect parameter indices, and the third with parameter indices for
#'  time-varying effect parameters}
#'
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
#' # e.g. penalize all parameters but the intercept
#' # note that prep$which.param includes 3 lists: the first with baseline model parameters, the
#' # second with main effect parameters, and the third with parameters for time-varying effects
#' penpars <- c(rep(TRUE, length(prep$which.param[[1]])),
#'              rep(TRUE, length(unlist(prep$which.param[2:3]))))
#' penpars[1] <- FALSE # do not penalize the intercept parameter
#'
#' # and use ridge for baseline parameters and lasso for all other parameters
#' l1l2 <- c(rep(0, length(unlist(prep$which.param[1]))),
#'           rep(1, length(unlist(prep$which.param[2:3]))))
#'
#' # fit model over the default lambda grid
#' mod <- regsurv(prep, penpars, l1l2, print=TRUE)
#' plot(mod)
regsurv <- function(prep, penpars, l1l2, lambda.grid=NULL, print=FALSE){

  if(class(prep) != "survprep"){
    stop("regsurv only takes objects of class survprep as a first argument")
  }

  X <- as.matrix(prep$sbt$d)
  p <- length(prep$parameters)
  beta <- CVXR::Variable(p)
  delta <- prep$delta
  tte <- prep$tte

  if(prep$model.scale == "loghazard"){

    z <- prep$z
    wl <-  prep$wl

    lambda <- pi
    obj <- loss.hazard(beta, X, delta, z, wl) -
      elastic_penalty(beta, penpars, lambda, l1l2)
    prob <- CVXR::Problem(CVXR::Maximize(obj))
    prob_data <- CVXR::get_problem_data(prob, solver="ECOS")

    # prob_data$data$G@x has all information that depends on lambda
    # stored in a dgCMatrix (https://slowkow.com/notes/sparse-matrix/#the-compressed-column-format-in-dgcmatrix)
    # the @x component contains all non-zero elements sorted by column

    ridge.index <- which(prob_data$data$G@x == -2 * sqrt(pi))
    lasso.index.plus <- which(prob_data$data$G@x == pi)
    lasso.index.min <- which(prob_data$data$G@x == -pi)

    sol <- list()
    if(is.null(lambda.grid)){
      lambda.grid <- exp(-6)
      cont <- TRUE
      i=1
      while(cont){
        prob_data$data$G@x[ridge.index] <- -2*sqrt(lambda.grid[i])
        prob_data$data$G@x[lasso.index.plus] <- lambda.grid[i]
        prob_data$data$G@x[lasso.index.min] <- -lambda.grid[i]
        solver_output <- ECOSolveR::ECOS_csolve(c = prob_data$data$c,
                                                G=  prob_data$data$G,
                                                h = prob_data$data$h,
                                                dims = list(l=as.numeric(prob_data$data$dims@nonpos),
                                                            q=as.numeric(prob_data$data$dims@soc),
                                                            e=as.numeric(prob_data$data$dims@exp)))
        sol[[i]] <- CVXR::unpack_results(prob, solver_output, prob_data$chain, prob_data$inverse_data)
        betahat <- sol[[i]]$getValue(beta)
        if(print){
          print(paste0("solved for lambda = log(", log(lambda.grid[i]), ")"))
        }
        cont <- !all(all(abs(betahat[penpars & l1l2]) < 1e-8), all(abs(betahat[penpars & !l1l2]) < 1e-2))
        if(cont){
          lambda.grid <- c(lambda.grid, exp(-6+i))
          i <- i+1
        }
      }
    } else {
      for(i in 1:length(lambda.grid)){
        prob_data$data$G@x[ridge.index] <- -2*sqrt(lambda.grid[i])
        prob_data$data$G@x[lasso.index.plus] <- lambda.grid[i]
        prob_data$data$G@x[lasso.index.min] <- -lambda.grid[i]
        solver_output <- ECOSolveR::ECOS_csolve(c = prob_data$data$c,
                                                G=  prob_data$data$G,
                                                h = prob_data$data$h,
                                                dims = list(l=as.numeric(prob_data$data$dims@nonpos),
                                                            q=as.numeric(prob_data$data$dims@soc),
                                                            e=as.numeric(prob_data$data$dims@exp)))
        sol[[i]] <- CVXR::unpack_results(prob, solver_output, prob_data$chain, prob_data$inverse_data)
        if(print){
          print(paste("solved for lambda", i, "out of", length(lambda.grid)))
        }
      }
    }
  }

  if(prep$model.scale == "logHazard"){
    ag.index <- c(prep$which.param[[1]], prep$which.param[[3]])
    Xd <- as.matrix(prep$dsbt$d)

    lambda <- pi
    if(prep$time.scale == "time"){
      obj <- loss.Hazard.time(beta, X, Xd, tte, delta, ag.index) -
        elastic_penalty(beta, penpars, lambda, l1l2)
    }
    if (prep$time.scale == "logtime"){
      obj <- loss.Hazard.logtime(beta, X, Xd, tte, delta, ag.index) -
        elastic_penalty(beta, penpars, lambda, l1l2)
    }

    prob <- CVXR::Problem(CVXR::Maximize(obj))
    prob_data <- CVXR::get_problem_data(prob, solver="ECOS")

    # prob_data$data$G@x has all information that depends on lambda
    # stored in a dgCMatrix (https://slowkow.com/notes/sparse-matrix/#the-compressed-column-format-in-dgcmatrix)
    # the @x component contains all non-zero elements sorted by column
    ridge.index <- which(prob_data$data$G@x == -2 * sqrt(pi))
    lasso.index.plus <- which(prob_data$data$G@x == pi)
    lasso.index.min <- which(prob_data$data$G@x == -pi)

    sol <- list()
    if(is.null(lambda.grid)){
      lambda.grid <- exp(-6)
      cont <- TRUE
      i=1
      while(cont){
        prob_data$data$G@x[ridge.index] <- -2*sqrt(lambda.grid[i])
        prob_data$data$G@x[lasso.index.plus] <- lambda.grid[i]
        prob_data$data$G@x[lasso.index.min] <- -lambda.grid[i]
        solver_output <- ECOSolveR::ECOS_csolve(c = prob_data$data$c,
                                                G=  prob_data$data$G,
                                                h = prob_data$data$h,
                                                dims = list(l=as.numeric(prob_data$data$dims@nonpos),
                                                            q=as.numeric(prob_data$data$dims@soc),
                                                            e=as.numeric(prob_data$data$dims@exp)))
        sol[[i]] <- CVXR::unpack_results(prob, solver_output, prob_data$chain, prob_data$inverse_data)
        betahat <- sol[[i]]$getValue(beta)
        if(print){
          print(paste0("solved for lambda = log(", log(lambda.grid[i]), ")"))
        }
        cont <- !all(all(abs(betahat[penpars & l1l2]) < 1e-8), all(abs(betahat[penpars & !l1l2]) < 1e-2))
        if(cont){
          lambda.grid <- c(lambda.grid, exp(-6+i))
          i <- i+1
        }
      }
    } else {
      for(i in 1:length(lambda.grid)){
        prob_data$data$G@x[ridge.index] <- -2*sqrt(lambda.grid[i])
        prob_data$data$G@x[lasso.index.plus] <- lambda.grid[i]
        prob_data$data$G@x[lasso.index.min] <- -lambda.grid[i]
        solver_output <- ECOSolveR::ECOS_csolve(c = prob_data$data$c,
                                                G=  prob_data$data$G,
                                                h = prob_data$data$h,
                                                dims = list(l=as.numeric(prob_data$data$dims@nonpos),
                                                            q=as.numeric(prob_data$data$dims@soc),
                                                            e=as.numeric(prob_data$data$dims@exp)))
        sol[[i]] <- CVXR::unpack_results(prob, solver_output, prob_data$chain, prob_data$inverse_data)
        if(print){
          print(paste("solved for lambda", i, "out of", length(lambda.grid)))
        }
      }
    }
  }

  status <- sapply(sol, function(x) x$status)
  optimal <- all(status == "optimal")
  if(!optimal){
    warning("solver did not find an optimal solution for all lambda's")
  }
  num.iters <- sapply(sol, function(x) x$num_iters)
  solve.times <- sapply(sol, function(x) x$solve_time)
  obj.value <- sapply(sol, function(x) x$value)
  betahat <- betahat.scaled <- sapply(sol, function(x) x$getValue(beta))

  if(prep$model.scale == "loghazard"){
    betahat.scaled <- betahat
    betahat[1, ] <- betahat.scaled[1, ] - as.vector(prep$shifts / prep$scales) %*% betahat.scaled[-1, ]
    betahat[-1, ] <- betahat.scaled[-1, ] / prep$scales
  }
  if(prep$model.scale == "logHazard"){
    betahat.scaled.excl.offset <- betahat
    betahat.scaled <- betahat
    betahat.scaled[2, ] <- betahat.scaled[2, ] + prep$scales["basis1"]
    betahat.scaled[1, ] <- betahat.scaled[1, ] + prep$shifts["basis1"]

    betahat <- betahat.scaled
    betahat[1, ] <- betahat.scaled[1, ] - as.vector(prep$shifts / prep$scales) %*% betahat.scaled[-1, ]
    betahat[-1, ] <- betahat.scaled[-1, ] / prep$scales
  }

  return(
    structure(
      list(optimal=optimal,
           lambda.grid=lambda.grid,
           obj.value=obj.value,
           betahat=betahat,
           betahat.scaled=betahat.scaled,
           num.iters=num.iters,
           solve.times=solve.times,
           which.param=prep$which.param),
      class="regsurv"))
}

loss.hazard <- function(beta, X, delta, z, wl){
  deltaloghazard <- X[delta==1, ] %*% beta
  cumhazard <- sum(exp(z %*% beta) * wl)
  sum(deltaloghazard) - sum(cumhazard)
}

loss.Hazard.time <- function(beta, X, Xd, tte, delta, ag.index){
  deltaloghazard <- X[delta == 1, ] %*% beta + tte[delta == 1] +  # NB the offset doesn't cancel in case of linear time
    CVXR::log1p(Xd[delta == 1, ] %*% beta[ag.index])
  cumhazard <- exp(X %*% beta + tte)
  sum(deltaloghazard) - sum(cumhazard)
}

loss.Hazard.logtime <- function(beta, X, Xd, tte, delta, ag.index){
  deltaloghazard <- X[delta == 1, ] %*% beta +
                          CVXR::log1p(Xd[delta == 1, ] %*% beta[ag.index])
  cumhazard <- exp(X %*% beta + tte)
  sum(deltaloghazard) - sum(cumhazard)
}

elastic_penalty <- function(beta, penpars, lambda, l1l2) {
  lambda <- penpars * lambda
  lambda[1] <- 0
  ridge <- .5 * CVXR::sum_squares(sqrt((1 - l1l2) * lambda) * beta)
  lasso <- CVXR::cvxr_norm((l1l2 * lambda) * beta, 1)
  lasso + ridge
}

