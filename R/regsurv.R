regsurv <- function(prep, penpars, l1l2, lambda.grid=NULL, print=FALSE, ...){

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
  betahat <- sapply(sol, function(x) x$getValue(beta))

  return(
    structure(
      list(optimal=optimal,
           lambda.grid=lambda.grid,
           obj.value=obj.value,
           betahat=betahat,
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

