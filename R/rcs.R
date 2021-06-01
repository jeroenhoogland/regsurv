#' Generate restricted cubic spline (rcs)  basis matrix
#'
#' Takes a numeric vector and a sequence of knots to generate restricted cubic splines basis
#' columns based on a truncated power series (as described in section 2.4.5 of the second edition of
#' Regression Modeling Strategies by Frank E. Harrrell, Jr).
#'
#' @param x numeric vector
#' @param knots numeric vector including at least 3 unique knot locations of increasing order,
#' with the outer ones indicating the boundary knots.
#'
#' @return matrix of basis columns
#' @export
#'
#' @examples
#' rcs(1:25)
rcs <- function(x, knots = stats::quantile(x, probs=c(0,0.1,0.5,0.9,1))){

  if(!is.numeric(x)){
    stop("x should be a numeric vector")
  }

  if(length(x) < 2){
    stop("x should have length >= 2")
  }

  if(!is.numeric(knots)){
    stop("k should be a numeric vector")
  }

  if(length(unique(knots)) < 3){
    stop("provide a minium of 3 unique knots")
  }

  if(!isTRUE(all.equal(order(knots), 1:length(knots)))){
    stop("provide knots in increasing order")
  }

  k <- length(knots)

  if(length(unique(knots)) != k){
    warning("Duplicate knots; reduced to unique knots")
    knots <- unique(knots)
    k <- length(knots)
  }

  X <- matrix(NA, nrow = length(x), ncol = k-1)
  X[ , 1] <- x

  for(j in 1:(k-2)){
    scale1 <- (knots[k] - knots[j]) / (knots[k] - knots[k-1])
    scale2 <- (knots[k-1] - knots[j]) / (knots[k] - knots[k-1])
    X[ ,j+1] <- pmax(0, x - knots[j])^3 -
      pmax(0, x - knots[k-1])^3 * scale1 +
      pmax(0, x - knots[k])^3 * scale2
  }
  tau <- (knots[k] - knots[1])^2
  X[ ,2:(k-1)] <- X[ ,2:(k-1)] / tau # scale such that all cols are on scale of X

  return(X)
}


# derivative of rcs()
drcs <- function(x, knots=stats::quantile(x, probs=c(0,0.1,0.5,0.9,1))){
  k <- length(knots)
  X <- matrix(NA, nrow = length(x), ncol = k-1)
  X[ , 1] <- 1
  for(j in 1:(k-2)){
    scale1 <- (knots[k] - knots[j]) / (knots[k] - knots[k-1])
    scale2 <- (knots[k-1] - knots[j]) / (knots[k] - knots[k-1])
    X[ ,j+1] <- 3*pmax(0, x - knots[j])^2 -
      3*pmax(0, x - knots[k-1])^2 * scale1 +
      3*pmax(0, x - knots[k])^2 * scale2
  }
  tau <- (knots[k] - knots[1])^2
  X[ ,2:(k-1)] <- X[ ,2:(k-1)] / tau
  return(X)
}
