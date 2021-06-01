# gaussquad function do not explicitly state dependencies, which
# causes trouble during development with load_all()
# this provides a temporary workaround
legendre.quadrature.rules <- function (n, normalized = FALSE)
{
  recurrences <- orthopolynom::legendre.recurrences(n, normalized)
  inner.products <- orthopolynom::legendre.inner.products(n)
  return(quadrature.rules(recurrences, inner.products))
}

quadrature.rules <- function (recurrences, inner.products)
{
  np1 <- nrow(recurrences)
  n <- np1 - 1
  rules <- as.list(rep(NULL, n))
  monic.recurrences <- orthopolynom::monic.polynomial.recurrences(recurrences)
  matrices <- orthopolynom::jacobi.matrices(monic.recurrences)
  matrix.eigens <- lapply(matrices, eigen)
  roots <- orthopolynom::polynomial.roots(monic.recurrences)
  h.0 <- inner.products[1]
  for (k in 1:n) {
    values <- matrix.eigens[[k]]$values
    vectors <- matrix.eigens[[k]]$vectors
    x <- values
    w <- rep(0, k)
    for (j in 1:k) {
      v.j <- vectors[1, j]
      w[j] <- h.0 * v.j * v.j
    }
    rule <- data.frame(cbind(x, w))
    names(rule) <- c("x", "w")
    rules[[k]] <- rule
  }
  return(rules)
}
