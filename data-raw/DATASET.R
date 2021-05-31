## code to prepare `DATASET` dataset goes here

library(simsurv)

set.seed(1)

# Weibul DGM with 4 main covariate effects, 2 of which also interact with time (i.e. are
# time-varying)
scalepar <- c(0.05)
shapepar <- c(1.8)
betas = c("x1" = 1, "x2"=.5, "x3"=.5, "x4"=-.5)
tde <- c("x1" = -0.1, "x2"=0.1)
tdefunction <- function(t) t

n <- 250   # sample size
q <- 5     # number of noise variables
maxt <- 10 # censoring time
ncols <- length(betas) + q

X <- matrix(stats::rnorm(n*ncols, 0, 1),
            nrow=n,
            ncol=ncols,
            dimnames=list(NULL, paste0("x", 1:ncols)))

simdata <- simsurv(dist = "weibull", lambdas = scalepar, gammas = shapepar,
                   mixture = FALSE, tde=tde, tdefunction = tdefunction,
                   x = data.frame(X),
                   betas=betas, maxt = maxt)

simdata <- cbind(simdata, X)

usethis::use_data(simdata, overwrite = TRUE)
