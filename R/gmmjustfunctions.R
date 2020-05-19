## just-identified gmm function
gmmjust <- function (beta, missing, x, N, N1, weights, estimand) {
  ps <- c(logistic(x %*% beta))
  if (estimand %in% c("ATT", "MO")) {
    g0cb <- c(((1 - missing) / (1 - ps)) - 1) * x * (N / N1) * weights
  } else { # ATE
    g0cb <- c(missing / ps - (1 - missing) / (1 - ps)) * x * weights
  }
  gmm.bread <- apply(g0cb, 2, mean)
  t(gmm.bread) %*% gmm.bread
}

## just-identified gmm gradient
gmmgradientjust <- function (beta, missing, x, N, N1, weights, estimand) {
  ps <- c(logistic(x %*% beta))
  if (estimand %in% c("ATT", "MO")) {
    g0cb <- c(((1 - missing) / (1 - ps)) - 1) * x * (N / N1) * weights    
    G0cb <- t(x) %*% (x * c(1 - missing) * ps / (1 - ps) * (N / N1) * weights)
  } else { # ATE
    g0cb <- c(missing / ps - (1 - missing) / (1 - ps)) * x * weights
    G0cb <- t(x) %*% (-x * c(missing - ps)^2 / (ps * (1 - ps)) * weights)
  }
  2 * apply(g0cb %*% G0cb, 2, sum) / N^2
}
