## cb log-likelihood function
cbll <- function (beta, missing, x, weights, estimand) {
  ps <- c(logistic(x %*% beta))
  missing <- c(missing)
  if (estimand %in% c("ATT", "MO")) {
    sum(((1 - missing) / (ps - 1) - missing * log(1 - ps) + 
         missing * log(ps)) * weights)
  } else { # ATEcombined
    sum(((1 - missing) / (ps - 1) - missing / ps + 
          (1 - 2 * missing) * (log(1 - ps) - log(ps))) * weights)
  }
}

## First derivative of log-likelihood function for the navigated weighting (cb)
cbllgradient <- function (beta, x, missing, estimand, weights) {
  ps <- c(logistic(x %*% beta))
  missing <- c(missing)
  if (estimand == "ATEcombined") {
    fde <- c(missing - ps) / (ps * (1 - ps)) * x * weights
  } else { # estimand == "ATT" or "MO"
    fde <- c(missing - ps) / (1 - ps) * x * weights
  }
  apply(fde, 2, sum)
}

## cb variance estimation
cbV <- function (missing, ps, x, outcome, est, N, N1, weights, 
                 estimand, est.weights) {
  missing <- c(missing)
  if (estimand == "MO") {
    sumw0 <- sum((1 - missing) / (1 - ps) * weights)
    G0cb <- t(x) %*% (-x * c(1 - missing) * ps / (1 - ps) * weights) / N
    g0cb <- c(missing - ps) / (1 - ps) * x * weights
    G0mu <- c(apply(ps / (1 - ps) * c(1 - missing) * c(outcome - est) * x * 
                      weights * N / sumw0, 2, sum) / N, 
              -1)
    g0mu <- (outcome - est) * est.weights * N
  } else if (estimand == "ATT") {
    G0cb <- t(x) %*% (-x * c(1 - missing) * ps / (1 - ps) * weights) / N
    g0cb <- -c(((1 - missing) / (1 - ps)) - 1) * x * weights
    est1 <- sum(missing * outcome * est.weights)
    est0 <- sum((1 - missing) * outcome * est.weights)
    sumw0 <- sum((1 - missing) * ps / (1 - ps) * weights)
    G0mub <- -apply(ps / (1 - ps) * c(1 - missing) * (outcome - est0) * x * 
                      weights * N / sumw0, 2, sum) / N
    G0mu <- c(G0mub, -1)
    g0mu <- (missing * (outcome - est1) * N / N1 - 
             (1 - missing) * (outcome - est0) / (1 - ps) * ps * N / sumw0) * 
            weights
  }  else if (estimand == "ATE") {
    ps1 <- ps[1, ]
    ps2 <- ps[2, ]
    est0 <- sum((1 - missing) * outcome * est.weights)
    est1 <- sum(missing * outcome * est.weights)
    sumw1 <- sum((1 - missing) / (1 - ps1) * weights)
    sumw2 <- sum(missing / ps2 * weights)
    G0cb1 <- t(x) %*% (-x * c(1 - missing) * ps1 / (1 - ps1) * weights) / N
    g0cb1 <- c(missing - ps1) / (1 - ps1) * x * weights
    G0mu1 <- c(apply(ps1 / (1 - ps1) * c(1 - missing) * c(outcome - est0) * x * 
                      weights * N / sumw1, 2, sum) / N, 
               -1)
    g0mu1 <- (1 - missing) * (outcome - est0) * est.weights * N
    G0cb2 <- t(x) %*% (-x * missing * (1 - ps2) / ps2 * weights) / N
    g0cb2 <- c(missing - ps2) / ps2 * x * weights
    G0mu2 <- c(-apply((1 - ps2) / ps2 * missing * c(outcome - est1) * x * 
                      weights * N / sumw2, 2, sum) / N, 
               -1)
    g0mu2 <- missing * (outcome - est1) * est.weights * N
  }  else { # ATEcombined
    G0cb <- t(x) %*% (-x * c(missing - ps)^2 / (ps * (1 - ps)) * weights) / N
    g0cb <- c(missing / ps - (1 - missing) / (1 - ps)) * x * weights
    est1 <- sum(missing * outcome * est.weights)
    est0 <- sum((1 - missing) * outcome * est.weights)
    sumw1 <- sum(missing / ps  * weights)
    sumw0 <- sum((1 - missing) / (1 - ps) * weights)
    G0mub <- -apply((1 - ps) / ps * missing * (outcome - est1) * x * 
                      weights * N / sumw1 + 
                    ps / (1 - ps) * (1 - missing) * (outcome - est0) * x * 
                      weights * N / sumw0, 2, sum) / N
    G0mu <- c(G0mub, -1)
    g0mu <- (missing * (outcome - est1) / ps * N / sumw1 - 
             (1 - missing) * (outcome - est0) / (1 - ps) * N / sumw0) * weights
  }
  if (estimand != "ATE") {
    G0 <- as.matrix(rbind(cbind(G0cb, 0), G0mu))
    g0 <- as.matrix(cbind(g0cb, g0mu))
    Sigma <- t(g0) %*% g0 / N
    MASS::ginv(G0) %*% Sigma %*% t(MASS::ginv(G0)) / N
  } else { # ATE
    zero <- matrix(0, nrow = ncol(x) + 1, ncol = ncol(x) + 1)
    G01 <- as.matrix(rbind(cbind(G0cb1, 0), G0mu1))
    G02 <- as.matrix(rbind(cbind(G0cb2, 0), G0mu2))
    G03 <- as.matrix(rbind(cbind(G01, zero), cbind(zero, G02)))
    G0 <- as.matrix(rbind(cbind(G03, 0), 
                          c(rep(0, ncol(x)), -1, rep(0, ncol(x)), 1, -1)))
    g0 <- as.matrix(cbind(g0cb1, g0mu1, g0cb2, g0mu2, 0))
    Sigma <- t(g0) %*% g0 / N
    (MASS::ginv(G0) %*% Sigma %*% t(MASS::ginv(G0)) / N)[
                                          -c(ncol(x) + 1, ncol(x) * 2 + 2), 
                                          -c(ncol(x) + 1, ncol(x) * 2 + 2)]
  }
}
