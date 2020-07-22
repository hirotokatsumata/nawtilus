## Log-likelihood function for the navigated weighting
navigatedll <- function (beta, x, missing, estimand, weights, alpha) {
  ps <- c(logistic(x %*% beta))
  missing <- c(missing)
  if (estimand == "ATEcombined") {
    if (alpha == 0) {
      sum((missing * log(ps) + (1 - missing) * log(1 - ps)) * 2 * weights)
    } else if (alpha == 1) {
      sum((missing * log(ps) + (1 - missing) * log(1 - ps)) * weights)
    } else if (alpha == 2) {
      sum((ps^2 - 2 * missing * ps + missing * log(ps) + (1 - missing) * log(1 - ps)) * weights)
    } else if (alpha == 3) {
      sum((ps^2 * 3 / 2 - 3 * missing * ps + missing * log(ps) + (1 - missing) * log(1 - ps)) * weights)
    } else if (alpha == 4) {
      sum((ps^4 / 2 - (missing + 1) * ps^3 * 2 / 3 + (missing + 2) * ps^2 - 4 * missing * ps + 
            missing * log(ps) + (1 - missing) * log(1 - ps)) * weights)
    } else if (alpha == 5) {
      sum((ps^4 * 5 / 4 - (missing + 1) * ps^3 * 5 / 3 + (missing + 1) * ps^2 * 5 / 2 - 5 * missing * ps + 
            missing * log(ps) + (1 - missing) * log(1 - ps)) * weights)
    } else if (alpha %in% seq(0.5, 4.5, by = 0.5)) {
      sum(((alpha >= 4.5) * ((ps^4.5 / 4.5 - missing * ps^3.5 / 3.5) + 
                             ((1 - ps)^4.5 / 4.5 - (1 - missing) * (1 - ps)^3.5 / 3.5)) +
           (alpha >= 3.5) * ((ps^3.5 / 3.5 - missing * ps^2.5 / 2.5) +
                              ((1 - ps)^3.5 / 3.5 - (1 - missing) * (1 - ps)^2.5 / 2.5)) +
           (alpha >= 2.5) * ((ps^2.5 / 2.5 - missing * ps^1.5 / 1.5) +
                              ((1 - ps)^2.5 / 2.5 - (1 - missing) * (1 - ps)^1.5 / 1.5)) +
           (alpha >= 1.5) * ((ps^1.5 / 1.5 - missing * ps^0.5 / 0.5) +
                              ((1 - ps)^1.5 / 1.5 - (1 - missing) * (1 - ps)^0.5 / 0.5)) +
           ps^0.5 / 0.5 + (1 - ps)^0.5 / 0.5 + 
            (1 - missing) * log(1 - ps^0.5) + missing * log(ps) -
            (1 - missing) * log(1 + ps^0.5) - 2 * missing * log(1 + (1 - ps)^0.5)) * weights)
    } else {
    sum((ps^alpha * (missing / alpha - 
                      (1 - missing) * ps * hypergeo::hypergeo(1, 1 + alpha, 2 + alpha, ps) / (alpha + 1)) +
         (1 - ps)^alpha * (-missing * (1 - ps) / ps * hypergeo::hypergeo(1, 1, 1 - alpha, 1 / ps) + 
             (1 - missing)) / alpha) * weights)
    }
  } else { # estimand == "ATT" or "MO"
    if (alpha == 0) {
      sum((missing * log(ps) + (1 - missing) * log(1 - ps)) * weights)
    } else if (alpha %in% 1:5) {
      sum(((alpha >= 5) * (ps^5 / 5 - missing * ps^4 / 4) +
           (alpha >= 4) * (ps^4 / 4 - missing * ps^3 / 3) +
           (alpha >= 3) * (ps^3 / 3 - missing * ps^2 / 2) +
           (alpha >= 2) * (ps^2 / 2 - missing * ps) +
           ps + (1 - missing) * log(1 - ps)) * weights)
    } else if (alpha %in% seq(0.5, 4.5, by = 0.5)) {
      sum(((alpha >= 4.5) * (ps^4.5 / 4.5 - missing * ps^3.5 / 3.5) +
           (alpha >= 3.5) * (ps^3.5 / 3.5 - missing * ps^2.5 / 2.5) +
           (alpha >= 2.5) * (ps^2.5 / 2.5 - missing * ps^1.5 / 1.5) +
           (alpha >= 1.5) * (ps^1.5 / 1.5 - missing * ps^0.5 / 0.5) +
           ps^0.5 / 0.5 + (1 - missing) * log(1 - ps^0.5) - (1 - missing) * log(1 + ps^0.5)) * weights)
    } else {
    sum((ps^alpha * (missing / alpha - 
                      (1 - missing) * ps * 
                        hypergeo::hypergeo(1, 1 + alpha, 2 + alpha, ps) / (alpha + 1))) * 
          weights)
    }
  }
}

## First derivative of log-likelihood function for the navigated weighting
navigatedgradient <- function (beta, x, missing, estimand, weights, alpha) {
  ps <- c(logistic(x %*% beta))
  if (estimand == "ATEcombined") {
    fde <- c(missing - ps) * (ps^alpha + (1 - ps)^alpha) * x * weights
  } else { # estimand == "ATT" or "MO"
    fde <- c(missing - ps) * ps^alpha * x * weights
  }
  apply(fde, 2, sum)
}

## navigated weighting (score): variance estimation
navigatedV <- function (missing, ps, x, outcome, est, N, N1, weights, estimand, 
                        est.weights, alpha) {
  missing <- c(missing)
  if (estimand == "MO") {
    sumw0 <- sum((1 - missing) / (1 - ps) * weights)
    G0score <- t(x) %*% (x * c(alpha * missing - (alpha + 1) * ps) * ps^alpha * (1 - ps) * weights) / N
    g0score <- c(missing - ps) * ps^alpha * x * weights
    G0mu <- c(apply(ps / (1 - ps) * c(1 - missing) * c(outcome - est) * x * weights * N / sumw0, 2, sum) / N, 
              -1)
    g0mu <- (outcome - est) * est.weights * N
  }  else if (estimand == "ATT") {
    G0score <- t(x) %*% (x * c(alpha * missing - (alpha + 1) * ps) * ps^alpha * (1 - ps) * weights) / N
    g0score <- c((1 - ps) - (1 - missing)) * ps^alpha * x * weights
    est1 <- sum(missing * outcome * est.weights)
    est0 <- sum((1 - missing) * outcome * est.weights)
    sumw0 <- sum((1 - missing) * ps / (1 - ps) * weights)
    G0mub <- -apply(ps / (1 - ps) * c(1 - missing) * (outcome - est0) * x * weights * N / sumw0, 2, sum) / N
    G0mu <- c(G0mub, -1)
    g0mu <- (missing * (outcome - est1) * N / N1 - (1 - missing) * (outcome - est0) / (1 - ps) * ps * N / sumw0) * weights
  }  else if (estimand == "ATE") { # ATE
    ps1 <- ps[1, ]
    ps2 <- ps[2, ]
    est0 <- sum((1 - missing) * outcome * est.weights)
    est1 <- sum(missing * outcome * est.weights)
    sumw1 <- sum((1 - missing) / (1 - ps1) * weights)
    sumw2 <- sum(missing / ps2) * weights
    G0score1 <- t(x) %*% (x * c(alpha * missing - (alpha + 1) * ps1) * ps1^alpha * (1 - ps1) * weights) / N
    g0score1 <- c(missing - ps1) * ps1^alpha * x * weights
    G0mu1 <- c(apply(ps1 / (1 - ps1) * c(1 - missing) * c(outcome - est0) * x * weights * N / sumw1, 2, sum) / N, 
               -1)
    g0mu1 <- (1 - missing) * (outcome - est0) * est.weights * N
    G0score2 <- t(x) %*% (x * c(alpha * (1 - missing) - (alpha + 1) * (1 - ps2)) * (1 - ps2)^alpha * ps2 * weights) / N
    g0score2 <- c(missing - ps2) * (1 - ps2)^alpha * x * weights
    G0mu2 <- c(-apply((1 - ps2) / ps2 * missing * c(outcome - est1) * x * weights * N / sumw2, 2, sum) / N, 
               -1)
    g0mu2 <- missing * (outcome - est1) * est.weights * N
  }  else { # ATEcombined
    G0score <- t(x) %*% (x * ((alpha * (missing - ps) - ps) * ps^alpha * (1 - ps) +
                          (alpha * (ps - missing) - (1 - ps)) * ps * (1 - ps)^alpha) * weights) / N
    g0score <- (missing - ps) * (ps^alpha + (1 - ps)^alpha) * x * weights
    est1 <- sum(missing * outcome * est.weights)
    est0 <- sum((1 - missing) * outcome * est.weights)
    sumw1 <- sum(missing / ps * weights)
    sumw0 <- sum((1 - missing) / (1 - ps) * weights)
    G0mub <- -apply((1 - ps) / ps * missing * (outcome - est1) * x * weights * N / sumw1 + 
                    ps / (1 - ps) * (1 - missing) * (outcome - est0) * x * weights * N / sumw0, 2, sum) / N
    G0mu <- c(G0mub, -1)
    g0mu <- (missing * (outcome - est1) / ps * N / sumw1 - 
             (1 - missing) * (outcome - est0) / (1 - ps) * N / sumw0) * weights
  }
  if (estimand != "ATE") {
    G0 <- as.matrix(rbind(cbind(G0score, 0), G0mu))
    g0 <- as.matrix(cbind(g0score, g0mu))
    Sigma <- t(g0) %*% g0 / N
    MASS::ginv(G0) %*% Sigma %*% t(MASS::ginv(G0)) / N
  } else { # ATE
    zero <- matrix(0, nrow = ncol(x) + 1, ncol = ncol(x) + 1)
    G01 <- as.matrix(rbind(cbind(G0score1, 0), G0mu1))
    G02 <- as.matrix(rbind(cbind(G0score2, 0), G0mu2))
    G03 <- as.matrix(rbind(cbind(G01, zero), cbind(zero, G02)))
    G0 <- as.matrix(rbind(cbind(G03, 0), c(rep(0, ncol(x)), -1, rep(0, ncol(x)), 1, -1)))
    g0 <- as.matrix(cbind(g0score1, g0mu1, g0score2, g0mu2, 0))
    Sigma <- t(g0) %*% g0 / N
    (MASS::ginv(G0) %*% Sigma %*% t(MASS::ginv(G0)) / N)[-c(ncol(x) + 1, ncol(x) * 2 + 2), 
                                                         -c(ncol(x) + 1, ncol(x) * 2 + 2)]
  }
}
