## "both" log-likelihood function
bothll <- function (beta, missing, x, weights, estimand, alpha, scratio) {
  score <- navigatedll(beta = beta, x = x, missing = missing, 
                       estimand = estimand, weights = weights, alpha = alpha)
  ps <- c(logistic(x %*% beta))
  missing <- c(missing)
  if (estimand %in% c("ATT", "MO")) {
    cb <- sum(((1 - missing) / (ps - 1) - 
               missing * log(1 - ps) + missing * log(ps)) * weights)
  } else { # ATEcombined
    cb <- sum(((1 - missing) / (ps - 1) - missing / ps + 
              (1 - 2 * missing) * (log(1 - ps) - log(ps))) * weights)
  }
  score * scratio + cb * (1 - scratio)
}

## First derivative of log-likelihood function for the navigated weighting ("both")
bothllgradient <- function (beta, x, missing, estimand, weights, alpha, scratio) {
  ps <- c(logistic(x %*% beta))
  missing <- c(missing)
  if (estimand == "ATEcombined") {
    fdescore <- c(missing - ps) * (ps^alpha + (1 - ps)^alpha) * x * weights
    fdecb <- c(missing - ps) / (ps * (1 - ps)) * x * weights
  } else { # estimand == "ATT" or "MO"
    fdescore <- c(missing - ps) * ps^alpha * x * weights
    fdecb <- c(missing - ps) / (1 - ps) * x * weights
  }
  fde <- scratio * fdescore + (1 - scratio) * fdecb
  apply(fde, 2, sum)
}

bothV <- function (missing, ps, x, outcome, est, N, N1, weights, 
                   estimand, est.weights, alpha, scratio) {
  missing <- c(missing)
  if (estimand == "MO") {
    sumw0 <- sum((1 - missing) / (1 - ps) * weights)
    G0sc <- 
      t(x) %*% 
        (x * 
          (c(1 - missing) * ps / (1 - ps) * (1 - scratio) - 
            c(alpha * missing - (alpha + 1) * ps) * ps^alpha * (1 - ps) * scratio) * 
         weights) / N
    g0sc <- (c(((1 - missing) / (1 - ps)) - 1) * (1 - scratio) - 
              c(missing - ps) * ps^alpha * scratio) * x * weights
    G0mu <- c(apply(ps / (1 - ps) * c(1 - missing) * c(outcome - est) * x * 
                      weights * N / sumw0, 2, sum) / N, 
              -1)
    g0mu <- (outcome - est) * est.weights * N
  } else if (estimand == "ATT") {
    G0sc <- 
      t(x) %*% 
        (-x * 
          (c(1 - missing) * ps / (1 - ps) * (1 - scratio) + 
            c((alpha + 1) * ps - alpha * missing) * ps^alpha * (1 - ps) * scratio) * 
         weights) / N
    g0sc <- -(c(((1 - missing) / (1 - ps)) - 1) * (1 - scratio) + 
                c((1 - missing) - (1 - ps)) * ps^alpha * scratio) * x * weights
    est1 <- sum(missing * outcome * est.weights)
    est0 <- sum((1 - missing) * outcome * est.weights)
    sumw0 <- sum((1 - missing) * ps / (1 - ps) * weights)
    G0mub <- -apply(ps / (1 - ps) * c(1 - missing) * (outcome - est0) * x * 
                      weights * N / sumw0, 2, sum) / N
    G0mu <- c(G0mub, -1)
    g0mu <- 
      (missing * (outcome - est1) * N / N1 - 
      (1 - missing) * (outcome - est0) / (1 - ps) * ps * N / sumw0) * weights
  }  else if (estimand == "ATE") {
    scratio1 <- scratio[1]
    scratio2 <- scratio[2]
    ps1 <- ps[1, ]
    ps2 <- ps[2, ]
    est0 <- sum((1 - missing) * outcome * est.weights)
    est1 <- sum(missing * outcome * est.weights)
    sumw1 <- sum((1 - missing) / (1 - ps1) * weights)
    sumw2 <- sum(missing / ps2 * weights)
    G0sc1 <- 
      t(x) %*% 
       (x * 
        (c(1 - missing) * ps1 / (1 - ps1) * (1 - scratio1) -
         c(alpha * missing - (alpha + 1) * ps1) * ps1^alpha * (1 - ps1) * scratio1) * 
        weights) / N
    g0sc1 <- (c(((1 - missing) / (1 - ps1)) - 1) * (1 - scratio1) -
               c(missing - ps1) * ps1^alpha * scratio1) * x * weights
    G0mu1 <- c(apply(ps1 / (1 - ps1) * c(1 - missing) * c(outcome - est0) * x * 
                      weights * N / sumw1, 2, sum) / N, 
               -1)
    g0mu1 <- (1 - missing) * (outcome - est0) * est.weights * N
    G0sc2 <- 
      t(x) %*% 
       (x * 
         (missing * (1 - ps2) / ps2 * (1 - scratio2) -
         c(alpha * (1 - missing) - 
            (alpha + 1) * (1 - ps2)) * (1 - ps2)^alpha * ps2 * scratio2) * 
         weights) / N
    g0sc2 <- 
      -(c((missing / ps2) - 1) * (1 - scratio2) - 
        c((1 - missing) - (1 - ps2)) * (1 - ps2)^alpha * scratio2) * x * weights
    G0mu2 <- c(-apply((1 - ps2) / ps2 * missing * c(outcome - est1) * x * 
                      weights * N / sumw2, 2, sum) / N, 
               -1)
    g0mu2 <- missing * (outcome - est1) * est.weights * N
  }  else { # ATEcombined
    G0sc <- 
      t(x) %*% 
       (-x * 
         (c(missing - ps)^2 / (ps * (1 - ps)) * (1 - scratio) -
          ((alpha * (missing - ps) - ps) * ps^alpha * (1 - ps) +
            (alpha * (ps - missing) - (1 - ps)) * ps * (1 - ps)^alpha) * scratio) * 
         weights) / N
    g0sc <- 
      (c(missing / ps - (1 - missing) / (1 - ps)) * (1 - scratio) +
       (missing - ps) * (ps^alpha + (1 - ps)^alpha) * scratio) * x * weights
    est1 <- sum(missing * outcome * est.weights)
    est0 <- sum((1 - missing) * outcome * est.weights)
    sumw1 <- sum(missing / ps * weights)
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
    G0 <- as.matrix(rbind(cbind(G0sc, 0), G0mu))
    g0 <- as.matrix(cbind(g0sc, g0mu))
    Sigma <- t(g0) %*% g0 / N
    MASS::ginv(G0) %*% Sigma %*% t(MASS::ginv(G0)) / N
  } else { # ATE
    zero <- matrix(0, nrow = ncol(x) + 1, ncol = ncol(x) + 1)
    G01 <- as.matrix(rbind(cbind(G0sc1, 0), G0mu1))
    G02 <- as.matrix(rbind(cbind(G0sc2, 0), G0mu2))
    G03 <- as.matrix(rbind(cbind(G01, zero), cbind(zero, G02)))
    G0 <- as.matrix(rbind(cbind(G03, 0), 
                          c(rep(0, ncol(x)), -1, rep(0, ncol(x)), 1, -1)))
    g0 <- as.matrix(cbind(g0sc1, g0mu1, g0sc2, g0mu2, 0))
    Sigma <- t(g0) %*% g0 / N
    (MASS::ginv(G0) %*% Sigma %*% t(MASS::ginv(G0)) / N)[
                                          -c(ncol(x) + 1, ncol(x) * 2 + 2), 
                                          -c(ncol(x) + 1, ncol(x) * 2 + 2)]
  }
}
