## gmm function
gmm <- function (beta, missing, x, N, N1, weights, estimand, alpha) {
  missing <- c(missing)
  ps <- c(logistic(x %*% beta))
  if (estimand %in% c("ATT", "MO")) {
    g0score <- c(missing - ps) * ps^alpha * x * (N / N1) * weights
    g0cb <- c(((1 - missing) / (1 - ps)) - 1) * x * (N / N1) * weights
    g <- as.matrix(cbind(g0score, g0cb))
  } else { # ATEcombined
    g0score <- (missing - ps) * (ps^alpha + (1 - ps)^alpha) * x * weights
    g0cb <- c(missing / ps - (1 - missing) / (1 - ps)) * x * weights
    g <- as.matrix(cbind(g0score, g0cb))
  }
  gmm.bread <- apply(g, 2, mean)
  invW <- gmmW(missing = missing, ps = ps, x = x, N = N, N1 = N1, 
               weights = weights, estimand = estimand, alpha = alpha)
  t(gmm.bread) %*% invW %*% gmm.bread * sqrt(N)
}

## gmm optimal covariance estimation
gmmW <- function (missing, ps, x, N, N1, weights, estimand, alpha) {
  if (estimand %in% c("ATT", "MO")) {
    Wul <- t(ps^(2 * alpha + 1) * (1 - ps) * (N / N1)^2 * weights * x) %*% x
    Wur <- t(ps^(alpha + 1) * (N / N1)^2 * weights * x) %*% x
    Wll <- Wur
    Wlr <- t(ps / (1 - ps) * (N / N1)^2 * weights * x) %*% x
    W <- as.matrix(rbind(cbind(Wul, Wur), cbind(Wll, Wlr))) / N
  } else { # ATEcombined
    Wul <- t(ps * (1 - ps) * (ps^alpha + (1 - ps)^alpha)^2 * weights * x) %*% x
    Wur <- t((ps^alpha + (1 - ps)^alpha) * weights * x) %*% x
    Wll <- Wur
    Wlr <- t(ps^{-1} / (1 - ps) * weights * x) %*% x
    W <- as.matrix(rbind(cbind(Wul, Wur), cbind(Wll, Wlr))) / N
  }
  MASS::ginv(W)
}

## gmm variance estimation
gmmV <- function (missing, ps, x, outcome, est, N, N1, weights, 
                  estimand, est.weights, alpha) {
  missing <- c(missing)
  if (estimand == "MO") {
    sumw0 <- sum((1 - missing) / (1 - ps) * weights)
    G0score <- 
      t(x) %*% 
        (x * c(alpha * missing - (alpha + 1) * ps) * ps^alpha * (1 - ps) * 
          weights) / N
    G0cb <- t(x) %*% (x * c(1 - missing) * ps / (1 - ps) * weights) / N
    g0score <- c(missing - ps) * ps^alpha * x * weights
    g0cb <- c(((1 - missing) / (1 - ps)) - 1) * x * weights    
    G0mub <- apply(ps / (1 - ps) * c(1 - missing) * c(outcome - est) * x * 
                   weights * N / sumw0, 2, sum) / N
    G0mu <- c(G0mub, -1)
    g0mu <- (outcome - est) * est.weights * N
  } else if (estimand == "ATT") {
    G0score <- 
      t(x) %*% 
        (-x * c(alpha * missing - (alpha + 1) * ps) * ps^alpha * (1 - ps) * 
          weights) / N
    G0cb <- t(x) %*% (-x * c(1 - missing) * ps / (1 - ps) * weights) / N
    g0score <- -c((1 - ps) - (1 - missing)) * ps^alpha * x * weights
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
  } else if (estimand == "ATE") {
    ps1 <- ps[1, ]
    ps2 <- ps[2, ]
    est0 <- sum((1 - missing) * outcome * est.weights)
    est1 <- sum(missing * outcome * est.weights)
    sumw1 <- sum((1 - missing) / (1 - ps1) * weights)
    sumw2 <- sum(missing / ps2) * weights
    G0score1 <- 
      t(x) %*% 
        (x * c(alpha * missing - (alpha + 1) * ps1) * ps1^alpha * (1 - ps1) * 
          weights) / N
    G0cb1 <- t(x) %*% (x * c(1 - missing) * ps1 / (1 - ps1) * weights) / N
    g0score1 <- c(missing - ps1) * ps1^alpha * x * weights
    g0cb1 <- c(((1 - missing) / (1 - ps1)) - 1) * x * weights    
    G0mu1 <- c(apply(ps1 / (1 - ps1) * c(1 - missing) * c(outcome - est0) * x * 
                weights * N / sumw1, 2, sum) / N, 
               -1)
    g0mu1 <- (1 - missing) * (outcome - est0) * est.weights * N
    G0score2 <- 
      t(x) %*% 
        (x * c(alpha * (1 - missing) - 
                (alpha + 1) * (1 - ps2)) * 
              (1 - ps2)^alpha * ps2 * weights) / N
    g0score2 <- c((1 - missing) - (1 - ps2)) * (1 - ps2)^alpha * x * weights
    G0cb2 <- t(x) %*% (x * missing * (1 - ps2) / ps2 * weights) / N
    g0cb2 <- c((missing / ps2) - 1) * x * weights
    G0mu2 <- c(apply((1 - ps2) / ps2 * missing * c(outcome - est1) * x * 
                      weights * N / sumw2, 2, sum) / N, 
               -1)
    g0mu2 <- missing * (outcome - est1) * est.weights * N
  } else { # ATEcombined
    G0score <- 
      t(x) %*% 
        (x * ((alpha * (missing - ps) - ps) * ps^alpha * (1 - ps) +
                (alpha * (ps - missing) - (1 - ps)) * ps * (1 - ps)^alpha) * 
          weights) / N
    G0cb <- t(x) %*% (-x * c(missing - ps)^2 / (ps * (1 - ps)) * weights) / N
    g0score <- (missing - ps) * (ps^alpha + (1 - ps)^alpha) * x * weights
    g0cb <- c(missing - ps) / (ps * (1 - ps)) * x * weights
    est1 <- sum(missing * outcome * est.weights)
    est0 <- sum((1 - missing) * outcome * est.weights)
    sumw1 <- sum(missing / ps * weights)
    sumw0 <- sum((1 - missing) / (1 - ps) * weights)
    G0mub <- 
      -apply((1 - ps) / ps * missing * (outcome - est1) * x * 
              weights * N / sumw1 + 
             ps / (1 - ps) * (1 - missing) * (outcome - est0) * x * 
               weights * N / sumw0, 2, sum) / N
    G0mu <- c(G0mub, -1)
    g0mu <- (missing * (outcome - est1) / ps * N / sumw1 - 
             (1 - missing) * (outcome - est0) / (1 - ps) * N / sumw0) * weights
  }
  if (estimand != "ATE") {
    G0 <- as.matrix(rbind(G0score, G0cb))
    G <- as.matrix(rbind(cbind(G0, 0), G0mu))
    g <- as.matrix(cbind(g0score, g0cb, g0mu))
    Sigma <- t(g) %*% g / N
    invW0 <- gmmW(missing = missing, ps = ps, x = x, N = N, N1 = N1, 
                  weights = weights, estimand = estimand, alpha = alpha)
    invW <- as.matrix(rbind(cbind(invW0, 0), 0))
    invW[nrow(invW), ncol(invW)] <- 1
    bread <- MASS::ginv(t(G) %*% invW %*% G)
    bread %*% t(G) %*% invW %*% Sigma %*% invW %*% G %*% bread / N
  } else { # ATE
    zero <- matrix(0, nrow = ncol(x) * 2 + 1, ncol = ncol(x) + 1)
    G01 <- as.matrix(rbind(cbind(rbind(G0score1, G0cb1), 0), G0mu1))
    G02 <- as.matrix(rbind(cbind(rbind(G0score2, G0cb2), 0), G0mu2))
    G03 <- as.matrix(rbind(cbind(G01, zero), cbind(zero, G02)))
    G <- as.matrix(
          rbind(cbind(G03, 0), c(rep(0, ncol(x)), -1, rep(0, ncol(x)), 1, -1)))
    g0 <- as.matrix(cbind(g0score1, g0cb1, g0mu1, g0score2, g0cb2, g0mu2, 0))
    Sigma <- t(g0) %*% g0 / N
    invW01 <- gmmW(missing = missing, ps = ps1, x = x, N = N, N1 = N1, 
                   weights = weights, estimand = "MO", alpha = alpha)
    invW1 <- as.matrix(rbind(cbind(invW01, 0), 0))
    invW1[nrow(invW1), ncol(invW1)] <- 1
    invW02 <- gmmW(missing = 1 - missing, ps = ps2, x = x, N = N, N1 = N - N1, 
                   weights = weights, estimand = "MO", alpha = alpha)
    invW2 <- as.matrix(rbind(cbind(invW02, 0), 0))
    invW2[nrow(invW2), ncol(invW2)] <- 1
    invWzero <- matrix(0, nrow = ncol(x) * 2 + 1, ncol = ncol(x) * 2 + 1)
    invW <- rbind(cbind(invW1, invWzero, 0), cbind(invWzero, invW2, 0), 0)
    invW[nrow(invW), ncol(invW)] <- 1
    bread <- MASS::ginv(t(G) %*% invW %*% G)
    (bread %*% t(G) %*% invW %*% Sigma %*% invW %*% G %*% bread / N)[
      -c(ncol(x) + 1, ncol(x) * 2 + 2), 
      -c(ncol(x) + 1, ncol(x) * 2 + 2)]
  }
}
