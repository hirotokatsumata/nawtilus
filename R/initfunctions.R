## initial values for estimation
init <- function (result, method, estimand, missing, x, alpha, N) {
  if (method %in% c("score", "both")) {
    if (alpha > 1) {
      alpha <- sqrt(alpha)
    }
  } else if (method == "cb") {
    alpha <- -1 / 2
  }  
  if (estimand == "ATEcombined") {
    if (method != "both") {
      weights <- init.omega(lp = result$linear.predictors, estimand = "ATEcombined", 
                            method = method, alpha = alpha, N = N)
    } else { # both
      weights1 <- init.omega(lp = result$linear.predictors, estimand = "ATEcombined", 
                             method = method, alpha = alpha, N = N)
      weights2 <- init.omega(lp = result$linear.predictors, estimand = "ATEcombined", 
                             method = method, alpha = -1 / 2, N = N)
      weights <- (weights1 + weights2) / 2
    }
    coef <- c(stats::glm(formula = missing ~ -1 + x, family = stats::quasibinomial, weights = weights)$coef)
  } else { # ATT and MO
    weights <- init.omega(lp = result$linear.predictors / 3, estimand = estimand, 
                          method = method, alpha = alpha, N = N)
    result <- stats::glm(formula = missing ~ -1 + x, family = stats::quasibinomial, weights = weights)
    weights <- init.omega(lp = result$linear.predictors * 2 / 3, estimand = estimand, 
                          method = method, alpha = alpha, N = N)
    coef <- c(stats::glm(formula = missing ~ -1 + x, family = stats::quasibinomial, weights = weights)$coef)
  }
  coef
}

## init internal
init.omega <- function(lp, estimand, method, alpha, N) {
  mintol <- 0.01
  maxtol <- N / 50
  if (estimand == "ATEcombined") {
    init <- (logistic(lp)^alpha + (1 - logistic(lp))^alpha) / 2
    init <- init / mean(init)
    init[init < mintol] <- mintol
    init[init > maxtol] <- maxtol
  } else { # ATT and MO
    if (method == "score") {
      init <- logistic(lp)^alpha
      init <- init / mean(init)
      init[init < mintol] <- mintol
      init[init > maxtol] <- maxtol
    } else if (method == "cb") {
      init <- (1 - logistic(lp))^(-1 / 2)
      init <- init / mean(init)
      init[init < mintol] <- mintol
      init[init > maxtol] <- maxtol
    } else { # both
      init1 <- logistic(lp)^alpha
      init1 <- init1 / mean(init1)
      init1[init1 < mintol] <- mintol
      init1[init1 > maxtol] <- maxtol
      init2 <- (1 - logistic(lp))^(-1 / 2)
      init2 <- init2 / mean(init2)
      init2[init2 < mintol] <- mintol
      init2[init2 > maxtol] <- maxtol
      init <- (init1 + init2) / 2
    }
  }
  init
}
