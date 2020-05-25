nawt0 <- function (outcome, estimand = "ATT", method = "score", 
                   missing = missing, x = x, N = N, weights = NULL, 
                   alpha = 2, twostep = TRUE, varcov = TRUE) {
  scratio <- NULL
  weights <- weights / mean(weights)
  N1 <- sum(missing * weights)
  names.x <- colnames(x)
  names.x[apply(x, 2, stats::sd) == 0] <- "(Intercept)"
  x0 <- x
  svdx <- svd(x)
  x <- svdx$u
  estimand0 <- estimand
  if (estimand == "ATC") {
    estimand <- "ATT"
    missing <- 1 - missing
  }
  result.vanilla <- stats::glm(formula = missing ~ -1 + x, 
                               family = stats::quasibinomial, weights = weights)
  est.weights <- numeric(N)
  if (method == "cb") {
    if (estimand == "ATE") {
      initval1 <- init(result = result.vanilla, method = method, estimand = "MO", 
                       missing = missing, x = x, alpha = 0, N = N)
      result.vanilla2 <- 
        stats::glm(formula = I(1 - missing) ~ -1 + x, 
                   family = stats::quasibinomial, weights = weights)
      initval2 <- init(result = result.vanilla2, method = method, estimand = "MO",
                       missing = 1 - missing, x = x, alpha = 0, N = N)
      result1 <- stats::optim(par = initval1, 
                              fn = cbll, 
                              gr = cbllgradient,
                              method = "BFGS", 
                              control = list(fnscale = -1, trace = FALSE), 
                              hessian = FALSE,
                              missing = missing, x = x, weights = weights, 
                              estimand = "MO")
      result2 <- stats::optim(par = initval2, 
                              fn = cbll, 
                              gr = cbllgradient,
                              method = "BFGS", 
                              control = list(fnscale = -1, trace = FALSE), 
                              hessian = FALSE,
                              missing = 1 - missing, x = x, weights = weights, 
                              estimand = "MO")
      coef1 <- c(result1$par)
      coef2 <- -c(result2$par)
      coef <- rbind(coef1, coef2)
      converged <- 
        as.logical((1 - result1$convergence) * (1 - result2$convergence))
      ps1 <- c(logistic(x %*% coef1))
      ps2 <- c(logistic(x %*% coef2))
      ps <- rbind(ps1, ps2)
      est.weights1 <- 
        invprob(ps = ps1, missing = missing, estimand = "MO", weights = weights)
      est.weights2 <- 
        invprob(ps = 1 - ps2, missing = 1 - missing, 
                estimand = "MO", weights = weights)
      est.weights <- est.weights1 + est.weights2
      est <- sum(est.weights * outcome * (2 * missing - 1))
    } else { # MO, ATT and ATEcombined
      initval <- 
        init(result = result.vanilla, method = method, estimand = estimand, 
             missing = missing, x = x, alpha = 0, N = N)
      result <- stats::optim(par = initval, 
                             fn = cbll, 
                             gr = cbllgradient,
                             method = "BFGS", 
                             control = list(fnscale = -1, trace = FALSE), 
                             hessian = FALSE,
                             missing = missing, x = x, weights = weights, 
                             estimand = estimand)
      coef <- c(result$par)
      converged <- as.logical(1 - result$convergence)
      ps <- c(logistic(x %*% coef))
      est.weights <- invprob(ps = ps, missing = missing, 
                             estimand = estimand, weights = weights)
      if (estimand == "MO") {
        est <- sum(est.weights * outcome * (1 - missing))
      } else { # ATT and ATEcombined
        est <- sum(est.weights * outcome * (2 * missing - 1))
      }
    }
    if (varcov == TRUE) {
      varcov <- cbV(missing = missing, ps = ps, x = x0, outcome = outcome, 
                    est = est, N = N, N1 = N1, weights = weights, 
                    estimand = estimand, est.weights = est.weights)
    }
  } else if (method == "score") {
    if (estimand == "ATE") {
      initval1 <- init(result = result.vanilla, method = method, estimand = "MO", 
                       missing = missing, x = x, alpha = 0, N = N)
      result.vanilla2 <- stats::glm(formula = I(1 - missing) ~ -1 + x, 
                                    family = stats::quasibinomial, 
                                    weights = weights)
      initval2 <- init(result = result.vanilla2, method = method, estimand = "MO", 
                       missing = 1 - missing, x = x, alpha = 0, N = N)
      result1 <- stats::optim(par = initval1, 
                              fn = navigatedll, 
                              gr = navigatedgradient,
                              method = "BFGS", 
                              control = list(fnscale = -1, trace = FALSE), 
                              hessian = FALSE,
                              missing = missing, x = x, weights = weights, 
                              estimand = "MO", alpha = alpha)
      result2 <- stats::optim(par = initval2, 
                              fn = navigatedll, 
                              gr = navigatedgradient,
                              method = "BFGS", 
                              control = list(fnscale = -1, trace = FALSE), 
                              hessian = FALSE,
                              missing = 1 - missing, x = x, weights = weights, 
                              estimand = "MO", alpha = alpha)
      coef1 <- c(result1$par)
      coef2 <- -c(result2$par)
      coef <- rbind(coef1, coef2)
      converged <- 
        as.logical((1 - result1$convergence) * (1 - result2$convergence))
      ps1 <- c(logistic(x %*% coef1))
      ps2 <- c(logistic(x %*% coef2))
      ps <- rbind(ps1, ps2)
      est.weights1 <- 
        invprob(ps = ps1, missing = missing, estimand = "MO", weights = weights)
      est.weights2 <- 
        invprob(ps = 1 - ps2, missing = 1 - missing, 
                estimand = "MO", weights = weights)
      est.weights <- est.weights1 + est.weights2
      est <- sum(est.weights * outcome * (2 * missing - 1))
    } else { # MO, ATT and ATEcombined
      initval <- 
        init(result = result.vanilla, method = method, estimand = estimand, 
             missing = missing, x = x, alpha = alpha, N = N)
      result <- stats::optim(par = initval, 
                             fn = navigatedll, 
                             gr = navigatedgradient,
                             method = "BFGS", 
                             control = list(fnscale = -1, trace = FALSE), 
                             hessian = FALSE,
                             missing = missing, x = x, weights = weights, 
                             estimand = estimand, alpha = alpha)
      coef <- c(result$par)
      converged <- as.logical(1 - result$convergence)
      ps <- c(logistic(x %*% coef))
      est.weights <- invprob(ps = ps, missing = missing, 
                             estimand = estimand, weights = weights)
      if (estimand == "MO") {
        est <- sum(est.weights * outcome * (1 - missing))
      } else { # ATT or ATEcombined
        est <- sum(est.weights * outcome * (2 * missing - 1))
      }
    }
    if (varcov == TRUE) {
      varcov <- navigatedV(missing = missing, ps = ps, x = x0, outcome = outcome, 
                           est = est, N = N, N1 = N1, weights = weights, 
                           estimand = estimand, est.weights = est.weights, 
                           alpha = alpha)
    }
  } else { # method = "both"
    if (twostep == FALSE) {
      if (estimand == "ATE") {
        initval1 <- init(result = result.vanilla, method = method, estimand = "MO", 
                         missing = missing, x = x, alpha = alpha, N = N)
        result.vanilla2 <- stats::glm(formula = I(1 - missing) ~ -1 + x, 
                                      family = stats::quasibinomial, 
                                      weights = weights)
        initval2 <- init(result = result.vanilla2, method = method, estimand = "MO",
                         missing = 1 - missing, x = x, alpha = alpha, N = N)
        result1 <- stats::optim(par = initval1, 
                                fn = gmm, 
                                method = "BFGS", 
                                control = list(fnscale = 1, trace = FALSE), 
                                hessian = FALSE,
                                missing = missing, x = x, 
                                N = N, N1 = N1, weights = weights, 
                                estimand = "MO", alpha = alpha)
        result2 <- stats::optim(par = initval2, 
                                fn = gmm, 
                                method = "BFGS", 
                                control = list(fnscale = 1, trace = FALSE), 
                                hessian = FALSE,
                                missing = 1 - missing, x = x, 
                                N = N, N1 = N - N1, weights = weights, 
                                estimand = "MO", alpha = alpha)
        coef1 <- c(result1$par)
        coef2 <- -c(result2$par)
        coef <- rbind(coef1, coef2)
        converged <- 
          as.logical((1 - result1$convergence) * (1 - result2$convergence))
        ps1 <- c(logistic(x %*% coef1))
        ps2 <- c(logistic(x %*% coef2))
        ps <- rbind(ps1, ps2)
        est.weights1 <- invprob(ps = ps1, missing = missing, 
                                estimand = "MO", weights = weights)
        est.weights2 <- invprob(ps = 1 - ps2, missing = 1 - missing, 
                                estimand = "MO", weights = weights)
        est.weights <- est.weights1 + est.weights2
        est <- sum(est.weights * outcome * (2 * missing - 1))
      } else { # MO, ATT and ATEcombined
        initval <- 
          init(result = result.vanilla, method = method, estimand = estimand, 
               missing = missing, x = x, alpha = alpha, N = N)
        result <- stats::optim(par = initval, 
                               fn = gmm, 
                               method = "BFGS", 
                               control = list(fnscale = 1, trace = FALSE), 
                               hessian = FALSE,
                               missing = missing, x = x, 
                               N = N, N1 = N1, weights = weights, 
                               estimand = estimand, alpha = alpha)
        coef <- c(result$par)
        converged <- as.logical(1 - result$convergence)
        ps <- c(logistic(x %*% coef))
        est.weights <- invprob(ps = ps, missing = missing, 
                               estimand = estimand, weights = weights)
        if (estimand == "MO") {
          est <- sum(est.weights * outcome * (1 - missing))
        } else { # ATT and ATEcombined
          est <- sum(est.weights * outcome * (2 * missing - 1))
        }
      }
      if (varcov == TRUE) {
        varcov <- gmmV(missing = missing, ps = ps, x = x0, outcome = outcome, 
                       est = est, N = N, N1 = N1, weights = weights, 
                       estimand = estimand, est.weights = est.weights, 
                       alpha = alpha)
      }
    } else { # twostep == TRUE
      if (estimand == "ATE") {
        initvalscore1 <- 
          init(result = result.vanilla, method = "score", estimand = "MO", 
               missing = missing, x = x, alpha = alpha, N = N)
        resultscore1 <- stats::optim(par = initvalscore1, 
                                     fn = navigatedll, 
                                     gr = navigatedgradient,
                                     method = "BFGS", 
                                     control = list(fnscale = -1, trace = FALSE), 
                                     hessian = FALSE,
                                     missing = missing, x = x, weights = weights, 
                                     estimand = "MO", alpha = alpha)
        initvalcb1 <- init(result = result.vanilla, method = "cb", estimand = "MO", 
                           missing = missing, x = x, alpha = 0, N = N)
        resultcb1 <- stats::optim(par = initvalcb1, 
                                  fn = cbll, 
                                  gr = cbllgradient,
                                  method = "BFGS", 
                                  control = list(fnscale = -1, trace = FALSE), 
                                  hessian = FALSE,
                                  missing = missing, x = x, weights = weights, 
                                  estimand = "MO")
        result.vanilla2 <- stats::glm(formula = I(1 - missing) ~ -1 + x, 
                                      family = stats::quasibinomial, 
                                      weights = weights)
        initvalscore2 <- 
          init(result = result.vanilla2, method = "score", estimand = "MO", 
               missing = missing, x = x, alpha = alpha, N = N)
        resultscore2 <- stats::optim(par = initvalscore2, 
                                     fn = navigatedll, 
                                     gr = navigatedgradient,
                                     method = "BFGS", 
                                     control = list(fnscale = -1, trace = FALSE), 
                                     hessian = FALSE,
                                     missing = 1 - missing, x = x, weights = weights, 
                                     estimand = "MO", alpha = alpha)
        initvalcb2 <- init(result = result.vanilla2, method = "cb", estimand = "MO", 
                           missing = missing, x = x, alpha = 0, N = N)
        resultcb2 <- stats::optim(par = initvalcb2, 
                                  fn = cbll, 
                                  gr = cbllgradient,
                                  method = "BFGS", 
                                  control = list(fnscale = -1, trace = FALSE), 
                                  hessian = FALSE,
                                  missing = 1 - missing, x = x, weights = weights, 
                                  estimand = "MO")
        ps1s <- c(logistic(x %*% c(resultscore1$par)))
        ps1c <- c(logistic(x %*% c(resultcb1$par)))
        ps2s <- c(logistic(x %*% c(-resultscore2$par)))
        ps2c <- c(logistic(x %*% c(-resultcb2$par)))
        varcovscore1 <- navigatedV(missing = missing, ps = ps1s, x = x0, 
                                   outcome = outcome, est = 0, N = N, N1 = N1, 
                                   weights = weights, estimand = "MO", 
                                   est.weights = rep(0, N), alpha = alpha)
        varcovcb1 <- cbV(missing = missing, ps = ps1c, x = x0, outcome = outcome, 
                         est = 0, N = N, N1 = N1, weights = weights, 
                         estimand = "MO", est.weights = rep(0 , N))
        varcovscore2 <- navigatedV(missing = 1 - missing, ps = 1 - ps2s, x = x0, 
                                   outcome = outcome, est = 0, N = N, N1 = N1, 
                                   weights = weights, estimand = "MO", 
                                   est.weights = rep(0, N), alpha = alpha)
        varcovcb2 <- cbV(missing = 1 - missing, ps = 1 - ps2c, x = x0, 
                         outcome = outcome, est = 0, N = N, N1 = N1, 
                         weights = weights, 
                         estimand = "MO", est.weights = rep(0 , N))
        scoreinvv1 <- 1 / apply(x0, 2, stats::sd) %*% sqrt(diag(varcovscore1)[1:ncol(x)])
        cbinvv1 <- 1 / apply(x0, 2, stats::sd) %*% sqrt(diag(varcovcb1)[1:ncol(x)])
        scoreinvv2 <- 1 / apply(x0, 2, stats::sd) %*% sqrt(diag(varcovscore2)[1:ncol(x)])
        cbinvv2 <- 1 / apply(x0, 2, stats::sd) %*% sqrt(diag(varcovcb2)[1:ncol(x)])
        scratio1 <- c(scoreinvv1 / (scoreinvv1 + cbinvv1))
        scratio2 <- c(scoreinvv2 / (scoreinvv2 + cbinvv2))
        scratio <- c(scratio1, scratio2)
        initval1 <- c(resultscore1$par) * scratio1 + 
                      c(resultcb1$par) * (1 - scratio1)
        initval2 <- c(resultscore2$par) * scratio2 + 
                      c(resultcb2$par) * (1 - scratio2)
        result1 <- stats::optim(par = initval1, 
                                fn = bothll,
                                gr = bothllgradient, 
                                method = "BFGS", 
                                control = list(fnscale = -1, trace = FALSE), 
                                hessian = FALSE,
                                missing = missing, x = x, weights = weights, 
                                estimand = "MO", alpha = alpha, scratio = scratio1)
        result2 <- stats::optim(par = initval2, 
                                fn = bothll,
                                gr = bothllgradient, 
                                method = "BFGS", 
                                control = list(fnscale = -1, trace = FALSE), 
                                hessian = FALSE,
                                missing = 1 - missing, x = x, weights = weights, 
                                estimand = "MO", alpha = alpha, scratio = scratio2)
        coef1 <- c(result1$par)
        coef2 <- -c(result2$par)
        coef <- rbind(coef1, coef2)
        converged <- 
          as.logical((1 - result1$convergence) * (1 - result2$convergence))
        ps1 <- c(logistic(x %*% coef1))
        ps2 <- c(logistic(x %*% coef2))
        ps <- rbind(ps1, ps2)
        est.weights1 <- invprob(ps = ps1, missing = missing, 
                                estimand = "MO", weights = weights)
        est.weights2 <- invprob(ps = 1 - ps2, missing = 1 - missing, 
                                estimand = "MO", weights = weights)
        est.weights <- est.weights1 + est.weights2
        est <- sum(est.weights * outcome * (2 * missing - 1))
      } else { # MO, ATT and ATEcombined
        initvalscore <- 
          init(result = result.vanilla, method = "score", estimand = estimand, 
               missing = missing, x = x, alpha = alpha, N = N)
        resultscore <- stats::optim(par = initvalscore, 
                                    fn = navigatedll, 
                                    gr = navigatedgradient,
                                    method = "BFGS", 
                                    control = list(fnscale = -1, trace = FALSE), 
                                    hessian = TRUE,
                                    missing = missing, x = x, weights = weights, 
                                    estimand = estimand, alpha = alpha)
        initvalcb <- 
          init(result = result.vanilla, method = "cb", estimand = estimand, 
               missing = missing, x = x, alpha = 0, N = N)
        resultcb <- stats::optim(par = initvalcb, 
                                 fn = cbll, 
                                 gr = cbllgradient,
                                 method = "BFGS", 
                                 control = list(fnscale = -1, trace = FALSE), 
                                 hessian = TRUE,
                                 missing = missing, x = x, weights = weights, 
                                 estimand = estimand)
        pssc <- c(logistic(x %*% c(resultscore$par)))
        pscb <- c(logistic(x %*% c(resultcb$par)))
        varcovscore <- navigatedV(missing = missing, ps = pssc, x = x0, 
                                  outcome = outcome, est = 0, N = N, N1 = N1, 
                                  weights = weights, estimand = "MO", 
                                  est.weights = rep(0, N), alpha = alpha)
        varcovcb <- cbV(missing = missing, ps = pscb, x = x0, outcome = outcome, 
                        est = 0, N = N, N1 = N1, weights = weights, 
                        estimand = "MO", est.weights = rep(0 , N))
        scoreinvv <- 1 / apply(x0, 2, stats::sd) %*% sqrt(diag(varcovscore)[1:ncol(x)])
        cbinvv <- 1 / apply(x0, 2, stats::sd) %*% sqrt(diag(varcovcb)[1:ncol(x)])
        scratio <- c(scoreinvv / (scoreinvv + cbinvv))
        initval <- c(resultscore$par) * scratio + c(resultcb$par) * (1 - scratio)
        result <- stats::optim(par = initval, 
                               fn = bothll,
                               gr = bothllgradient, 
                               method = "BFGS", 
                               control = list(fnscale = -1, trace = FALSE), 
                               hessian = FALSE,
                               missing = missing, x = x, weights = weights, 
                               estimand = estimand, alpha = alpha, scratio = scratio)
        coef <- c(result$par)
        converged <- as.logical(1 - result$convergence)
        ps <- c(logistic(x %*% coef))
        est.weights <- invprob(ps = ps, missing = missing, 
                               estimand = estimand, weights = weights)
        if (estimand == "MO") {
          est <- sum(est.weights * outcome * (1 - missing))
        } else { # ATT and ATEcombined
          est <- sum(est.weights * outcome * (2 * missing - 1))
        }
      }
      if (varcov == TRUE) {
        varcov <- bothV(missing = missing, ps = ps, x = x0, outcome = outcome, 
                        est = est, N = N, N1 = N1, weights = weights, 
                        estimand = estimand, est.weights = est.weights, 
                        alpha = alpha, scratio = scratio)
      }
    }
  }
  naive_coef <- c(result.vanilla$coefficients)
  ps.naive <- c(logistic(x %*% naive_coef))
  naive_weights <- invprob(ps = ps.naive, missing = missing, 
                           estimand = estimand, weights = weights)
  if (estimand != "ATE") {
    coef <- svdtranscoef(coef = coef, svdx = svdx)
    names(coef) <- names.x
    if (varcov[1] != FALSE) {
      colnames(varcov) <- c(names.x, "est")
      rownames(varcov) <- c(names.x, "est")
    }
  } else { # ATE
    coef1 <- svdtranscoef(coef = coef1, svdx = svdx)
    coef2 <- svdtranscoef(coef = coef2, svdx = svdx)
    names(coef1) <- names.x
    names(coef2) <- names.x
    coef <- rbind(coef1, coef2)
    if (varcov[1] != FALSE) {
      colnames(varcov) <- 
        c(paste0(rep(c("ps1_", "ps2_"), each = length(names.x)), 
                 rep(names.x, 2)), "est")
      rownames(varcov) <- 
        c(paste0(rep(c("ps1_", "ps2_"), each = length(names.x)), 
                 rep(names.x, 2)), "est")
    }
  }
  naive_coef <- svdtranscoef(coef = naive_coef, svdx = svdx)
  names(naive_coef) <- names.x
  if (estimand0 == "ATC") {
    estimand <- "ATC"
    est <- -est
    ps <- 1 - ps
    coef <- -coef
    naive_coef <- -naive_coef
    missing <- 1 - missing
  }
  list(est = est,
       weights = est.weights,
       ps = ps,
       coefficients = coef,
       varcov = varcov,
       converged = converged,
       naive_weights = naive_weights,
       naive_coef = naive_coef,
       scratio = scratio,
       estimand = estimand,
       method = method,
       outcome = outcome,
       alpha = alpha,
       names.x = names.x,
       prior.weights = weights,
       treat = c(missing))
}
