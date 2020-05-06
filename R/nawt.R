## library(MASS) for ginv()
## library(hypergeo) for hypergeo()
## weights: user-specified weights; e.g. sampling weights

## navigated weighitng: estimation
nawt <- function (formula, estimand = "ATT", method = "score", outcome, alpha = 2,
									twostep = TRUE, weights = NULL, 
									boot = FALSE, B = 2000, clevel = 0.95, 
									message = TRUE, data) {
	# Check
	if (estimand %in% c("ATE", "ATT", "MO", "ATEcombined") == 0) {
		stop("estimand must be \"ATE\", \"ATT\", \"MO\", or \"ATEcombined\"")
	}
	if (method %in% c("score", "cb", "both") == 0) {
		stop("method must be \"score\", \"cb\", or \"both\"")
	}
	stopifnot(length(alpha) == 1)
	if (alpha < 0) {
		stop("alpha must be equal to or larger than 0")
	}
	if (min(clevel) <= 0 | max(clevel) >= 1) {
		stop("clevel must be between 0 and 1")
	}
	call <- match.call()
	if (message == TRUE) {
		if (method == "score") {
			printmethod <- "estimation by (weighted) score conditions"
		} else if (method == "cb") {
			printmethod <- "estimation by covariate balancing"
		} else { # method == "both" 
			printmethod <- "estimation by (weighted) score and covariate balancing conditions"
		}
		print(paste("Estimate weights for the", estimand, printmethod))
	}
	if (boot == TRUE) {
		stopifnot(length(B) == 1)
		result <- nawt0(formula = formula, estimand = estimand, method = method, 
										outcome = outcome, alpha = alpha, twostep = twostep, weights = weights, 
										varcov = FALSE, boot = boot, data = data)
		result2 <- matrix(, nrow = B, ncol = length(result$coefficients) + 1)
		for (b in 1:B) {
			bs.sample <- sample(1:nrow(data), nrow(data), replace = TRUE)
			result0 <- nawt0(formula = formula, estimand = estimand, method = method, 
											 outcome = outcome, alpha = alpha, twostep = twostep, weights = weights[bs.sample], 
											 varcov = FALSE, boot = boot, data = data[bs.sample, ])
			if (estimand != "ATE") {
				result2[b, ] <- c(result0$coefficients, result0$est)
			} else { # ATE
				result2[b, ] <- c(result0$coefficients[1, ], result0$coefficients[2, ], result0$est)
			}
		}
		result$varcov <- cov(result2)
		lower <- sort(result2[, ncol(result2)])[floor(signif(B * (1 - clevel) / 2, digits = 5))]
		upper <- sort(result2[, ncol(result2)])[ceiling(signif(B * (1 + clevel) / 2, digits = 5)) + 1]
		result$ci <- cbind(lower, upper)
		if (estimand != "ATE") {
			colnames(result$varcov) <- c(result$names.x, "est")
			rownames(result$varcov) <- c(result$names.x, "est")
		} else { # ATE
			colnames(result$varcov) <- c(paste0(rep(c("ps1_", "ps2_"), each = length(result$names.x)), 
																					rep(result$names.x, 2)), "est")
			rownames(result$varcov) <- c(paste0(rep(c("ps1_", "ps2_"), each = length(result$names.x)), 
																					rep(result$names.x, 2)), "est")
		}
	} else {
		result <- nawt0(formula = formula, estimand = estimand, method = method,
										outcome = outcome, alpha = alpha, twostep = twostep, weights = weights, 
										varcov = TRUE, boot = boot, data = data)
		cilength <- sqrt(diag(result$varcov)[nrow(result$varcov)]) * qnorm((1 + clevel) / 2)
		result$ci <- cbind(lower = result$est - cilength, upper = result$est + cilength)
	}
	rownames(result$ci) <- paste0(clevel * 100, "%")
	result$omega <- weight_omega(result = result)
	if (estimand != "ATE") {
		result$effN_ps <- sum(result$prior.weights * result$omega)^2 / 
											sum((result$prior.weights * result$omega)^2)
	} else { # ATE
		result$effN_ps <- c(sum(result$prior.weights * result$omega[1, ])^2 / 
												sum(result$prior.weights * result$omega[1, ]^2),
											 sum(result$prior.weights * result$omega[2, ])^2 / 
												sum(result$prior.weights * result$omega[2, ]^2))
	}
	if (estimand != "MO") {
		result$effN_est <- c(treat = sum(result$weights[result$treat == 1])^2 / 
																	sum(result$weights[result$treat == 1]^2),
												control = sum(result$weights[result$treat == 0])^2 / 
																	sum(result$weights[result$treat == 0]^2))
	} else { # MO
		result$effN_est <- c(sum(result$weights[result$treat == 0])^2 / sum(result$weights[result$treat == 0]^2))
	}
	effN_original <- sum(result$prior.weights)^2 / sum((result$prior.weights)^2)
	if (method != "both" | twostep == TRUE) {
		if (min(result$effN_ps) < effN_original / 4) {
			warning(paste("Propensity score estimates may be unstable (Effective N = ", result$effN_ps,
										"). Try smaller alpha. "))
		}
	}
	result$effN_original <- effN_original
	result$call <- call
	class(result) <- "nawt"
	result
}
