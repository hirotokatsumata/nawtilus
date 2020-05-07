#' Navigated Weighting (NAWT) Estimation
#'
#' \code{nawt} estimates a pre-specified parameter of interest (e.g., the 
#' average treatment effects or the average treatment effects on the treated) 
#' with the inverse probability weighting where propensity scores are estimated
#' using estimating equations suitable for the parameter of interest. It 
#' includes the covariate balancing propensity score proposed by Imai and 
#' Ratkovic (2014), which uses covariate balancing conditions in propensity 
#' score estimation. \code{nawt} can also be used to estimate average outcomes
#' in missing outcome cases.
#'
#' DETAILS
#'
#' @param formula an object of class \code{\link[stats]{formula}} (or one that
#'   can be coerced to that class): a symbolic description of the model to be 
#'   fitted.
#' @param outcome a character string specifying the name of outcome values
#'   in \code{data}
#' @param estimand a character string specifying a parameter of interest. Choose
#'   "ATT" for the average treatment effect on treated estimation, "ATE" for the
#'   average treatment effect estimation, or "MO" for the average outcomes 
#'   estimation in missing outcome cases. You can choose "ATEcombined" for the  
#'   combined estimation for the average treatment effect estimation.
#' @param method a character string specifying a type of weignting functions in 
#'   propensity score estimation (\eqn{\omega(\pi)}). Choose "score" for a power
#'   function of propensity scores (need to specify the value for alpha), "cb"
#'   for a covariate balancing weighting function, or "both" to use both the
#'   above weighting functions (need to specify the value for alpha).
#' @param data a data frame containing the outcomes and the variables in the 
#'   model.
#' @param weights an optional vector of ‘prior weights’ (e.g. sampleing weights)
#'   to be used in the fitting process. Should be NULL or a numeric vector.
#' @param alpha a positive value for an exponent in a power weighting function.
#'   Default is 2. Set to 0 to use the standard logistic regression for  
#'   propensity score estimation. Note that \code{nawt} with \eqn{m / 2} where 
#'   \eqn{0 \le m \le 10} runs substantially faster than with any other values.
#' @param twostep a logical value indicating whether to use a two step estimator
#'   when \code{method = "both"}. Default is \code{TRUE}. Set to \code{FALSE} 
#'   to use a continuously-updating GMM estimator, which is substantially
#'   computationally intensive.
#' @param boot a logical value indicating whether to use a non-parametric 
#'   bootstrapping method to estimate the variance-covariance matrix and 
#'   confidence intervals for parameters. Default is \code{FALSE}. Set to 
#'   \code{FALSE} to use a sandwich-type asymptotic covariance estimator.
#' @param B the number of bootstrap replicates. Default is 2,000.
#' @param clevel confidence level. Default is 0.95.
#' @param message a logical value indicating whether messages are shown or not.
nawt <- function (formula, outcome, estimand = "ATT", method = "score", 
									data, weights = NULL, alpha = 2, twostep = TRUE, 
									boot = FALSE, B = 2000, clevel = 0.95, message = TRUE) {
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
			printmethod <- 
				"estimation by (weighted) score and covariate balancing conditions"
		}
		print(paste("Estimate weights for the", estimand, printmethod))
	}
	if (boot == TRUE) {
		stopifnot(length(B) == 1)
		result <- nawt0(formula = formula, outcome = outcome, estimand = estimand, 
										method = method, data = data, weights = weights, 
										alpha = alpha, twostep = twostep, varcov = FALSE)
		result2 <- matrix(, nrow = B, ncol = length(result$coefficients) + 1)
		for (b in 1:B) {
			bs.sample <- sample(1:nrow(data), nrow(data), replace = TRUE)
			result0 <- nawt0(formula = formula, outcome = outcome, estimand = estimand, 
											 method = method, data = data[bs.sample, ], 
											 weights = weights[bs.sample], alpha = alpha, 
											 twostep = twostep, varcov = FALSE)
			if (estimand != "ATE") {
				result2[b, ] <- c(result0$coefficients, result0$est)
			} else { # ATE
				result2[b, ] <- 
					c(result0$coefficients[1, ], result0$coefficients[2, ], result0$est)
			}
		}
		result$varcov <- cov(result2)
		lower <- sort(result2[, ncol(result2)])[
								floor(signif(B * (1 - clevel) / 2, digits = 5))]
		upper <- sort(result2[, ncol(result2)])[
								ceiling(signif(B * (1 + clevel) / 2, digits = 5)) + 1]
		result$ci <- cbind(lower, upper)
		if (estimand != "ATE") {
			colnames(result$varcov) <- c(result$names.x, "est")
			rownames(result$varcov) <- c(result$names.x, "est")
		} else { # ATE
			colnames(result$varcov) <- 
				c(paste0(rep(c("ps1_", "ps2_"), each = length(result$names.x)), 
								 rep(result$names.x, 2)), "est")
			rownames(result$varcov) <- 
				c(paste0(rep(c("ps1_", "ps2_"), each = length(result$names.x)), 
								 rep(result$names.x, 2)), "est")
		}
	} else {
		result <- nawt0(formula = formula, outcome = outcome, estimand = estimand, 
										method = method, data = data, weights = weights, 
										alpha = alpha, twostep = twostep, varcov = TRUE)
		cilength <- 
			sqrt(diag(result$varcov)[nrow(result$varcov)]) * qnorm((1 + clevel) / 2)
		result$ci <- 
			cbind(lower = result$est - cilength, upper = result$est + cilength)
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
		result$effN_est <- c(sum(result$weights[result$treat == 0])^2 / 
													sum(result$weights[result$treat == 0]^2))
	}
	effN_original <- sum(result$prior.weights)^2 / sum((result$prior.weights)^2)
	if (method != "both" | twostep == TRUE) {
		if (min(result$effN_ps) < effN_original / 4) {
			warning(paste("Propensity score estimates may be unstable (Effective N = ", 
										result$effN_ps, "). Try smaller alpha. "))
		}
	}
	result$effN_original <- effN_original
	result$call <- call
	class(result) <- "nawt"
	result
}
