#' Summarizing Navigated Weighting Estimation
#'
#' Prints a summary of a fitted \code{nawt} object.
#'
#' Prints a summmary of a \code{nawt} object, in a format similar to glm.
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
## S3 method for class 'CBPS'
summary.nawt <- function (result) {
	est <- signif(result$est, digits = getOption("digits") - 3)
	k <- length(result$coef)
	se <- sqrt(diag(result$varcov))
	if (result$estimand == "ATE") {
		coef <- c(result$coef[1, ], result$coef[2, ])
	} else { # ATT, ATEcombined, and MO
		coef <- result$coef
	}
	coef.table <- cbind(as.vector(c(result$est, coef)),
											as.vector(se[c(k + 1, 1:k)]),
											as.vector(c(result$est, coef) / se[c(k + 1, 1:k)]),
											as.vector(2 - 2 * pnorm(abs(c(result$est, coef) / 
																										se[c(k + 1, 1:k)]))))
	colnames(coef.table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
	if (result$estimand == "ATE") {
		rownames(coef.table) <- c("est", result$names.x, result$names.x)
	} else { # ATT, ATEcombined, and MO
		rownames(coef.table) <- c("est", result$names.x)
	}
	pval <- coef.table[, 4]
	symp <- symnum(pval, corr = FALSE,
								 cutpoints = c(0, .001, .01, .05, .1, 1),
								 symbols = c("***", "**", "*", ".", " "))
	coef.print <- data.frame(coef.table, as.vector(symp))
	colnames(coef.print)[5] <- ""
	if (result$estimand == "ATE") {
		rownames(coef.print) <- 
			c("est", paste0(rep(c("ps1_", "ps2_"), each = k / 2), result$names.x))
	} else { # ATT, ATEcombined, and MO
		rownames(coef.print) <- c("est", result$names.x)
	}
	cat("\nCall: \n", paste(deparse(result$call), sep = "\n", collapse = "\n"), 
			"\n", sep = "")
	cat("\n", paste(result$estimand, "estimates: ", est), "\n", sep = "")
	cat("\nCoefficients:\n")
	print(coef.print, digits = getOption("digits") - 3)
	cat("---\n")
	cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
	cat("\nEffective N for propensity score estimation:", 
			round(result$effN_ps, digits = 2))
	if (result$estimand == "MO") {
		cat("\nEffective N for the", result$estimand, "estimation: ", 
				round(result$effN_est, digits = 2), "\n")
	} else {
		cat("\nEffective N for the", result$estimand, "estimation:\n")
		cat("  treatment:", round(result$effN_est[1], digits = 2), 
				"\n  control:  ", round(result$effN_est[2], digits = 2), "\n")
	}
	out <- list("call" = result$call, "coefficients" = coef.table)
	invisible(out)
}
