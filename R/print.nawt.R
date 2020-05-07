## Print (S3 method)
print.nawt <- function (result) {
	est <- signif(result$est, digits = getOption("digits") - 3)
	if (result$estimand == "ATE") {
		coef <- c(result$coef[1, ], result$coef[2, ])
		names(coef) <- c(paste0(rep(c("ps1_", "ps2_"), each = length(coef) / 2), 
										 				result$names.x))
	} else { # ATT, ATEcombined, and MO
		coef <- result$coef
	}
	cat("\nCall: \n", paste(deparse(result$call), sep = "\n", collapse = "\n"), 
			"\n", sep = "")
	cat("\n", paste(result$estimand, "estimates: ", est), "\n", sep = "")
	cat("\nCoefficients:\n")
	print(coef, digits = getOption("digits") - 3)
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
}
