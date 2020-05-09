## Print (S3 method)
print.nawt <- function (object, ...) {
	est <- signif(object$est, digits = getOption("digits") - 3)
	if (object$estimand == "ATE") {
		coef <- c(object$coef[1, ], object$coef[2, ])
		names(coef) <- c(paste0(rep(c("ps1_", "ps2_"), each = length(coef) / 2), 
										 				object$names.x))
	} else { # ATT, ATEcombined, and MO
		coef <- object$coef
	}
	cat("\nCall: \n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
			"\n", sep = "")
	cat("\n", paste(object$estimand, "estimates: ", est), "\n", sep = "")
	cat("\nCoefficients:\n")
	print(coef, digits = getOption("digits") - 3)
	cat("\nEffective N for propensity score estimation:", 
			round(object$effN_ps, digits = 2))
	if (object$estimand == "MO") {
		cat("\nEffective N for the", object$estimand, "estimation: ", 
				round(object$effN_est, digits = 2), "\n")
	} else {
		cat("\nEffective N for the", object$estimand, "estimation:\n")
		cat("  treatment:", round(object$effN_est[1], digits = 2), 
				"\n  control:  ", round(object$effN_est[2], digits = 2), "\n")
	}
}
