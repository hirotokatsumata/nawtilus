## covariate balance check
## standardized mean differences for continuous variables
## (unstandardized) proportion differences for categorical variables
cbcheck <- function (object, addcov = NULL, standardize = TRUE, 
										 plot = TRUE, absolute = TRUE, threshold = 0, sort = TRUE) {
	formula <- as.formula(object$formula)
	model <- model.frame(formula, data = object$data)
	missing <- c(model.extract(model, "response"))
	x <- as.matrix(model.matrix(formula, model))
	w.original <- object$prior.weights
	if (is.null(addcov) == 0) {
		formula2 <- as.formula(addcov)
		model2 <- model.frame(formula2, data = object$data)
		x2 <- as.matrix(model.matrix(formula2, model2))
		colnamesx <- colnames(x)
		colnamesx2 <- colnames(x2)[-1]
		x <- as.matrix(cbind(x, x2[, -1]))
		colnames(x) <- c(colnamesx, colnamesx2)
	}
	type <- apply(x, 2, function(cov) bincont(cov))
	type <- factor(type, levels = c("continuous", "binary"))
	if (object$estimand == "MO") {
		diff.adj <- apply(w.original * x, 2, sum) / sum(w.original) - 
									apply(object$weights * x * (1 - missing), 2, sum)
		diff.un <- apply(w.original * x, 2, sum) / sum(w.original) - 
								apply(w.original * x * (1 - missing), 2, sum) / 
												sum(w.original * (1 - missing))
	} else { # ATT, ATE, and ATEcombined
		diff.adj <- apply(object$weights * x * (2 * missing - 1), 2, sum)
		diff.un <- apply(w.original * x * missing, 2, sum) / 
											sum(w.original * missing) - 
								apply(w.original * x * (1 - missing), 2, sum) / 
												sum(w.original * (1 - missing))
	}
	covariates <- factor(colnames(x), levels = colnames(x))
	res <- data.frame(covariates = covariates, 
										type = type,
										diff.adj = diff.adj,
										diff.un = diff.un)
	if (standardize == TRUE) {
		if (object$estimand == "ATT") {
			std <- 
				apply(x[missing == 1, ], 2, 
							function(cov) weighted_sd(cov, 
									 											weights = w.original[missing == 1]))
		} else { # ATE and MO
			std <- apply(x, 2, function(cov) weighted_sd(cov, weights = w.original))
		}
		std[type == "binary"] <- 1
		std[1] <- 1
		res$diff.adj <- res$diff.adj / std
		res$diff.un <- res$diff.un / std
	}
	rownames(res) <- NULL
	if (plot == TRUE) {
		plot_balance(result = res, standardize = standardize, 
								 absolute = absolute, threshold = threshold, sort = sort)
	}
	invisible(res)
}
