## covariate balance check
## standardized mean differences for continuous variables
## (unstandardized) proportion differences for categorical variables
cbcheck <- function (result, addcov = NULL, standardize = TRUE, 
										 plot = TRUE, absolute = TRUE, threshold = 0, sort = TRUE) {
	formula <- as.formula(result$formula)
	model <- model.frame(formula, data = result$data)
	missing <- c(model.extract(model, "response"))
	x <- as.matrix(model.matrix(formula, model))
	w.original <- result$prior.weights
	if (is.null(addcov) == 0) {
		formula2 <- as.formula(addcov)
		model2 <- model.frame(formula2, data = result$data)
		x2 <- as.matrix(model.matrix(formula2, model2))
		colnamesx <- colnames(x)
		colnamesx2 <- colnames(x2)[-1]
		x <- as.matrix(cbind(x, x2[, -1]))
		colnames(x) <- c(colnamesx, colnamesx2)
	}
	type <- apply(x, 2, function(cov) bincont(cov))
	type <- factor(type, levels = c("continuous", "binary"))
	if (result$estimand == "MO") {
		diff.adj <- apply(w.original * x, 2, sum) / sum(w.original) - 
									apply(result$weights * x * (1 - missing), 2, sum)
		diff.un <- apply(w.original * x, 2, sum) / sum(w.original) - 
								apply(w.original * x * (1 - missing), 2, sum) / sum(w.original * (1 - missing))
	} else { # ATT, ATE, and ATEcombined
		diff.adj <- apply(result$weights * x * (2 * missing - 1), 2, sum)
		diff.un <- apply(w.original * x * missing, 2, sum) / sum(w.original * missing) - 
								apply(w.original * x * (1 - missing), 2, sum) / sum(w.original * (1 - missing))
	}
	covariates <- factor(colnames(x), levels = colnames(x))
	res <- data.frame(covariates = covariates, 
										type = type,
										diff.adj = diff.adj,
										diff.un = diff.un)
	if (standardize == TRUE) {
		if (result$estimand == "ATT") {
			std <- apply(x[missing == 1, ], 2, function(cov) weighted_sd(cov, weights = w.original[missing == 1]))
		} else { # ATE and MO
			std <- apply(x, 2, function(cov) weighted_sd(cov, weights = w.original))
		}
		std[type == "binary"] <- 1
		std[1] <- 1
		res <- res %>% mutate(diff.adj = diff.adj / std,
													diff.un = diff.un / std)
	}
	rownames(res) <- NULL
	plotcb <- plot_balance(result = res, standardize = standardize, 
												 absolute = absolute, threshold = threshold, sort = sort)
	if (plot == TRUE) {
		print(plotcb)
	}
	res <- list(balance = res,
							plot = plotcb)
	invisible(res)
}
