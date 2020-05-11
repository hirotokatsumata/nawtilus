#' Summarize and Plot Covariate Balance
#'
#' Summarizes and plots covariate balance between treatment and control groups 
#'   before and after the navigated weighting.
#'
#' @export
#'
#' @param object an object of class “nawt”, usually, a result of a call to \code{\link{nawt}}.
#' @param addcov a one-sided formula specifying additional covariates whose 
#'   balance will be checked.
#' @param standardize a logical value indicating whether weighted mean 
#'   differences are standardiezed or not.
#' @param plot a logical value indicating whether a cavariate balance plot is
#'   displyaed. 
#' @param absolute a logical value indicating whether the absolute values of 
#'   differeces in weighted means are used in the covariate balance plot.
#' @param threshold an optional numeric vector used as threshold markers in the 
#'   covariate balance plot.
#' @param sort a logical value indicating whether covariates in the covariate 
#'   balance plot are sorted by the values of differences in the weighted means 
#'   before the navigated weighting.
#'
#' @return A matrix whose rows are the covariates and columns are the 
#'   differences in the (un)standardized weighted mean between the treatment and
#'   control groups before and after the navigated weighting. The standardized 
#'   weighted mean is the weighted mean divided by the standard deviation of the 
#'   covariate for the target population (the treatment group for the average 
#'   treatment effects on the treated estimation and the whole population for 
#'   the other quantity of interest. The differences in the categorical 
#'   variables are not standardized.
#'
#' @author Hiroto Katsumata.
#' 
#' @examples 
#' # Simulation from Kang and Shafer (2007) and Imai and Ratkovic (2014)
#' # ATT estimation
#' # True ATT is 10
#' tau <- 10
#' set.seed(12345)
#' n <- 1000
#' X <- matrix(rnorm(n * 4, mean = 0, sd = 1), nrow = n, ncol = 4)
#' prop <- 1 / (1 + exp(X[, 1] - 0.5 * X[, 2] + 0.25 * X[, 3] + 0.1 * X[, 4]))
#' treat <- rbinom(n, 1, 1 - prop)
#' y <- 210 + 27.4 * X[, 1] + 13.7 * X[, 2] + 13.7 * X[, 3] + 13.7 * X[, 4] + 
#'      tau * treat + rnorm(n)
#' df <- data.frame(X, treat, y)
#' colnames(df) <- c("x1", "x2", "x3", "x4", "treat", "y")
#'
#' # A misspecified model
#' Xmis <- data.frame(x1mis = exp(X[, 1] / 2), 
#'                    x2mis = X[, 2] * (1 + exp(X[, 1]))^(-1) + 10,
#'                    x3mis = (X[, 1] * X[, 3] / 25 + 0.6)^3, 
#'                    x4mis = (X[, 2] + X[, 4] + 20)^2)
#'
#' # Data frame and a misspecified formula for propensity score estimation
#' df <- data.frame(df, Xmis)
#' formula_m <- as.formula(treat ~ x1mis + x2mis + x3mis + x4mis)
#'
#' # Misspecified propensity score model
#' # Power weighting function with alpha = 2
#' fits2m <- nawt(formula = formula_m, outcome = "y", estimand = "ATT", 
#'                method = "score", data = df, alpha = 2)
#' cbcheck(fits2m, addcov = ~ x1 + x2 + x3 + x4)
#' 
#' # Covariate balancing weighting function
#' fitcbm <- nawt(formula = formula_m, outcome = "y", estimand = "ATT", 
#'                method = "cb", data = df)
#' cbcheck(fitcbm, addcov = ~ x1 + x2 + x3 + x4)
#' 
#' # Standard logistic regression
#' fits0m <- nawt(formula = formula_m, outcome = "y", estimand = "ATT", 
#'                method = "score", data = df, alpha = 0)
#' cbcheck(fits0m, addcov = ~ x1 + x2 + x3 + x4)
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
