## DGP for simulation studies
dgp <- function (coef.ps, coef.y, N, ybinary = FALSE, 
								 formula, formula.ps) {
	covariate <- rnorm(n = N * 4, mean = 0, sd = 1)
	data <- data.frame(id = 1:N,
					  				 x1 = covariate[1:N], 
					  				 x2 = covariate[(N + 1):(N * 2)],
					  				 x3 = covariate[(N * 2 + 1):(N * 3)],
					  				 x4 = covariate[(N * 3 + 1):(N * 4)]) %>%
					mutate(x001 = exp(x1 / 2),
							   x002 = x2 / (1 + exp(x1)) + 10,
							   x003 = (x1 * x3 / 25 + 0.6)^3,
							   x004 = (x2 + x4 + 20)^2) %>%
					mutate(x101 = (x001 - mean(x001)) / sd(x001),
							   x102 = (x002 - mean(x002)) / sd(x002),
							   x103 = (x003 - mean(x003)) / sd(x003),
							   x104 = (x004 - mean(x004)) / sd(x004))
	assign.ps <- dgp0(formula = formula.ps, 
								  	coef = coef.ps,
								  	data = data,
								  	binary = TRUE)
	data$treat <- assign.ps$outcome
	data$tps <- assign.ps$ps
	data$outcome <- dgp0(formula = formula, 
											 coef = coef.y,
					  					 data = data,
					  					 binary = ybinary)$outcome
	data
}

## Data generating functions
dgp0 <- function (formula, coef, data, binary = FALSE) {
	N <- nrow(data)
	attach(data)
	X <- as.matrix(model.matrix(formula))
	detach()
	y_star <- c(X %*% coef)
	if (binary == TRUE) {
		y_star <- logistic(y_star)
		outcome <- rbinom(rep(1, N), rep(1, N), prob = y_star)
	} else {
		outcome <- y_star + rnorm(N, mean = 0, sd = 1)
	}
	list(outcome = outcome, ps = y_star)
}
