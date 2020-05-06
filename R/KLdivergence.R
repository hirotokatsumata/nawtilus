## Kullback-Leibler divergence
KLdivergence <- function (weights = est.weights, estimand, missing) {
	KLd <- sum(weights[missing == 0] * log(weights[missing == 0] / mean(weights[missing == 0])))
	if (estimand %in% c("ATE", "ATEcombined")) {
		KLd <- KLd + sum(weights[missing == 1] * log(weights[missing == 1] / mean(weights[missing == 1])))
	}
	KLd
}
