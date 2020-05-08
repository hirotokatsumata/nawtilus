## Classify covariates into binary or continuous (internal only)
bincont <- function (cov) {
	ifelse(length(unique(cov)) == 2, "binary", "continuous")
}

## Weighted standard deviation (internal only)
weighted_sd <- function (cov, weights) {
	weighted.mean <- sum(cov * weights) / sum(weights)
	sqrt(1 / (sum(weights) - 1) * sum(weights * (cov - weighted.mean)^2))
}

## Covariate balance plot (internal only)
plot_balance <- function (result, standardize = TRUE, absolute = TRUE, threshold = 0, sort = TRUE) {
	if (standardize == TRUE) {
		std <- "Standardized"
	} else {
		std <- "Unstandardized"
	}
	diff <- paste(std, "Mean Differences")
	result$covariates <- factor(result$covariates, levels = result$covariates[nrow(result):1])
	if (absolute == TRUE) {
		result$diff.adj <- abs(result$diff.adj)
		result$diff.un <- abs(result$diff.un)
		diff <- paste("Absolute", std, "Mean Differences")
	} else { # absolute == FALSE
		threshold <- c(-threshold, threshold)
	}
	if (sort == TRUE) {
		order.un <- order(result$diff.un)
		if (absolute == FALSE) {
			order.un <- order(result$diff.un, decreasing = TRUE)			
		}
		result$covariates <- factor(result$covariates, levels = result$covariates[order.un])
	}
	ggplot(result, aes(x = diff.un, y = covariates, colour = "unadjusted", shape = type)) +
		geom_point(size = 4.5, stroke = 2.2) +
		geom_point(aes(x = diff.adj, y = covariates, colour = "adjusted", shape = type), 
									 size = 4.5, stroke = 2.2) +
		scale_shape_manual(values=c(1, 2)) +
		geom_vline(xintercept = threshold, colour = "grey50", linetype = "dashed", size = 0.6) +
		geom_vline(xintercept = 0, colour = "grey10", linetype = "solid", size = 0.6) +
		xlab(diff) +
		ylab("Covariates") +
		theme_classic() +
		theme(axis.title.x = element_text(size = 20), 
				  axis.title.y = element_text(size = 20),
				  axis.text.x = element_text(size = 22),
				  axis.text.y = element_text(size = 22),
				  legend.title = element_text(size = 0),
				  legend.text = element_text(size = 24),
				  plot.title = element_text(hjust = 0.5, size = 26)) +
		ggtitle("Covariate balance")
}
