## Comparison between weights: scattered plot
plot.nawt <- function (result) {
	maxweight <- max(c(result$naive_weights, result$weights))
	minweight <- min(c(result$naive_weights, result$weights))
	par(pty = "s")
	plot(result$naive_weights[result$treat == 0], 
			 result$weights[result$treat == 0], 
			 col = rgb(39/ 255, 139/ 255, 210 / 255, alpha = 0.75),
			 xlim = c(0, maxweight), ylim = c(0, maxweight),
			 xlab = "", ylab = "")
	par(new = TRUE)
	plot(result$naive_weights[result$treat == 1], 
			 result$weights[result$treat == 1], 
			 col = rgb(220 / 255, 50 / 255, 46 / 255, alpha = 0.75),
			 xlim = c(0, maxweight), ylim = c(0, maxweight),
			 xlab = "Naive weights", ylab = "Navigated weights")
	text(x = minweight, y = maxweight * 0.9, labels = "control", 
			 col = rgb(39/ 255, 139/ 255, 210 / 255, alpha = 0.85), 
			 adj = 0, cex = 1.2)
	text(x = minweight, y = maxweight * 0.95, labels = "treat", 
			 col = rgb(220 / 255, 50 / 255, 46 / 255, alpha = 0.85), 
			 adj = 0, cex = 1.2)
	abline(0, 1)
}
