## Comparison between weights: scattered plot
plot.nawt <- function (object, ...) {
	maxweight <- max(c(object$naive_weights, object$weights))
	minweight <- min(c(object$naive_weights, object$weights))
	par(pty = "s")
	plot(object$naive_weights[object$treat == 0], 
			 object$weights[object$treat == 0], 
			 col = rgb(39/ 255, 139/ 255, 210 / 255, alpha = 0.75),
			 xlim = c(0, maxweight), ylim = c(0, maxweight),
			 xlab = "", ylab = "")
	par(new = TRUE)
	plot(object$naive_weights[object$treat == 1], 
			 object$weights[object$treat == 1], 
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
