## Plot weights for the score omega(hat(pi)) and estimated propensity score distribution
plot_omega <- function (result, relative = TRUE) {
	omega <- result$omega
	if (relative == TRUE) {
		ylab <- expression(paste("Relative weights for the score ", omega(hat(pi))))
	} else {
		ylab <- expression(paste("Absolute weights for the score ", omega(hat(pi))))
	}
	if (result$estimand != "ATE") {
		dens <- density(result$ps, from = 0, to = 1)
		limx <- c(0, 1)
		plot(dens,
				 xlim = limx,
				 ylim = c(0, max(dens$y)),
				 lwd = 1.2,
				 ann = FALSE,
				 axes = FALSE)
	} else { # ATE
		dens1 <- density(result$ps[1, ], from = 0, to = 1)
		dens2 <- density(result$ps[2, ], from = 0, to = 1)
		limx <- c(0, 1)
		plot(dens1,
				 lty = 1,
				 col = rgb(0.5, 0.5, 0.5),
				 xlim = limx,
				 ylim = c(0, max(c(dens1$y, dens2$y))),
				 lwd = 1.2,
				 ann = FALSE,
				 axes = FALSE)
		par(new = TRUE)
		plot(dens2,
				 lty = 2,
				 col = rgb(0.2, 0.2, 0.2),
				 xlim = limx,
				 ylim = c(0, max(c(dens1$y, dens2$y))),
				 lwd = 1.2,
				 ann = FALSE,
				 axes = FALSE)
	}
	axis(side = 4)
	mtext("Density", side = 4, line = 3)
	par(new = TRUE)
	if (result$estimand == "ATE") {
		ps1 <- result$ps[1, ]
		ps2 <- result$ps[2, ]
		omega1 <- result$omega[1, ]
		omega2 <- result$omega[2, ]
		if (relative == TRUE) {
			omega1 <- omega1 / mean(omega1)
			omega2 <- omega2 / mean(omega2)
		}
		psc <- c(ps1[result$treat == 0], ps2[result$treat == 0])
		pst <- c(ps1[result$treat == 1], ps2[result$treat == 1])
		omegac <- c(omega1[result$treat == 0], omega2[result$treat == 0])
		omegat <- c(omega1[result$treat == 1], omega2[result$treat == 1])
		plot(psc, omegac,
				 col = rgb(39/ 255, 139/ 255, 210 / 255, alpha = 0.3),
				 xlim = limx,
				 ylim = c(0, max(c(omega1, omega2))),
				 xlab = "", 
				 ylab = "",
				 axes = FALSE)
		par(new = TRUE)
		plot(pst, omegat,
				 col = rgb(220 / 255, 50 / 255, 46 / 255, alpha = 0.3),
				 xlim = limx,
				 ylim = c(0, max(c(omega1, omega2))),
				 xlab = expression(paste("Estimated propensity score ", hat(pi))), 
				 ylab = ylab,
				 axes = FALSE)
	} else { # ATT and MO
		if (relative == TRUE) {
			omega <- omega / mean(omega)
		}
		omegac <- omega[result$treat == 0]
		omegat <- omega[result$treat == 1]
		plot(result$ps[result$treat == 0], omegac,
				 col = rgb(39/ 255, 139/ 255, 210 / 255, alpha = 0.3),
				 xlim = limx,
				 ylim = c(0, max(c(omegac, omegat))),
				 xlab = "", 
				 ylab = "",
				 axes = FALSE)
		par(new = TRUE)
		plot(result$ps[result$treat == 1], omegat,
				 col = rgb(220 / 255, 50 / 255, 46 / 255, alpha = 0.3),
				 xlim = limx,
				 ylim = c(0, max(c(omegac, omegat))),
				 xlab = expression(paste("Estimated propensity score ", hat(pi))), 
				 ylab = ylab,
				 axes = FALSE)
	}
	axis(side = 1)
	axis(side = 2)
	box(bty = "u")
}
