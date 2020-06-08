#' Plot weights for propensity score estimation in the navigated weighting
#'
#' Plots weight of each observation in the score condition \eqn{\omega(\pi)} for
#'   propensity score estimation and estimated propensity score distribution 
#'   in the navigated weighting.
#'
#' The x-axis shows estimated propensity scores, and the y-axis shows weight of 
#' each observation in propensity score estimation. When \code{estimand = "ATE"},
#' the navigated weighting estimates two propensity scores for each observation;
#' one for estimating the average of the potential outcomes with treatment and  
#' the other for estimating the average of the potential outcomes without 
#' treatment. Therefore, there are two weighting functions for estimating two 
#' sets of propensity scores and two propensity score distributions. Points 
#' rising to the right and a solid curve represent the weighting functions and 
#' distribution of propensity scores for estimating the average of the potential
#' outcomes without treatment whereas points rising to the left and a dashed 
#' curve represent the weighting functions and distribution of propensity scores
#' for estimating the average of the potential outcomes with treatment.
#'
#' Position of the legend is determined internally.
#'
#' @export
#'
#' @param object an object of class “nawt”, usually, a result of a call to \code{\link{nawt}}.
#'   Note that it cannot be used when the \code{object} is a result of a call to
#'   \code{nawt} where \code{method = "both"} and \code{twostep = FALSE}.
#' @param relative a logical value indicating whether or not relative weights 
#'   standardized to have mean one are shown.
#'
#' @return No retrun value, called for side effects.
#'
#' @author Hiroto Katsumata
#' 
#' @examples
#' # Simulation from Kang and Shafer (2007) and Imai and Ratkovic (2014)
#' tau <- 10
#' set.seed(12345)
#' n <- 1000
#' X <- matrix(rnorm(n * 4, mean = 0, sd = 1), nrow = n, ncol = 4)
#' prop <- 1 / (1 + exp(X[, 1] - 0.5 * X[, 2] + 0.25 * X[, 3] + 0.1 * X[, 4]))
#' treat <- rbinom(n, 1, prop)
#' y <- 210 + 27.4 * X[, 1] + 13.7 * X[, 2] + 13.7 * X[, 3] + 13.7 * X[, 4] + 
#'      tau * treat + rnorm(n)
#'
#' # Data frame and formulas for propensity score estimation
#' df <- data.frame(X, treat, y)
#' colnames(df) <- c("x1", "x2", "x3", "x4", "treat", "y")
#' formula_c <- as.formula(treat ~ x1 + x2 + x3 + x4)
#'
#' # Power weighting function with alpha = 2
#' # ATT estimation
#' fitatt <- nawt(formula = formula_c, outcome = "y", estimand = "ATT", 
#'                method = "score", data = df, alpha = 2)
#' plot_omega(fitatt)
#'
#' # ATE estimation
#' fitate <- nawt(formula = formula_c, outcome = "y", estimand = "ATE", 
#'                method = "score", data = df, alpha = 2)
#' plot_omega(fitate)
#'
#' # Use method = "both"
#' # Two-step estimation
#' fitateb2s <- nawt(formula = formula_c, outcome = "y", estimand = "ATE", 
#'                   method = "both", data = df, alpha = 2, twostep = TRUE)
#' plot_omega(fitateb2s)
#'
#' # Continuously-updating GMM estimation
#' \dontrun{
#' fitatebco <- nawt(formula = formula_c, outcome = "y", estimand = "ATE", 
#'                   method = "both", data = df, alpha = 2, twostep = FALSE)
#' plot_omega(fitatebco) # error}
plot_omega <- function (object, relative = TRUE) {
  if (length(object$omega) == 0) {
    stop("Use twostep = TRUE to use this function")
  }
  omega <- object$omega
  limx <- c(0, 1)
  if (relative == TRUE) {
    ylab <- expression(paste("Relative weights for the score ", omega(hat(pi))))
  } else {
    ylab <- expression(paste("Absolute weights for the score ", omega(hat(pi))))
  }
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  if (object$estimand != "ATE") {
    dens <- stats::density(object$ps, from = 0, to = 1)
    graphics::plot(dens,
                   xlim = limx,
                   ylim = c(0, max(dens$y)),
                   lwd = 1.2,
                   ann = FALSE,
                   axes = FALSE)
  } else { # ATE
    dens1 <- stats::density(object$ps[1, ], from = 0, to = 1)
    dens2 <- stats::density(object$ps[2, ], from = 0, to = 1)
    graphics::plot(dens1,
                   lty = 1,
                   col = grDevices::rgb(0.5, 0.5, 0.5),
                   xlim = limx,
                   ylim = c(0, max(c(dens1$y, dens2$y))),
                   lwd = 1.2,
                   ann = FALSE,
                   axes = FALSE)
    graphics::par(new = TRUE)
    graphics::plot(dens2,
                   lty = 2,
                   col = grDevices::rgb(0.2, 0.2, 0.2),
                   xlim = limx,
                   ylim = c(0, max(c(dens1$y, dens2$y))),
                   lwd = 1.2,
                   ann = FALSE,
                   axes = FALSE)
  }
  graphics::axis(side = 4)
  graphics::mtext("Density", side = 4, line = 3)
  graphics::par(new = TRUE)
  if (object$estimand == "ATE") {
    ps1 <- object$ps[1, ]
    ps2 <- object$ps[2, ]
    omega1 <- object$omega[1, ]
    omega2 <- object$omega[2, ]
    if (relative == TRUE) {
      omega1 <- omega1 / mean(omega1)
      omega2 <- omega2 / mean(omega2)
    }
    psc <- c(ps1[object$treat == 0], ps2[object$treat == 0])
    pst <- c(ps1[object$treat == 1], ps2[object$treat == 1])
    omegac <- c(omega1[object$treat == 0], omega2[object$treat == 0])
    omegat <- c(omega1[object$treat == 1], omega2[object$treat == 1])
    graphics::plot(psc, omegac,
                   col = grDevices::rgb(39/ 255, 139/ 255, 210 / 255, 
                                        alpha = 0.3),
                   xlim = limx,
                   ylim = c(0, max(c(omega1, omega2))),
                   xlab = "", 
                   ylab = "",
                   axes = FALSE)
    graphics::par(new = TRUE)
    graphics::plot(pst, omegat,
                   col = grDevices::rgb(220 / 255, 50 / 255, 46 / 255, 
                                        alpha = 0.3),
                   xlim = limx,
                   ylim = c(0, max(c(omega1, omega2))),
                   xlab = expression(paste("Estimated propensity score ", 
                                           hat(pi))), 
                   ylab = ylab,
                   axes = FALSE)
  } else { # ATT, MO,  and ATEcombined
    if (relative == TRUE) {
      omega <- omega / mean(omega)
    }
    omegac <- omega[object$treat == 0]
    omegat <- omega[object$treat == 1]
    graphics::plot(object$ps[object$treat == 0], omegac,
                   col = grDevices::rgb(39/ 255, 139/ 255, 210 / 255, 
                                        alpha = 0.3),
                   xlim = limx,
                   ylim = c(0, max(c(omegac, omegat))),
                   xlab = "", 
                   ylab = "",
                   axes = FALSE)
              graphics::par(new = TRUE)
    graphics::plot(object$ps[object$treat == 1], omegat,
                   col = grDevices::rgb(220 / 255, 50 / 255, 46 / 255, 
                                        alpha = 0.3),
                   xlim = limx,
                   ylim = c(0, max(c(omegac, omegat))),
                   xlab = expression(paste("Estimated propensity score ", 
                                           hat(pi))), 
                   ylab = ylab,
                   axes = FALSE)
  }
  if (length(unique(object$omega)) == 1) {
    legendy <- 0.2
  } else {
    legendy <- max(omegat, omegac)
  }
  graphics::axis(side = 1)
  graphics::axis(side = 2)
  graphics::box(bty = "u")
  graphics::legend(x = 0.1, y = legendy,
                   legend = c("treat", "control"), 
                   col = c(grDevices::rgb(220 / 255, 50 / 255, 46 / 255, 
                                          alpha = 0.75), 
                            grDevices::rgb(39/ 255, 139/ 255, 210 / 255, 
                                           alpha = 0.75)), 
                   pch = 21, yjust = 0.55, pt.lwd = 1.5, cex = 1.1,
                   x.intersp = 0.7, y.intersp = 0.9,
                   bty = "n",
                   bg = "transparent")
}
