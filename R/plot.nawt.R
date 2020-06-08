#' Plot a scattered plot comparing the navigated weighting and naive estimation
#'
#' Plots a scattered plot comparing the resulting inverse probability weights 
#' estimated by the navigated weighting and the standard logistic regression.
#'
#' The x-axis shows the inverse probability weights estimated by estimating 
#' propensity scores with the standard logistic regression whereas the y-axis
#' shows those with the navigated weighting. Excessively heavy weights on only a
#' few observations in the navigated weighting may indicate the failure of 
#' the estimation.
#'
#' Position of the legend is determined internally.
#'
#' @export
#'
#' @param x an object of class “nawt”, usually, a result of a call to \code{\link{nawt}}.
#' @param ... additional arguments to be passed to plot.
#'
#' @return No retrun value, called for side effects.
#'
#' @author Hiroto Katsumata
#' 
#' @seealso \code{\link{nawt}}, \code{\link[graphics]{plot}}
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
#' plot(fitatt)
#'
#' # ATE estimation
#' fitate <- nawt(formula = formula_c, outcome = "y", estimand = "ATE", 
#'                method = "score", data = df, alpha = 2)
#' plot(fitate)
plot.nawt <- function (x, ...) {
  maxweight <- max(c(x$naive_weights, x$weights))
  minweight <- min(c(x$naive_weights, x$weights))
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(pty = "s")
  graphics::plot(x = x$naive_weights[x$treat == 0], 
                 y = x$weights[x$treat == 0], 
                 col = grDevices::rgb(39/ 255, 139/ 255, 210 / 255, 
                                      alpha = 0.75),
                 xlim = c(0, maxweight), ylim = c(0, maxweight),
                 xlab = "", ylab = "")
  graphics::par(new = TRUE)
  graphics::plot(x = x$naive_weights[x$treat == 1], 
                 y = x$weights[x$treat == 1], 
                 col = grDevices::rgb(220 / 255, 50 / 255, 46 / 255, 
                                      alpha = 0.75),
                 xlim = c(0, maxweight), ylim = c(0, maxweight),
                 xlab = "Naive weights", ylab = "Navigated weights")
  graphics::text(x = minweight, y = maxweight * 0.9, labels = "control", 
                 col = grDevices::rgb(39/ 255, 139/ 255, 210 / 255, 
                                      alpha = 0.85), 
                 adj = 0, cex = 1.2)
  graphics::text(x = minweight, y = maxweight * 0.95, labels = "treat", 
                 col = grDevices::rgb(220 / 255, 50 / 255, 46 / 255, 
                                      alpha = 0.85), 
                 adj = 0, cex = 1.2)
  graphics::abline(0, 1)
}
