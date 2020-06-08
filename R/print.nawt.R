#' Print navigated weighting estimation
#'
#' Prints a fitted \code{nawt} object.
#'
#' @export
#'
#' @param x an object of class “nawt”, usually, a result of a call to \code{\link{nawt}}.
#' @param ... additional arguments to be passed to print.
#'
#' @return No retrun value, called for side effects.
#'
#' @author Hiroto Katsumata
#' 
#' @seealso \code{\link{nawt}}, \code{\link[base]{print}}
#'
print.nawt <- function (x, ...) {
  est <- signif(x$est, digits = getOption("digits") - 3)
  if (x$estimand == "ATE") {
    coef <- c(x$coef[1, ], x$coef[2, ])
    names(coef) <- c(paste0(rep(c("ps1_", "ps2_"), each = length(coef) / 2), 
                            x$names.x))
  } else { # ATT, ATEcombined, and MO
    coef <- x$coef
  }
  cat("\nCall: \n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n", sep = "")
  cat("\n", paste(x$estimand, "estimates: ", est), "\n", sep = "")
  cat("\nCoefficients:\n")
  print(coef, digits = getOption("digits") - 3)
  cat("\nEffective N for propensity score estimation:", 
      round(x$effN_ps, digits = 2))
  if (x$estimand == "MO") {
    cat("\nEffective N for the", x$estimand, "estimation: ", 
        round(x$effN_est, digits = 2), "\n")
  } else {
    cat("\nEffective N for the", x$estimand, "estimation:\n")
    cat("  treatment:", round(x$effN_est[1], digits = 2), 
        "\n  control:  ", round(x$effN_est[2], digits = 2), "\n")
  }
}
