#' Summarizing navigated weighting estimation
#'
#' Prints a summary of a fitted \code{nawt} object.
#'
#' Prints a summary of a \code{nawt} object, in a format similar to glm.
#'
#' @export
#'
#' @param object an object of class “nawt”, usually, a result of a call to \code{\link{nawt}}.
#' @param ... additional arguments to be passed to summary.
#'
#' @return 
#' \item{call}{the matched call.}
#' \item{est}{the point estimate of the parameter of interest.}
#' \item{coefficients}{a table including coefficients, standard errors, 
#'   z-values, and two-sided p-values.}
#' \item{effN_ps}{the effective sample size for the propensity score estimation.}
#' \item{effN_est}{the effective sample size for the parameter of interest
#'   estimation.}
#'
#' @author Hiroto Katsumata
#' 
#' @seealso \code{\link{nawt}}, \code{\link[base]{summary}}
#'
#' @examples # For examples see example(nawt)
summary.nawt <- function (object, ...) {
  est <- signif(object$est, digits = getOption("digits") - 3)
  k <- length(object$coef)
  se <- sqrt(diag(object$varcov))
  if (object$estimand == "ATE") {
    coef <- c(object$coef[1, ], object$coef[2, ])
  } else { # ATT, ATEcombined, and MO
    coef <- object$coef
  }
  coef.table <- cbind(as.vector(c(object$est, coef)),
                      as.vector(se[c(k + 1, 1:k)]),
                      as.vector(c(object$est, coef) / se[c(k + 1, 1:k)]),
                      as.vector(2 - 2 * stats::pnorm(abs(c(object$est, coef) / 
                                                         se[c(k + 1, 1:k)]))))
  colnames(coef.table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  if (object$estimand == "ATE") {
    rownames(coef.table) <- c("est", object$names.x, object$names.x)
  } else { # ATT, ATEcombined, and MO
    rownames(coef.table) <- c("est", object$names.x)
  }
  pval <- coef.table[, 4]
  symp <- stats::symnum(pval, corr = FALSE,
                        cutpoints = c(0, .001, .01, .05, .1, 1),
                        symbols = c("***", "**", "*", ".", " "))
  coef.print <- data.frame(coef.table, as.vector(symp))
  colnames(coef.print)[5] <- ""
  if (object$estimand == "ATE") {
    rownames(coef.print) <- 
      c("est", paste0(rep(c("ps1_", "ps2_"), each = k / 2), object$names.x))
  } else { # ATT, ATEcombined, and MO
    rownames(coef.print) <- c("est", object$names.x)
  }
  cat("\nCall: \n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
      "\n", sep = "")
  cat("\n", paste(object$estimand, "estimates: ", est), "\n", sep = "")
  cat("\nCoefficients:\n")
  print(coef.print, digits = getOption("digits") - 3)
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
  cat("\nEffective N for the propensity score estimation:", 
      round(object$effN_ps, digits = 2))
  if (object$estimand == "MO") {
    cat("\nEffective N for the", object$estimand, "estimation: ", 
        round(object$effN_est, digits = 2), "\n")
  } else {
    cat("\nEffective N for the", object$estimand, "estimation:\n")
    cat("  treatment:", round(object$effN_est[1], digits = 2), 
        "\n  control:  ", round(object$effN_est[2], digits = 2), "\n")
  }
  out <- list("call" = object$call, 
              "est" = object$est,
              "coefficients" = coef.table,
              "effN_ps" = object$effN_ps,
              "effN_est" = object$effN_est)
  invisible(out)
}
