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
  result$type <- ((result$type == "binary") * 3 + 21)
  result <- result[nrow(result):1, ]
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
    result <- result[order.un, ]
  }
  mindiff <- min(c(0, min(c(result$diff.un, result$diff.adj))))
  maxdiff <- max(c(result$diff.un, result$diff.adj))
  if (absolute == TRUE) {
    legendx <- (maxdiff * 5 / 8)
  } else {
    legendx <- mindiff * 1.1
  }
  cols0 <- c(grDevices::rgb(0 / 255, 184 / 255, 148 / 255), 
             grDevices::rgb(225 / 255, 112 / 255, 85 / 255))
  cols <- rep(cols0, each = 2)
  pchs <- rep(c(21, 24), 2)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::plot(x = result$diff.un, 
                 y = result$covariates,
                 pch = result$type,
                 col = cols0[1],
                 cex = 1.7,
                 lwd = 2.2,
                 xlim = c(mindiff, maxdiff),
                 xlab = "", ylab = "",
                 axes = FALSE)
  graphics::par(new = TRUE)
  graphics::plot(x = result$diff.adj, 
                 y = result$covariates, 
                 pch = result$type,
                 col = cols0[2],
                 cex = 1.7,
                 lwd = 2.2,
                 xlim = c(mindiff, maxdiff),
                 xlab = diff, ylab = "",
                 yaxt = "n",
                 main = "Covariate balance")
  graphics::abline(v = 0, col = "grey10", lty = "solid")
  graphics::abline(v = threshold, col = "grey50", lty = "dashed", lwd = 1.2)
  graphics::axis(2, at = c(1:nrow(result)), labels = result$covariates, las = 1)
  graphics::par(xpd = TRUE)
  if (sort == FALSE) {
    legendx <- graphics::par()$usr[2]
  }
  graphics::legend(x = legendx, y = 1,
                   legend = c("Adjusted: continuous", "Adjusted: binary",
                               "Unadjusted: continuous", "Unadjusted: binary"), 
                   col = cols[4:1], pch = pchs, pt.cex = 1.5, pt.lwd = 2, yjust = 0,
                   x.intersp = 0.7, y.intersp = 0.9,
                   bty = "n",
                   bg = "transparent")
}
