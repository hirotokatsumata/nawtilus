## inverse probability weights
invprob <- function (ps, missing, estimand, weights) {
  invp <- numeric(length(missing))
  if (estimand == "ATT") {
    invp[missing == 0] <- ps[missing == 0] / (1 - ps[missing == 0]) * weights[missing == 0]
    invp[missing == 1] <- weights[missing == 1] / sum(weights[missing == 1])
    invp[missing == 0] <- invp[missing == 0] / sum(invp[missing == 0])
  } else if (estimand == "MO") {
    invp[missing == 0] <- 1 / (1 - ps[missing == 0]) * weights[missing == 0]
    invp[missing == 1] <- 0
    invp <- invp / sum(invp)
  } else { # ATE
    invp[missing == 0] <- 1 / (1 - ps[missing == 0]) * weights[missing == 0]
    invp[missing == 1] <- 1 / ps[missing == 1] * weights[missing == 1]
    invp[missing == 1] <- invp[missing == 1] / sum(invp[missing == 1])
    invp[missing == 0] <- invp[missing == 0] / sum(invp[missing == 0])
  }
  invp
}
