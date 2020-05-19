## Calculate weights for the score omega(hat(pi))
weight_omega <- function (result) {
  if (result$estimand == "ATEcombined") {
    if (result$method == "score") {
      omega <- result$ps^result$alpha + (1 - result$ps)^result$alpha
    } else if (result$method == "cb") {
      omega <- 1 / (1 - result$ps) + 1 / result$ps
    } else if (result$method == "both") {
      omega <- (result$ps^result$alpha + (1 - result$ps)^result$alpha) * 
                 result$scratio +
                (1 / (1 - result$ps) + 1 / result$ps) * (1 - result$scratio)
    }
  } else if (result$estimand == "ATE") {
    ps1 <- result$ps[1, ]
    ps2 <- result$ps[2, ]
    if (result$method == "score") {
      omega1 <- ps1^result$alpha
      omega2 <- (1 - ps2)^result$alpha
    } else if (result$method == "cb") {
      omega1 <- 1 / (1 - ps1)
      omega2 <- 1 / ps2
    } else if (result$method == "both") {
      omega1 <- ps1^result$alpha * result$scratio[1] +
                 1 / (1 - ps1) * (1 - result$scratio[1])
      omega2 <- (1 - ps2)^result$alpha * result$scratio[2] +
                 1 / ps2 * (1 - result$scratio[2])
    }
    omega <- rbind(omega1 = omega1, omega2 = omega2)
  } else { # ATT and MO
    if (result$method == "score") {
      omega <- result$ps^result$alpha
    } else if (result$method == "cb") {
      omega <- 1 / (1 - result$ps)
    } else if (result$method == "both") {
      omega <- result$ps^result$alpha * result$scratio +
                1 / (1 - result$ps) * (1 - result$scratio)
    }
  }
  omega
}
