## Transform results based on sigular value decomposition: coefficients
svdtranscoef <- function (coef, svdx) {
  c(MASS::ginv(diag(svdx$d) %*% t(svdx$v)) %*% coef)
}

## Transform results based on sigular value decomposition: variance
svdtransvcov <- function (varcov, svdx) {
  Kp1 <- nrow(varcov)
  res1 <- MASS::ginv(diag(svdx$d) %*% t(svdx$v)) %*% varcov[-Kp1, -Kp1] %*% 
            t(MASS::ginv(diag(svdx$d) %*% t(svdx$v)))
  res2 <- c(MASS::ginv(diag(svdx$d) %*% t(svdx$v)) %*% varcov[Kp1, -Kp1])
  rbind(as.matrix(cbind(res1, res2)), c(res2, varcov[Kp1, Kp1]))
}
