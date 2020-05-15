context("Variance-covariance matrix")

test_that("standard errors are correct", {
  skip_on_cran() # Takes long time
  skip("skip: manually comment out to run the test") # Takes long time.
  # Simulation
  # True treatment effect is 10
  tau <- 10
  set.seed(12345)
  n <- 10000
  X <- matrix(runif(n * 4, -1.5, 1.5), nrow = n, ncol = 4)
  prop <- 1 / (1 + exp(X[, 1] - 0.5 * X[, 2] + 0.25 * X[, 3] + 0.1 * X[, 4]))
  treat <- rbinom(n, 1, prop)
  y <- 210 + 27.4 * X[, 1] + 13.7 * X[, 2] + 13.7 * X[, 3] + 13.7 * X[, 4] +
    tau * treat + rnorm(n)

  # Data frame and formulas for propensity score estimation
  df <- data.frame(X, treat, y)
  colnames(df) <- c("x1", "x2", "x3", "x4", "treat", "y")
  formula_c <- as.formula(treat ~ x1 + x2 + x3 + x4)

  # tolerance
  tol <- 2e-2
  ones1 <- rep(1, 6)
  names(ones1) <- c("(Intercept)", "x1", "x2", "x3", "x4", "est")
  ones2 <- rep(1, 11)
  names(ones2) <- c(paste0(rep(c("ps1_", "ps2_"), each = 5), 
                           c("(Intercept)", "x1", "x2", "x3", "x4")), 
                    "est")

  # Bootstrap replicates
  B <- 2000

  # ATT estimation
  # Power weighting function with alpha = 2
  fits2a <- nawt(formula = formula_c, outcome = "y", estimand = "ATT",
                 method = "score", data = df, alpha = 2)
  fits2b <- nawt(formula = formula_c, outcome = "y", estimand = "ATT",
                 method = "score", data = df, alpha = 2, boot = TRUE, B = B)
  expect_equal(sqrt(diag(fits2a$varcov)) / sqrt(diag(fits2b$varcov)), ones1,
               tolerance = tol)

  # Covariate balancing weighting function
  fitcba <- nawt(formula = formula_c, outcome = "y", estimand = "ATT",
                 method = "cb", data = df)
  fitcbb <- nawt(formula = formula_c, outcome = "y", estimand = "ATT",
                 method = "cb", data = df, boot = TRUE, B = B)
  expect_equal(sqrt(diag(fitcba$varcov)) / sqrt(diag(fitcbb$varcov)), ones1,
               tolerance = tol)

  # Standard logistic regression
  fits0a <- nawt(formula = formula_c, outcome = "y", estimand = "ATT",
                 method = "score", data = df, alpha = 0)
  fitlog <- vcov(glm(formula = formula_c, family = binomial, data = df))
  expect_equal(sqrt(diag(fits0a$varcov)[-6]) / sqrt(diag(fitlog)), ones1[-6],
               tolerance = tol)

  # ATE estimation
  # Power weighting function with alpha = 2
  fits2a2 <- nawt(formula = formula_c, outcome = "y", estimand = "ATE",
                  method = "score", data = df, alpha = 2)
  fits2b2 <- nawt(formula = formula_c, outcome = "y", estimand = "ATE",
                  method = "score", data = df, alpha = 2, boot = TRUE, B = B)
  expect_equal(sqrt(diag(fits2a2$varcov)) / sqrt(diag(fits2b2$varcov)), ones2,
               tolerance = tol)

  # Covariate balancing weighting function
  fitcba2 <- nawt(formula = formula_c, outcome = "y", estimand = "ATE",
                  method = "cb", data = df)
  fitcbb2 <- nawt(formula = formula_c, outcome = "y", estimand = "ATE",
                  method = "cb", data = df, boot = TRUE, B = B)
  expect_equal(sqrt(diag(fitcba2$varcov)) / sqrt(diag(fitcbb2$varcov)), ones2,
               tolerance = tol)
})
