devtools::load_all()

# Simulate simple data
set.seed(23948)
n <- 1000
# Generate formative indicators v1 and v2 with correlation r
r <- 0.6
v1 <- rnorm(n)
v2 <- r * v1 + rnorm(n, sd = sqrt(1 - r^2))

# Form A from v1 and v2
A <- v1 + 0.8 * v2

# Form v3 and v4 from A
v3 <- A + rnorm(n, sd = sqrt(0.2))
v4 <- 0.8 * A + rnorm(n, sd = sqrt(0.2))

data <- data.frame(v1, v2, v3, v4)

# Minimal syntax with both reflective and formative indicators
model <- '
  # Reflective indicators
  A =~ v3 + v4

  # Formative/Composite indicators
  A <~ v1 + v2
'

testthat::expect_no_error({
  fit <- pls(model, data = data)
  unstandardized_estimates(fit)
  pls_predict(fit, newdata = data)
  pls_predict(fit)
})
