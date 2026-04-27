devtools::load_all()


# Third-order higher-order model with a structural equation mixing
# first-, second-, and third-order predictors.

set.seed(812341)
N <- 1200L

std <- function(x) (x - mean(x)) / stats::sd(x)

indicator <- function(x, lambda) {
  x <- std(x)
  sigma <- sqrt(1 - lambda^2)
  lambda * x + stats::rnorm(N, mean = 0, sd = sigma)
}

# Third-order latent variables
ABCD <- std(stats::rnorm(N))

# Second-order latent variables
AB <- indicator(ABCD, 0.9)
CD <- indicator(ABCD, 0.8)
EF <- std(stats::rnorm(N))

# First-order latent variables
A <- indicator(AB, 0.9)
B <- indicator(AB, 0.8)
C <- indicator(CD, 0.9)
D <- indicator(CD, 0.8)
E <- indicator(EF, 0.9)
F <- indicator(EF, 0.8)

Yproj <- 0.5 * ABCD + 0.4 * EF + 0.3 * A
Y <- Yproj + stats::rnorm(N, mean = 0, sd = sqrt(1 - stats::var(Yproj)))

# Observed indicators
a1 <- indicator(A, 0.9)
a2 <- indicator(A, 0.8)
a3 <- indicator(A, 0.7)

b1 <- indicator(B, 0.9)
b2 <- indicator(B, 0.8)
b3 <- indicator(B, 0.7)

c1 <- indicator(C, 0.9)
c2 <- indicator(C, 0.8)
c3 <- indicator(C, 0.7)

d1 <- indicator(D, 0.9)
d2 <- indicator(D, 0.8)
d3 <- indicator(D, 0.7)

e1 <- indicator(E, 0.9)
e2 <- indicator(E, 0.8)
e3 <- indicator(E, 0.7)

f1 <- indicator(F, 0.9)
f2 <- indicator(F, 0.8)
f3 <- indicator(F, 0.7)

y1 <- indicator(Y, 0.9)
y2 <- indicator(Y, 0.8)
y3 <- indicator(Y, 0.7)

data <- data.frame(
  a1, a2, a3,
  b1, b2, b3,
  c1, c2, c3,
  d1, d2, d3,
  e1, e2, e3,
  f1, f2, f3,
  y1, y2, y3
)


syntax <- '
  A =~ a1 + a2 + a3
  B =~ b1 + b2 + b3
  C =~ c1 + c2 + c3
  D =~ d1 + d2 + d3
  E =~ e1 + e2 + e3
  F =~ f1 + f2 + f3

  AB =~ A + B
  CD =~ C + D
  EF =~ E + F
  ABCD =~ AB + CD

  Y =~ y1 + y2 + y3

  # Mix 3rd-, 2nd-, and 1st-order predictors
  Y ~ ABCD + EF + A
'


testthat::expect_no_error({
  fit <- pls(syntax, data = data, consistent = TRUE)
  summary(fit)

  pt <- parameter_estimates(fit)
  testthat::expect_true(any(pt$lhs == "Y" & pt$op == "~" & pt$rhs == "ABCD"))
  testthat::expect_true(any(pt$lhs == "Y" & pt$op == "~" & pt$rhs == "EF"))
  testthat::expect_true(any(pt$lhs == "Y" & pt$op == "~" & pt$rhs == "A"))
})
