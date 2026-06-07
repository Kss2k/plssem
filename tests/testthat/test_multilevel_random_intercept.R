devtools::load_all()

syntax <- '
  f =~ y1 + y2 + y3
  f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)
'

testthat::expect_no_error({
  fit <- pls(syntax, data = randomIntercepts, bootstrap = FALSE)
  summary(fit)
})


syntax <- '
  f =~ y1 + y2 + y3
  f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)
'

testthat::expect_no_error({
  fit <- pls(syntax, data = randomInterceptsOrdered, bootstrap = TRUE)
  summary(fit)
})


df <- randomIntercepts

df$w_y1 <- pls_within_values(df$y1, df$cluster)
df$w_y2 <- pls_within_values(df$y2, df$cluster)
df$w_y3 <- pls_within_values(df$y3, df$cluster)
df$b_y1 <- pls_between_values(df$y1, df$cluster)
df$b_y2 <- pls_between_values(df$y2, df$cluster)
df$b_y3 <- pls_between_values(df$y3, df$cluster)

syntax <- '
  fw =~ w_y1 + w_y2 + w_y3
  fb =~ b_y1 + b_y2 + b_y3
  fw ~ x1 + x2 + x3 + (1 | cluster)
  fb ~ w1 + w2 + (1 | cluster)
'

fit <- pls(syntax, data = df)
summary(fit)

syntax.lav <- '
    level: 1
        fw =~ y1 + y2 + y3
        fw ~ x1 + x2 + x3
    level: 2
        fb =~ y1 + y2 + y3
        fb ~ w1 + w2
'

library(lavaan)
summary(sem(syntax.lav, df, cluster = "cluster"))
