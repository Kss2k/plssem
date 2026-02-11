devtools::load_all()

rthreshold <- \(k, offset = runif(1, min = -0.7, max = 0.7), sigma = 0.2) {
  t <- seq_len(k) - mean(seq_len(k)) + offset
  t <- t + runif(k, min = -sigma, max = sigma)
  c(-Inf, t, Inf)
}

set.seed(2308257)
K <- 6

cut_data <- function(data, k = K, choose = NULL) {
  standardize <- \(x) (x - mean(x)) / sd(x)

  if (is.null(choose))
    choose <- setdiff(colnames(data), "cluster")

  means <- data

  for (i in seq_along(choose)) {
    var <- choose[[i]]
    x <- data[[var]]
    breaks <- rthreshold(k = k)
    y <- cut(standardize(x), breaks = breaks)
    y <- as.ordered(as.integer(y))

    data[[var]] <- y
    z <- rep(NA_real_, length(y))

    for (v in unique(y))
        z[y == v] <- mean(x[y == v])

    means[[var]] <- z
  }

  # list(data = data, means = as.data.frame(means))
  data
}


randomInterceptsOrdered <- cut_data(randomIntercepts)
randomSlopesOrdered <- cut_data(randomSlopes)
TPB_Ordered <- cut_data(modsem::TPB)
oneIntOrdered <- cut_data(modsem::oneInt)


save(randomSlopesOrdered, file = "data/randomSlopesOrdered.rda")
save(randomInterceptsOrdered, file = "data/randomInterceptsOrdered.rda")
save(TPB_Ordered, file = "data/TPB_Ordered.rda")
save(oneIntOrdered, file = "data/oneIntOrdered.rda")
