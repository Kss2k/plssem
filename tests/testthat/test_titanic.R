devtools::load_all()

head(titanic)

m <- "Survived ~ Age + Sex + Age:Sex"
fit <- pls(m, data = titanic, ordered = c("Survived", "Sex"),
           bootstrap = TRUE, boot.parallel = "multicore", 
           boot.ncpus = 4, boot.R = 50)

pls_predict(fit, benchmark = c("acc"))

fit <- pls(m, data = titanic, ordered = c("Survived"),
           bootstrap = TRUE, boot.parallel = "multicore", 
           boot.ncpus = 4, boot.R = 50, missing = "kNN", knn.k = 10)

testthat::expect_error(pls("Survived ~ Embarked", data = titanic),
                       regexp = "Please recode .*")
