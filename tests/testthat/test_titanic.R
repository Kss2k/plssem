devtools::load_all()

head(titanic)

m <- "Survived ~ Age + Sex + Age:Sex"

fit <- pls(m, data = titanic, ordered = c("Survived", "Sex"),
           bootstrap = TRUE, boot.parallel = "multicore", 
           boot.ncpus = 4, boot.R = 500)

pls_predict(fit, benchmark = c("acc"))
