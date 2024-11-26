library(modsem)
library(cSEM)
devtools::load_all()

tpb <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC 
'

# cSEM -------------------------------------------------------------------------
start1 <- Sys.time()
est_cSEM <- csem(.data = modsem::TPB, .model = tpb, 
                 .approach_weights = "PLS-PM", 
                 .disattenuate = FALSE, 
                 .resample_method = "bootstrap", .R = 3000)
time1 <- Sys.time() - start1
summarize(est_cSEM)

start2 <- Sys.time()
est <- pls(tpb, modsem::TPB, consistent = TRUE, 
           standardize = FALSE, convergence = 1e-10, 
           bootstrap = FALSE, sample = 3000)
time2 <- Sys.time() - start2
est$fit
as.data.frame(est$params)
est$factorScores


# Interaction ------------------------------------------------------------------
tpbInt <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC 
  BEH ~ INT:PBC  
'
#estInt <- pls(tpbInt, modsem::TPB, consistent = TRUE, 
#              standardize = TRUE, convergence = 1e-10, 
#              bootstrap = TRUE)
#getReliabilityCoefs(estInt)
#estInt$fit
#estInt$params |> as.data.frame()


# Using product indicator approach ---------------------------------------------
modsem(tpbInt, data = modsem::TPB, standardize = TRUE)$newData -> tpbInt
tpbIndProds <- '
ATT =~ att1 + att2 + att3 + att4 + att5
SN =~ sn1 + sn2
PBC =~ pbc1 + pbc2 + pbc3
INT =~ int1 + int2 + int3
BEH =~ b1 + b2
INT ~ ATT + SN + PBC
BEH ~ INT + PBC + INTPBC
INTPBC =~ int1pbc1 + int2pbc1 + int3pbc1 + 
  int1pbc2 + int2pbc2 + int3pbc2 +
  int1pbc3 + int2pbc3 + int3pbc3
'

estIndProdsCsem <- csem(.data = tpbInt, .model = tpbIndProds, 
                        .approach_weights = "PLS-PM", 
                        .disattenuate = TRUE)
summarize(estIndProdsCsem)
estIndProds <- pls(tpbIndProds, tpbInt, standardize = TRUE, 
                   consistent = TRUE, convergence = 1e-5)
estIndProds$fit

covProdInds(c("pbc1", "pbc2", "pbc3"), c("int1", "int2", "int3"), 
            modsem::TPB)
