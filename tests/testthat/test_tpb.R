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
  BEH ~ INT + PBC + INT:PBC
'

fit <- pls(tpb, modsem::TPB, bootstrap = TRUE, sample = 50)
summary(fit)


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

fit <- pls(tpb, TPB_Ordered, bootstrap = TRUE)
summary(fit)


TPB_UK <- modsem::TPB_UK

tpb_uk <- "
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att3 + att2 + att1 + att4
  SN =~ sn4 + sn2 + sn3 + sn1
  PBC =~ pbc2 + pbc1 + pbc3 + pbc4
  INT =~ int2 + int1 + int3 + int4
  BEH =~ beh3 + beh2 + beh1 + beh4

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC
"

mcpls(tpb_uk, TPB_UK, ordered = colnames(TPB_UK),
      max.iter.mc = 500, consistent = FALSE)
