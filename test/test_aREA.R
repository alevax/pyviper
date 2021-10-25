##
# Test Functions for Enrichment Analysis in R
# -------------------------------------------
# source("../vaxtools/R/gsea-plot.R")

source("libs/aREA_functions.R")
library(viper)

data(bcellViper, package = 'bcellViper')

d1 <- Biobase::exprs(dset)
res <- viper::rowTtest(dset, "description", "CB", "N")

ges <- res$statistic[,1]

gsr <- gsea_regulon(ges, regulon$MYB)
plot(gsr, plotSignature = TRUE,
     signatureNames = c('CB', 'N'),
     signatureType = 't-statistic')