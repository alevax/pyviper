# ------------------------------------------------------------
## Libraries
## ------------------------------------------------------------
library(viper)
library(Biobase)

## ------------------------------------------------------------
# Set up directory
## ------------------------------------------------------------
## setwd("/Users/friva/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/Desktop/SBP_CODE/PleiotropyProject/BRCA")
setwd("/Users/rcassius/projects/pyviper/tests/resources")

## ------------------------------------------------------------
## Helper: table -> regulon (same shape viper expects)
## ------------------------------------------------------------
TableToInteractome <- function(net_table){
  blank_reg <- list(c(), c()); names(blank_reg) <- c("tfmode", "likelihood")
  my_reg <- list()
  pb <- txtProgressBar(min = 0, max = nrow(net_table), initial = 0, style = 3)
  for(i in 1:nrow(net_table)){
    my_reg[[net_table[i,"regulator"]]] <- blank_reg
    setTxtProgressBar(pb, i)
  }
  close(pb)
  pb <- txtProgressBar(min = 0, max = nrow(net_table), initial = 0, style = 3)
  for(i in 1:nrow(net_table)){
    new_w <- net_table[i,"likelihood"]; new_m <- net_table[i,"mor"]
    names(new_w) <- net_table[i,"target"]; names(new_m) <- net_table[i,"target"]
    my_reg[[net_table[i,"regulator"]]][["likelihood"]] <- c(my_reg[[net_table[i,"regulator"]]][["likelihood"]], new_w)
    my_reg[[net_table[i,"regulator"]]][["tfmode"]]     <- c(my_reg[[net_table[i,"regulator"]]][["tfmode"]],     new_m)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  my_reg
}

## ------------------------------------------------------------
## 1) Load data
## ------------------------------------------------------------
# expr <- read.delim("DEG-BRCA.tsv", row.names = 1, check.names=FALSE)
expr <- read.csv("ges.csv", row.names = 1, check.names=FALSE)
# expr_matrix <- t(as.matrix(expr))  # samples in columns
expr_matrix <- as.matrix(expr)  # samples in columns
dset <- ExpressionSet(assayData = expr_matrix)

# interactome <- read.delim("Interactome_prunedBRCA1.tsv", stringsAsFactors = FALSE)
interactome <- read.delim("test_net1.tsv", stringsAsFactors = FALSE)
regulon <- TableToInteractome(interactome)


## ------------------------------------------------------------
## 2) VIPER pre-pleiotropy NES
## ------------------------------------------------------------
tf_nes_original <- viper(
  exprs(dset), regulon,
  pleiotropy = FALSE,
  method = "none",
  nes = TRUE,
  eset.filter = TRUE,
  minsize = 0,    #disables regulon-size filtering
  cores = 1,
  verbose = TRUE
)
# write.csv(tf_nes_original, "nes_RVIPER_precorrectionBRCA1.csv", row.names = TRUE)
write.csv(tf_nes_original, "viper_nes_R_output.csv", row.names = TRUE)

tf_nes_corrected <- viper(
  exprs(dset), regulon,
  pleiotropy = TRUE,
  method = "none",
  nes = TRUE,
  eset.filter = TRUE,
  minsize = 0,    #disables regulon-size filtering
  cores = 1,
  verbose = TRUE
)
# write.csv(tf_nes_original, "nes_RVIPER_precorrectionBRCA1.csv", row.names = TRUE)
write.csv(tf_nes_original, "viper_nes_R_output.csv", row.names = TRUE)



dd_default <- viperSimilarity( x = tf_nes_original, nn = NULL, ws = c(4,2), method = c("two.sided")) 
write.csv(dd_default, file = "viper_similarity_R_output_two_sided.csv", row.names = TRUE)


dd_default_greater <- viperSimilarity( x = tf_nes_original, nn = NULL, ws = c(4,2), method = c("greater"))   
write.csv(dd_default_greater, file = "viper_similarity_R_output_greater.csv", row.names = TRUE)


dd_default_less <- viperSimilarity( x = tf_nes_original, nn = NULL, ws = c(4,2), method = c("less"))   
write.csv(dd_default_less, file = "viper_similarity_R_output_less.csv", row.names = TRUE)

dd_default_50 <- viperSimilarity( x = tf_nes_original, nn = 50, ws = c(4,2), method = c("less"))              # returns a signatureDistance (square matrix)
write.csv(dd_default_50, file = "viper_similarity_R_output_50.csv", row.names = TRUE)



