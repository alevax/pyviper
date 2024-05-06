## Import packages
library(tidyverse)
library(patchwork)
library(PISCES)
library(viper)
library(patchwork)

# Create an option parser
option_list <- list(
  make_option(
    c("--input", "-i"), # Option name and short form
    type = "character",   # Option type 
    help = "input matrix"  # Option description
  ),
  make_option(
    c("--network1", "-n1"),
    type = "character",   # Option type 
    help = "network1 to use"  # Option description
  ),
  make_option(
    c("--network2", "-n2"),
    type = "character",   # Option type 
    help = "network2 to use"  # Option description
  ),
  make_option(
    c("--output", "-o"),
    type = "character",
    default = '',# Option type (integer in this case)
    help = "output name"  # Option description
  )
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

## Help function

duplicate_Removal_conversion<-function(ARACNe3_regulon){
  
  elements_to_remove<-names(ARACNe3_regulon)[which(duplicated(names(ARACNe3_regulon))==TRUE)]
  
  ARACNe3_regulon_NoDuplicates<-ARACNe3_regulon
  
  ARACNe3_regulon_NoDuplicates[elements_to_remove]<-NULL
  
  ARACNe3_regulon_NoDuplicates_converted<-list()
  
  for (i in 1:length(ARACNe3_regulon_NoDuplicates)){
    
    tfmode<-ARACNe3_regulon_NoDuplicates[[i]]$am
    
    likelihood<-unname(ARACNe3_regulon_NoDuplicates[[i]]$aw)
    
    tmp_list<-list(tfmode, likelihood)
    
    names(tmp_list)<-c("tfmode", "likelihood")
    
    ARACNe3_regulon_NoDuplicates_converted[[i]]<-tmp_list
    
    names(ARACNe3_regulon_NoDuplicates_converted)[i]<-names(ARACNe3_regulon_NoDuplicates)[i]
    
  }
  
  class(ARACNe3_regulon_NoDuplicates_converted)<-"regulon"
  
  return(ARACNe3_regulon_NoDuplicates_converted)
  
}

## 
ges = readRDS(opt$ges)
test_net1 = readRDS(opt$network1) 
test_net2 = readRDS(opt$network2) 

net.list = list(test_net1,test_net2)

op_narnea = meta_narnea(ges, net.list, sample.weights = TRUE)


net.list = lapply(net.list, duplicate_Removal_conversion)
op_viper = viper(ges, net.list, method = 'none', minsize = 30)

### load pydata

readdata = function(thepath){
  df = read_csv(thepath)
  indexes =  df %>% pull(...1)
  df = df %>%
    select(-...1)
  rownames(df) = indexes
  df = t(df)
  
  return(df)
}



op_nes_py = readdata('testpyviper_narnea_nes.csv') 
op_pes_py = readdata('testpyviper_narnea_pes.csv') 
viper_py = readdata("testpyviper_area_nes.csv")


### calculate diff

caldiff = function(r_score, py_score){
  
  r_score = r_score[,colnames(py_score)]
  r_score = r_score[rownames(py_score),]
  
  cor = cor(as.vector(r_score), as.vector(py_score))
  cor = sprintf("%.2f", round(cor, 3))
  ass_dif = mean(abs(as.vector(r_score)-as.vector(py_score))/abs(as.vector(r_score)),na.rm = TRUE)
  max_dif = max(abs(as.vector(r_score)-as.vector(py_score)))

  return(c(cor, ass_dif, max_dif))
}

narnea_nes = caldiff(op_narnea$NES, op_nes_py)
narnea_pes = caldiff(op_narnea$PES, op_pes_py)
viper_nes = caldiff(op_viper,viper_py)

op_table = rbind(narnea_nes, narnea_pes, viper_nes)
op_table = data.frame(op_table)
colnames(op_table) = c("cor", "ass_dif","max_dif")

print(op_table)

write.table(op_table, file = "test_op.txt", sep = ",", col.names = TRUE, row.names = TRUE)