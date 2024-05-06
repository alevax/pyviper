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
    default = 'test',# Option type (integer in this case)
    help = "output name"  # Option description
  )
)

# Parse the command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


## 
#let's set csv to be both input
#test_exp = read.csv(opt$input,row.names = 1)
#use this to saveRDS(mtcars, "mtcars.rds")
ges = readRDS(opt$ges)
test_net1 = readRDS(opt$network1) 
test_net2 = readRDS(opt$network2) 
#prepare data
write.csv(ges, file = "ges.csv")

InteractomeToTable = function(net.obj, out.file, netfrom = c('ARACNe3','ARACNe_ap' )) {
  # make df
  if (netfrom == 'ARACNe3'){
    tfmode = 'am'
    likelihood = 'aw'
  }
  if (netfrom == 'ARACNe_ap'){
    tfmode = 'tfmode'
    likelihood = 'likelihood'
  }
  
  
  
  
  net.df <- do.call(rbind, lapply(names(net.obj), function(x) {
    reg.obj <- net.obj[[x]] 
    if (is.null(reg.obj$subnets))
    {
      reg.df <- data.frame('regulator' = rep(x, length(reg.obj[[tfmode]])),
                           'target' = names(reg.obj[[tfmode]]),
                           'mor' = reg.obj[[tfmode]],
                           'likelihood' = reg.obj[[likelihood]])
    } else {
      reg.df <- data.frame('regulator' = rep(x, length(reg.obj[[tfmode]])),
                           'target' = names(reg.obj[[tfmode]]),
                           'mor' = reg.obj[[tfmode]],
                           'likelihood' = reg.obj[[likelihood]],
                           'subnets' = reg.obj$subnets)
    }
    
    rownames(reg.df) <- NULL
    return(reg.df)
  }))
  # if specified, write to file
  if (!missing(out.file)) {
    write.table(net.df, file = out.file, quote = FALSE, sep = '\t',
                row.names = FALSE, col.names = TRUE)
  } else {
    return(net.df)
  }
}


test_net1_csv = InteractomeToTable(test_net1, netfrom = 'ARACNe3', 'test_net1.tsv')
test_net2_csv = InteractomeToTable(test_net2, netfrom = 'ARACNe3', 'test_net2.tsv')

