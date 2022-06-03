## Libraries
# print("Loading optparse and VIPER libraries and functions...")
print("Running test_aREA_func.R...")
print("Loading optparse library and functions...")
suppressMessages(suppressWarnings(library(optparse)))
# suppressMessages(suppressWarnings(library(viper)))

## Functions
source("../../libs/area_fn.R")

## Options
option_list <- list(
  make_option(
    c('-q', '--quiet'),
    action = 'store_false',
    default = TRUE,
    help = 'Suppresses status updates.',
    dest = 'print'
  ),
  make_option(c('-g', '--ges'), type = "character", help = 'Gene expression signature'),
  make_option(c('-t', '--network'), type = "character", help = 'Network regulon'),
  make_option(c('-n', '--out_name_project'), type = "character", help = 'Output files name convention: project name.'),
  make_option(c('-o', '--out_dir'), type = "character", help = 'Directory for output.')
)
opt <- parse_args(OptionParser(option_list = option_list))

print("Retrieving input files...")
ges <- readRDS(opt$ges) #readRDS("test_files_1/LNCaPWT_gExpr_GES.rds")
network <- readRDS(opt$network) #readRDS("test_files_1/LNCaPWT_pruned.rds")
print("Running aREA...")
vip <- aREA(ges, network)
print("Saving results...")
outfile <- paste0(opt$out_dir, opt$out_name_project, "_aREA_PAct.rds")
saveRDS(vip$nes, outfile)
print("Done.")
