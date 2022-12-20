## Libraries
print("Running test_InteractomeToTable_func.R...")
print("Loading optparse and viper libraries...")
suppressMessages(suppressWarnings(library(optparse)))

# suppressMessages(suppressWarnings(library(viper)))
source("../../libs/r-helper.R")

## Options
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-n', '--net_obj'), type="character", help='Transcription factor and cotranscription factor network regulon'),
  make_option(c('-o', '--out_name_project'), type="character", help='Naming convention'),
  make_option(c('-d', '--out_dir'), type="character", help='Directory for output.')
)
opt <- parse_args(OptionParser(option_list = option_list))

print("Retrieving input files...")
net_obj <- readRDS(opt$net_obj)
outfile <- paste0(opt$out_dir, opt$out_name_project, "_network.tsv")

print("Converting interactome into a table and saving...")
InteractomeToTable(net_obj, outfile) 

print("Done.")
