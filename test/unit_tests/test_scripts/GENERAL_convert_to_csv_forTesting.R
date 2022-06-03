## Libraries
print("Running GENERAL_convert_to_csv_forTesting.R...")
print("Loading optparse, tidyverse and stringr libraries...")
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(stringr)))

## Functions
parseOptCharacter <- function(opt_char, default = NULL){
  if(testMissing(opt_char)){
    opt_char <- default
  }
  opt_char
}
testMissing <- function(opt_value){
  length(opt_value) == 0 || is.na(opt_value) || is.null(opt_value) || opt_value == ""
}
getInputFileFolder <- function(mat_file_path){
  require(tidyverse)
  file_folder <- mat_file_path %>%
    str_split("[/]") %>%
    unlist() %>%
    head(-1) %>%
    paste(collapse = "/") %>%
    paste("/", sep="")
  file_folder
}
getInputFileName <- function(mat_file_path, ext = FALSE){
  require(tidyverse)
  pattern <- ifelse(ext, "/", "[/\\.]")
  n <- ifelse(ext, -1, -2)
  file_name <- mat_file_path %>%
    str_split(pattern) %>%
    unlist() %>%
    nth(n)
  file_name
}

## Arguments
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-f', '--file'), type="character", help='RDS file to make a CSV'),
  make_option(c('-e', '--ext'), type="character", help='csv, tsv, or other (default = csv)'),
  make_option(c('-o', '--out_dir'), type="character", help='Out directory for file. Default is file location')
)

opt <- parse_args(OptionParser(option_list = option_list))
file_ext <- tolower(parseOptCharacter(opt$ext, default = "csv"))
out_dir <- parseOptCharacter(opt$out_dir, default =  getInputFileFolder(opt$file))
if(file_ext == "tsv"){
  sep_char = "\t"
} else if(file_ext == "txt"){
  sep = " "
} else{ #if(file_ext == "csv"){
  sep_char = ","
}

print("Retrieving input RDS file...")
file <- readRDS(opt$file)
outfile <- paste0(out_dir,
                  getInputFileName(opt$file),
                  ".",
                  file_ext)
print(paste0("Saving as ", toupper(file_ext), "..."))
write.table(file, file = outfile, sep = sep_char)
print("Done.")
