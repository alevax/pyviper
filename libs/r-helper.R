
InteractomeToTable <- function(net.obj, out.file) {
  # make df
  net.df <- do.call(rbind, lapply(names(net.obj[1:3]), function(x) {
    reg.obj <- net.obj[[x]] 
    reg.df <- data.frame('regulator' = rep(x, length(reg.obj$tfmode)),
                         'target' = names(reg.obj$tfmode),
                         'mor' = reg.obj$tfmode,
                         'likelihood' = reg.obj$likelihood)
    rownames(reg.df) <- NULL
    return(reg.df)
  }))
  # if specified, write to file
  if (!missing(out.file)) {
    write.table(net.df, file = out.file, quote = FALSE, sep = '\t',
                row.names = FALSE, col.names = TRUE)
  }
}

InteractomeToTable(net.obj, 'C://Users/lvlah/linux/ac_lab/test-net.tsv')

ges.mat <- readRDS('C://Users/lvlah/linux/ac_lab/data/GTEx/subTissue_tpm/GTEx_artery-aorta_tpm.rds')
ges.mat <- 