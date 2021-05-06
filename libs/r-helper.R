
#' Generates a .tsv from an ARACNe network file.
#' @param net.obj Interactome object.
#' @param out.file File path to write out as .tsv.
#' @return If no outfile, returns the network as a data frame.
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
  } else {
    return(reg.df)
  }
}


### random code snippets, ignore below ###



InteractomeToTable(net.obj, 'C://Users/lvlah/linux/ac_lab/test-net.tsv')

net.obj <- readRDS('C://Users/lvlah/linux/ac_lab/GTEx-Networks/pruned_networks/brain-anterior-cingulate-cortex-ba24_pruned.rds')
ges.mat <- readRDS('C://Users/lvlah/linux/ac_lab/data/GTEx/subTissue_tpm/GTEx_artery-aorta_tpm.rds')
ges.mat <- GESTransform(ges.mat)

library(viper)

write.table(ges.mat, file = 'C://Users/lvlah/linux/ac_lab/data/pyther_data/gtex-aorta_ges.tsv',
            sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)


### MOR and weight matrices
targets <- unique(unlist(lapply(net.obj[1:3], function(x) names(x$tfmode)), use.names = FALSE))
# 1.2 Create the Mode of Regulation matrix from the regulon object.
mor <- sapply(net.obj[1:3], function(x, genes) {
  return(x$tfmode[match(genes, names(x$tfmode))])
}, genes = targets)
wts <- sapply(net.obj[1:3], function(x, genes) {
  tmp <- x$likelihood[match(genes, names(x$tfmode))]
  tmp[is.na(match(genes, names(x$tfmode)))] <- NA
  return(tmp/max(tmp, na.rm = T))
}, genes = targets)
mor[is.na(mor)] <- 0
wts[is.na(wts)] <- 0
wtss <- scale(wts, center = FALSE, scale = colSums(wts))

wts.vec <- net.obj$ENSG00000000003$likelihood
wts.vec <- wts.vec / max(wts.vec)
wts.vec <- wts.vec / sum(wts.vec)


### GES adjustment
t2.obj <- apply(ges.mat[,1:2], 2, rank) / (nrow(ges.mat) + 1)
t2q.obj <- qnorm(t2.obj)

t1.obj <- abs(t2.obj - 0.5) * 2
t1.obj <- t1.obj + (1 - max(t1.obj)) / 2
t1q.obj <- qnorm(t1.obj)

pos <- match(names(net.obj$ENSG00000000003$tfmode), names(t2.obj))

## 2-tail
t2.es <- t(mor*wtss) %*% t2q.obj[targets,]

## 1-tail
t1.es <- t((1 - abs(mor)) * wtss) %*% t1q.obj[targets,]

## integrate
ss <- sign(t2.es)
ss[ss == 0] <- 1
tot.es <- (abs(t2.es) + t1.es * (t1.es > 0)) * ss

## make NES
lwt <- sqrt(colSums(wts**2))
nes <- tot.es * lwt

viper(ges.mat[,1:2], net.obj[1:3], eset.filter = FALSE)

