##
# Test Functions for Enrichment Analysis in R
# -------------------------------------------
# source("../vaxtools/R/gsea-plot.R")


library(viper)

aREA_test_function <- function(ges, regulon, wm = NULL){
        
        # Make sure gene expression is in matrix format
        if (is.null(ncol(ges))) {
                ges <- matrix(ges, length(ges), 1, dimnames = list(names(ges), NULL))
        }
        
        # t2, t1 and weight matrices
        t2 <- t(t(apply(ges, 2, rank, na.last = "keep"))/(colSums(!is.na(ges)) + 1))
        t1 <- abs(t2 - 0.5) * 2
        t1 <- t(t(t1) + (1 - apply(t1, 2, max, na.rm = TRUE))/2)
        t1 <- qnorm(t1)
        t2 <- qnorm(t2)
        
        if (is.null(wm)) {
                wm <- matrix(1, nrow(ges), ncol(ges),
                             dimnames = list(rownames(ges), colnames(ges)))
        }
        
        # prepare regulon, prune and match
        x <- regulon
        common <- intersect(names(x$tfmode), rownames(ges))
        if(length(common) == 0) stop('None of the targets are in the gene expression matrix!', call. = FALSE)
        keep <- names(x$tfmode) %in% common
        x <- lapply(x, function(i) i[keep])
        pos <- match(names(x$tfmode), rownames(t1))
        
        # Calculations with t2 and t1 matrices
        sum1 <- matrix(x$tfmode * x$likelihood, 1, length(x$tfmode)) %*% (t2[pos, ] * wm[pos, ])
        ss <- sign(sum1)
        if (ss == 0) ss <- 1
        sum2 <- matrix((1 - abs(x$tfmode)) * x$likelihood, 1, length(x$tfmode)) %*% (t1[pos, ] * wm[pos, ])
        
        # derive ES, NES and p-value
        es <- (abs(sum1) + sum2 * (sum2 > 0))/sum(x$likelihood * wm[pos, ]) * ss
        nes <- es * sqrt(sum((x$likelihood/max(x$likelihood))^2))
        pval <- stats::pnorm(abs(nes), lower.tail = FALSE)*2
        lapply(list(es = es, nes = nes, pval = pval), as.vector)
        
}

