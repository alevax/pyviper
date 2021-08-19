aREA <- function (eset, regulon, method = c("auto", "matrix", 
                                            "loop"), minsize = 20, cores = 1, wm = NULL, verbose = FALSE) 
{
  method <- match.arg(method)
  if (is.null(ncol(eset))) 
    eset <- matrix(eset, length(eset), 1, dimnames = list(names(eset), 
                                                          NULL))
  if (minsize > 0) {
    regulon <- lapply(regulon, function(x, genes) {
      pos <- names(x$tfmode) %in% genes
      list(tfmode = x$tfmode[pos], likelihood = x$likelihood[pos])
    }, genes = rownames(eset))
    regulon <- regulon[sapply(regulon, function(x) length(x$tfmode)) >= 
                         minsize]
    if (length(regulon) == 0) 
      stop(paste("There are no regulons with size of at least ", 
                 minsize, " targets.", sep = ""), 
           call. = FALSE)
    class(regulon) <- "regulon"
  }
  targets <- unique(unlist(lapply(regulon, function(x) names(x$tfmode)), 
                           use.names = FALSE))
  if (length(which(is.na(eset))) > 0) 
    method <- "loop"
  if (method == "auto") {
    method <- "matrix"
    if (length(targets) > 1000) 
      method <- "loop"
  }
  if (length(wm) > 0 | length(which(is.na(eset))) > 0) 
    method <- "loop"
  switch(method, matrix = {
    mor <- sapply(regulon, function(x, genes) {
      return(x$tfmode[match(genes, names(x$tfmode))])
    }, genes = targets)
    wts <- sapply(regulon, function(x, genes) {
      tmp <- x$likelihood[match(genes, names(x$tfmode))]
      tmp[is.na(match(genes, names(x$tfmode)))] <- NA
      return(tmp/max(tmp, na.rm = T))
    }, genes = targets)
    mor[is.na(mor)] <- 0
    wts[is.na(wts)] <- 0
    nes <- sqrt(colSums(wts^2))
    wts <- scale(wts, center = FALSE, scale = colSums(wts))
    pos <- match(targets, rownames(eset))
    print(length(pos))
    t2 <- apply(eset, 2, rank)/(nrow(eset) + 1)
    t1 <- abs(t2 - 0.5) * 2
    #print(head(t1))
    t1 <- t1 + (1 - max(t1))/2
    #print(head(t1))
    t1 <- qnorm(filterRowMatrix(t1, pos))
   # print(head(t1))
    t2 <- qnorm(filterRowMatrix(t2, pos))
    print(t2[c('NTMT1', 'AGTRAP', 'CDC20'),])
    print(t1[c('NTMT1', 'AGTRAP', 'CDC20'),])
    sum1 <- t(mor * wts) %*% t2
    print('DES:')
    print(sum1[1:10,])
    print(tail(sum1, 10))
    #print(sum1[c('NTMT1', 'AGTRAP', 'CDC20'),])
    sum2 <- t((1 - abs(mor)) * wts) %*% t1
    print(head(sum2))
    #print(sum2[c('NTMT1', 'AGTRAP', 'CDC20'),])
    ss <- sign(sum1)
    ss[ss == 0] <- 1
    tmp <- (abs(sum1) + sum2 * (sum2 > 0)) * ss
    print(head(tmp))
    tmp <- list(es = tmp, nes = tmp * nes)
  }, loop = {
    t2 <- t(t(apply(eset, 2, rank, na.last = "keep"))/(colSums(!is.na(eset)) + 
                                                         1))
    t1 <- abs(t2 - 0.5) * 2
    t1 <- t(t(t1) + (1 - apply(t1, 2, max, na.rm = TRUE))/2)
    t1 <- qnorm(t1)
    t2 <- qnorm(t2)
    if (is.null(wm)) wm <- matrix(1, nrow(eset), ncol(eset), 
                                  dimnames = list(rownames(eset), colnames(eset))) else {
                                    if (is.null(ncol(wm))) wm <- matrix(wm, length(wm), 
                                                                        ncol(eset), dimnames = list(names(wm), colnames(eset)))
                                    wm <- filterRowMatrix(wm, match(rownames(eset), rownames(wm)))
                                    rownames(wm) <- rownames(eset)
                                    wm[is.na(wm)] <- 0
                                  }
    wm[is.na(t1)] <- 0
    t1[is.na(t1)] <- 0
    t2[is.na(t2)] <- 0
    pb <- NULL
    if (cores > 1) {
      temp <- mclapply(1:length(regulon), function(i, regulon, 
                                                   t1, t2, ws) {
        x <- regulon[[i]]
        pos <- match(names(x$tfmode), rownames(t1))
        sum1 <- matrix(x$tfmode * x$likelihood, 1, length(x$tfmode)) %*% 
          (filterRowMatrix(t2, pos) * filterRowMatrix(ws, 
                                                      pos))
        ss <- sign(sum1)
        ss[ss == 0] <- 1
        sum2 <- matrix((1 - abs(x$tfmode)) * x$likelihood, 
                       1, length(x$tfmode)) %*% (filterRowMatrix(t1, 
                                                                 pos) * filterRowMatrix(ws, pos))
        return(as.vector(abs(sum1) + sum2 * (sum2 > 0))/colSums(x$likelihood * 
                                                                  filterRowMatrix(ws, pos)) * ss)
      }, regulon = regulon, t1 = t1, t2 = t2, mc.cores = cores, 
      ws = wm)
      temp <- sapply(temp, function(x) x)
    } else {
      if (verbose) {
        pb <- txtProgressBar(max = length(regulon), style = 3)
      }
      temp <- sapply(1:length(regulon), function(i, regulon, 
                                                 t1, t2, pb, ws) {
        x <- regulon[[i]]
        pos <- match(names(x$tfmode), rownames(t1))
        sum1 <- matrix(x$tfmode * x$likelihood, 1, length(x$tfmode)) %*% 
          (filterRowMatrix(t2, pos) * filterRowMatrix(ws, 
                                                      pos))
        ss <- sign(sum1)
        ss[ss == 0] <- 1
        sum2 <- matrix((1 - abs(x$tfmode)) * x$likelihood, 
                       1, length(x$tfmode)) %*% (filterRowMatrix(t1, 
                                                                 pos) * filterRowMatrix(ws, pos))
        if (is(pb, "txtProgressBar")) setTxtProgressBar(pb, 
                                                        i)
        return(as.vector(abs(sum1) + sum2 * (sum2 > 0))/colSums(x$likelihood * 
                                                                  filterRowMatrix(ws, pos)) * ss)
      }, regulon = regulon, t1 = t1, t2 = t2, pb = pb, 
      ws = wm)
    }
    if (is.null(ncol(temp))) temp <- matrix(temp, 1, length(temp))
    colnames(temp) <- names(regulon)
    rownames(temp) <- colnames(eset)
    if (length(which(wm < 1)) > 0) {
      w <- sapply(regulon, function(x, ws) {
        tmp <- x$likelihood * filterRowMatrix(ws, match(names(x$tfmode), 
                                                        rownames(ws)))
        sqrt(colSums(apply(tmp, 2, function(x) x/max(x))^2))
      }, ws = wm)
      if (is.null(ncol(w))) w <- t(w)
      w <- t(w)
    } else {
      w <- sapply(regulon, function(x) sqrt(sum((x$likelihood/max(x$likelihood))^2)))
    }
    return(list(es = t(temp), nes = t(temp) * w))
  })
  return(tmp)
}
