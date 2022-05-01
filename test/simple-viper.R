
# 1.1 Load VIPER library
library(viper) # 1.2 Load data
data(bcellViper, package="bcellViper",verbose = TRUE)
# 1.3 Get expression matrix
eset <- exprs(dset)
str(regulon,1)

# Step 1 - Filter the expression data.
# 1.1 Get a list of all targets of every regulator in the network.
all_targets <- c(names(regulon), unlist(lapply(regulon, function(x) names(x$tfmode)), use.name = FALSE ))
# 1.2 Remove expression data of all the genes which aren't in the list identified abov
eset <- eset[rownames(eset) %in% unique(all_targets),]
# 1.3 Scale the expression matrix on rows [i.e. row = (row - mean(row)) / sd(row) ].
tt <- t(scale(t(eset)))

# Step 2 - Filter the regulon object (i.e. the regulatory network).
# 2.1 Remove targets in the regulon which are not in the rownames of the expression matrix
regulon <- lapply(regulon, function(x, genes) {
	filtro <- names(x$tfmode) %in% genes 
	x$tfmode <- x$tfmode[filtro]
	if (length(x$likelihood) == length(filtro))
		x$likelihood <- x$likelihood[filtro] 
	return(x)
}, genes = rownames(eset))
# 2.2 Define minimum regulon size for filtering (default is 20).
# The 'minsize' parameter is specified in the function parameters.
# 2.3 Remove regulators with a regulon size below the 'minsize' parameter.
minsize <- 50
regulon <- regulon[sapply(regulon, function(x) length(x$tfmode)) >= minsize] # aREA
length(regulon)
# nes <- simpleaREA(tt, regulon)

regulon <- pruneRegulon(regulon,cutoff = 50,adaptive = F,eliminate = T)

simpleaREA <- function (tt, regulon)
{
	## Step 1 - Create the 'Mode of Regulation' and 'Weights' matrices.
	# 1.1 Get a list of all targets of every regulator in the network.
	targets <- unique(unlist(lapply(regulon, function(x) names(x$tfmode)), use.names = FALSE) )
	# 1.2 Create the Mode of Regulation matrix from the regulon object.
	mor <- sapply(regulon, function(x, genes) { return(x$tfmode[match(genes, names(x$tfmode))])
	}, genes = targets)
	
	dim(mor)
	tmp <- ifelse(is.na(mor),0,mor)
	colSums(tmp != 0)
	
	dim(mor)
	head(rownames(mor),20)
	head(targets,20)
	tail(rownames(mor),20)
	tail(targets,20)	
	length(targets)
	
	# 1.2 Create the Weights matrix from the regulon object.
	wts <- sapply(regulon, function(x, genes) {
		tmp <- x$likelihood[match(genes, names(x$tfmode))]
		tmp[is.na(match(genes, names(x$tfmode)))] <- NA
		return(tmp/max(tmp, na.rm = T)) }, genes = targets)
	
	tmp <- ifelse(is.na(wts),0,wts)
	colSums(tmp != 0)
	
	# 1.3 For each regulator, assign values of 0 to genes which are not listed as its targets
	mor[is.na(mor)] <- 0 
	wts[is.na(wts)] <- 0
	# 1.4 Scale the columns of the 'Weights' matrix to the sum of the weights.
	wtss <- scale(wts, center = FALSE, scale = colSums(wts))
	
	## Step 2 - Calculate the two-tail enrichment scores.
	# 2.1 Calculate the 'T2 rank' matrix from the expression dataset. 
	t1 <- apply(tt, 2, rank)/(nrow(tt) + 1)
	pos <- match(targets, rownames(tt)) # This line of code is necessary to match the order of genes # for the subsequent matrix multiplication steps.
	# 2.2 Transform T2 ranks to Gaussian values.
	t2 <- filterRowMatrix(t1, pos)
	dim(t1)
	dim(t2)
	t2q <- qnorm(t2) 
	dim(t2q)
	
	#genewise
	plot(density(t2[1,]))
	plot(density(t2q[1,]))
	
	#samplewise
	plot(density(t2[,1]))
	plot(density(t2q[,1]))
	
	tmp[1:5,1:5] # t2
	t2q[1:5,1:5]
	
	rowSums(t2[1:5,])
	rowSums(t2q[1:5,])
	# colSums(t2[,1:5])
	# colSums(t2q[,1:5])
	plot(density(eset[1,]))
	plot(density(tt[1,]))
	plot(density(t1[1,]))
	plot(density(t2[1,]))
	plot(density(t2q[1,]))
	
	
	# 2.3 Matrix multiplication.
	dim(mor)
	dim(wtss)
	dim(t2q)
	
	sum(!is.na(rownames(mor)))
	interactome_matrix <- t(mor * wtss) # mode of action * likelihood
	colSums(tmp != 0)
	
	ges_matrix <- t2q
	
	dim(interactome_matrix)
	dim(ges_matrix)
	sum1 <- interactome_matrix %*% ges_matrix
	sum(!is.na(colnames(interactome_matrix)))
	sum(!is.na(rownames(interactome_matrix)))
	dim(sum1)
	
	## Step 3 - Calculate the one-tail score matrix.
	# 3.1 Calculate the 'T1 Rank' matrix from the 'T2 Rank' matrix. 
	t1 <- abs(t2 - 0.5) * 2
	t1 <- t1 + (1 - max(t1))/2
	# 3.2 Get qnorm values
	pos <- match(targets, rownames(t1))
	t1q <- qnorm(filterRowMatrix(x = t1, filter = pos))
	# 3.3 Matrix multiplication.
	sum2 <- t((1 - abs(mor)) * wtss) %*% t1q
	
	## Step 4 - Calculate the Three-tail enrichment score.
	# 4.1 Extract the signs of the Two-tail enrichment scores 
	ss <- sign(sum1)
	ss[ss == 0] <- 1
	# 4.2 Combine the Two-tail and One-tail enrichment score matrices.
	sum3 <- (abs(sum1) + sum2 * (sum2 > 0)) * ss
	# Step 5 - Calculate the Normalized Enrichment Scores.
	# 5.1 For each regulator, calculate an index proportional to the likelihood value # of all its regulatory interactions.
	lwt <- sqrt(colSums(wts^2))
	# 5.2 Adjust 3T enrichment scores proportionally to the weights calculated above.
	nes <- sum3 * lwt
}

x <- simpleaREA(tt,regulon)

# aREA_single(ges = ges, regulon = regulon$MYB)
y <- aREA_single(ges = tt[,1], regulon = regulon$MYB)

x["MYB",1:5]
y$nes[1:5]

