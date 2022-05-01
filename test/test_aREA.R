##
# Test Functions for Enrichment Analysis in R
# -------------------------------------------
source("../vaxtools/R/gsea-plot.R")

source("../pyther/libs/area_fn.R")
library(viper)

data(bcellViper, package = 'bcellViper')

d1 <- Biobase::exprs(dset)
res <- viper::rowTtest(dset, "description", "CB", "N")

ges <- res$statistic[,1]

gsr <- gsea_regulon(ges, regulon$MYB)
plot(gsr, plotSignature = TRUE,
     signatureNames = c('CB', 'N'),
     signatureType = 't-statistic')

aREA_single(ges = ges, regulon = regulon$MYB)


## My Area ----

my_area <- function(ges,interactome)
{
	create_mor_matrix <- function(interactome,target_genes) {
		target_genes <- sort(target_genes)
		mat <- sapply( interactome , function(x) { return( x$tfmode[ match(target_genes, names(x$tfmode)) ] ) })	
		mat[is.na(mat)] <- 0 
		rownames(mat) <- target_genes
		return(mat)
	}
	
	create_likelihood_matrix <- function(interactome,target_genes) {
		target_genes <- sort(target_genes)
		mat <- sapply( interactome , function(x) { 
			res <- x$likelihood[ match(target_genes, names(x$tfmode)) ] 
			res[is.na(match(target_genes, names(x$tfmode)))] <- 0
			res <- res/max(res,na.rm = T) # If the max is < 1, then it increases the actual score because its div by 0.x
			return(res)	
		})	
		rownames(mat) <- target_genes
		return(mat)
	}
	
	targets <- unique(unlist(lapply(interactome, function(x) names(x$tfmode)), use.names = FALSE) )
	targets <- sort(targets)
	
	## Create Mode Of Regulation Matrix ----
	mor_mat <- create_mor_matrix(interactome,targets)
	summary(colSums(mor_mat != 0))
	## Create Likelihood/weights Matrix ----
	lik_mat <- create_likelihood_matrix(interactome = interactome,target_genes = targets)
	summary(colSums(lik_mat != 0))
	
	stopifnot(identical(colnames(mor_mat),colnames(lik_mat)))
	stopifnot(identical(rownames(mor_mat),rownames(lik_mat)))
	
	lik_mat_scaled <- scale(lik_mat, center = FALSE, scale = colSums(lik_mat))
	str(lik_mat_scaled,1)
	
	## Step 2 - Calculate the two-tail enrichment scores.
	# 2.1 Calculate the 'T2 rank' matrix from the expression dataset. 
	t1 <- apply(ges, 2, rank)/(nrow(ges) + 1) # t1 scores should be contained in (0,1)
	pos <- match(targets, rownames(ges)) # This line of code is necessary to match the order of genes # for the subsequent matrix multiplication steps.
	# 2.2 Transform T2 ranks to Gaussian values.
	t2 <- filterRowMatrix(t1, pos)
	dim(t1)
	dim(t2)
	t2q <- qnorm(t2) # Converting uniform distribution to a normal distribution?!? Bah...
	dim(t2q)
	ges_matrix <- t2q
	dim(ges_matrix)
	
	interactome_matrix <- t(mor_mat * lik_mat_scaled) # mode of action * likelihood
	
	dim(interactome_matrix)
	dim(ges_matrix)
	
	sum1 <- interactome_matrix %*% ges_matrix
	dim(sum1)
	
	
	foxm1_targets <- interactome_matrix["FOXM1", interactome_matrix["FOXM1",] != 0 ]
	foxm1_targets <- names(foxm1_targets)
	ges_matrix[foxm1_targets,"GSM44170"]
	mean(ges_matrix[foxm1_targets,"GSM44170"]) # WHAAAAAT?!?!?!?
	sum1["FOXM1","GSM44170"] # WHAAAAAT?!?!?!? WHAAAAAT?!?!?!?
	
	## my little test on FOXM1 and Sample ID:1 ----
	x <- interactome_matrix["FOXM1",]
	x <- x[x!=0]
	y <- ges_matrix[,1]
	y <- y[names(x)]
	sum(x*y)
	# --
	sum1["FOXM1",1]
	
	## Step 3 - Calculate the one-tail score matrix.
	# 3.1 Calculate the 'T1 Rank' matrix from the 'T2 Rank' matrix. 
	t1 <- abs(t2 - 0.5) * 2 # This (t2) is not yet converted to a gaussian
	t1 <- t1 + (1 - max(t1))/2 # WHY this? it adds a very low constant to each value in the matrix !?!?
	# 3.2 Get qnorm values
	# pos <- match(targets, rownames(t1))
	# t1q <- qnorm(filterRowMatrix(x = t1, filter = pos))
	t1q <- qnorm(t1)
	# 3.3 Matrix multiplication.
	sum2 <- t((1 - abs(mor_mat)) * lik_mat_scaled) %*% t1q
	
	sum2["FOXM1",1]
	
	## Step 4 - Calculate the Three-tail enrichment score.
	# 4.1 Extract the signs of the Two-tail enrichment scores 
	ss <- sign(sum1)
	ss[ss == 0] <- 1
	# 4.2 Combine the Two-tail and One-tail enrichment score matrices.
	sum3 <- (abs(sum1) + sum2 * (sum2 > 0)) * ss # This should be the es, not yet nes
	
	sum3["FOXM1",1] 
	sum3["FOXM1",1] * lwt["FOXM1"] # This is the final score AKA nes
	
	# Step 5 - Calculate the Normalized Enrichment Scores.
	# 5.1 For each regulator, calculate an index proportional to the likelihood value # of all its regulatory interactions.
	lwt <- sqrt(colSums(lik_mat^2))
	# 5.2 Adjust 3T enrichment scores proportionally to the weights calculated above.
	nes <- sum3 * lwt
	
	return(nes)
	
}

library(viper) # 1.2 Load data
data(bcellViper, package="bcellViper",verbose = TRUE)
# 1.3 Get expression matrix
eset <- exprs(dset)
# str(regulon,1)

ges <- t(scale(t(eset)))
interactome <- pruneRegulon(regulon,cutoff = 50,adaptive = FALSE,eliminate = TRUE)

my_vp <- my_area(ges,interactome)
my_sa <- simpleaREA(ges,interactome)

identical( round(my_vp,3) , round(my_sa,3) )

my_vp[1:5,1:5]
my_sa[1:5,1:5]

my_vp["FOXM1",1] # sum3["FOXM1",1] * lwt["FOXM1"] # This is the final score AKA nes
