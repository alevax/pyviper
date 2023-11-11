#' Generates a .tsv from an ARACNe network file.
#' @param net.obj Interactome object.
#' @param out.file File path to write out as .tsv.
#' @return If no outfile, returns the network as a data frame.
InteractomeToTable <- function(net.obj, out.file, is_a3_network=FALSE) {
	# make df

	if( is_a3_network == TRUE) {
		net.df <- do.call(rbind, lapply(names(net.obj), function(x) {
			reg.obj <- net.obj[[x]]
			reg.df <- data.frame('regulator' = rep(x, length(reg.obj$am)),
													 'target' = names(reg.obj$am),
													 'mor' = reg.obj$am,
													 'likelihood' = reg.obj$aw)
			rownames(reg.df) <- NULL
			return(reg.df)
		}))
	} else {
		net.df <- do.call(rbind, lapply(names(net.obj), function(x) {
			reg.obj <- net.obj[[x]]
			reg.df <- data.frame('regulator' = rep(x, length(reg.obj$tfmode)),
													 'target' = names(reg.obj$tfmode),
													 'mor' = reg.obj$tfmode,
													 'likelihood' = reg.obj$likelihood)
			rownames(reg.df) <- NULL
			return(reg.df)
		}))
	}

	# if specified, write to file
	if (!missing(out.file)) {
		write.table(net.df, file = out.file, quote = FALSE, sep = '\t',
								row.names = FALSE, col.names = TRUE)
	} else {
		return(reg.df)
	}
}

