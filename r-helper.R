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


TableToInteractome <- function(net_table){
  blank_reg <- list(
    c(),
    c()
  )
  names(blank_reg) <- c("tfmode", "likelihood")
  my_reg <- list()
  pb = txtProgressBar(min = 0, max = nrow(net_table), initial = 0, style = 3)
  for(i in 1:nrow(net_table)){
    my_reg[[net_table[i,"regulator"]]] <- blank_reg
    setTxtProgressBar(pb, i)
  }
  close(pb)

  pb = txtProgressBar(min = 0, max = nrow(net_table), initial = 0, style = 3)
  for (i in 1:nrow(net_table)){
    new_w <- net_table[i,"likelihood"]
    new_m <- net_table[i,"mor"]
    names(new_w) <- net_table[i,"target"]
    names(new_m) <- net_table[i,"target"]
    my_reg[[net_table[i,"regulator"]]][["likelihood"]] <- c(
      my_reg[[net_table[i,"regulator"]]][["likelihood"]],
      new_w
    )
    my_reg[[net_table[i,"regulator"]]][["tfmode"]] <- c(
      my_reg[[net_table[i,"regulator"]]][["tfmode"]],
      new_m
    )
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(my_reg)
}
