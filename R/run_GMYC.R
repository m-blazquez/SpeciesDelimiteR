
run_GMYC <- function(phylogeny, outgroup = NULL, method = NULL) {

  # Make sure that the tree is ultrametric. If not, return an error message.
  if(ape::is.ultrametric(phylogeny) == FALSE) {

    stop("ERROR: The tree is not ultrametric.")

  }

  # Ensure a valid method is specified
  if (is.null(method) || !method %in% c("single", "multiple")) {

    stop("ERROR: You must specify a valid method ('single' or 'multiple').")

  }


  # Remove outgroup sequences (if provided)
  if (!is.null(outgroup) && length(outgroup) > 0) {

    for (i in seq_along(outgroup)) {

      if (outgroup[i] %in% phylogeny$tip.label) {

        phylogeny <- drop.tip(phylogeny, outgroup[i])

      } else {

        stop(paste("ERROR: Outgroup", outgroup[i], "not found in the phylogeny."))

      }
    }
  } else {

    message("No outgroup provided or outgroup vector is empty. Proceeding with the full tree.")

    }

  # Runs GMYC using the splits package, using the single- or multiple- threshold options.
  if(method == "single") {

    GMYC <- splits::gmyc(phylogeny, method = "single")
    GMYC <- spec.list(GMYC)
    return(GMYC)

  }

  if(method == "multiple") {

    GMYC <- splits::gmyc(phylogeny, method = "multiple")
    GMYC <- spec.list(GMYC)
    return(GMYC)

  }

  # Save the putative species file
  if(method == "single") {

    write.table(GMYC,
                file = "GMYCs_output.txt",
                sep = "\t",
                row.names = FALSE)

    # Tell the user where the new file is
    message(paste0("The file 'GMYCs_output.txt' has been created at '", getwd(), "'"))

  }

  if(method == "multiple") {

    write.table(GMYC,
                file = "GMYCm_output.txt",
                sep = "\t",
                row.names = FALSE)

    # Tell the user where the new file is
    message(paste0("The file 'GMYCm_output.txt' has been created at '", getwd(), "'"))

  }

}

