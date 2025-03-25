
#' run_GMYC
#'
#' @description
#' This function takes an ultrametric tree and performs several checks to ensure its suitability for GMYC (Generalized Mixed Yule Coalescent) analysis.
#' Specifically, it verifies that the input is a valid phylogenetic tree, ensures there are no duplicated tip labels, confirms that the tree has at least three tips, removes any outgroup sequences, checks for negative branch lengths and verifies that the tree is ultrametric.
#' The function then acts as a wrapper for splits::gmyc, running the analysis using either the "single" or "multiple" method.
#' After the analysis is complete, it uses splits::spec.list to extract the putative species assignments for each specimen. The function creates a new directory (either "GMYC/Single" or "GMYC/Multiple") and saves the species list as a .txt file. Additionally, it returns the species list for further use.
#'
#' @param phylogeny An ultrametric phylogenetic tree, as an object of class "phylo" read using read.BEAST or ape::read.nexus.
#' @param outgroup A string vector containing the name of the outgroup sequences, if they were included in the phylogenetic tree.
#' @param method A string specifying the GMYC method to use. It must be either "single" for single-threshold analysis or "multiple" for multiple-threshold analysis.
#'
#' @returns A `data.frame` with columns `GMYC_spec` (the putative species assigned by GMYC) and `sample_name` (the name of each specimen).
#' @export
#'
#' @examples
#' ultrametric_tree <- read.BEAST(file = "BEAST.nexus")
#' GMYC_single <- run_GMYC(phylogeny = ultrametric_tree, outgroup = c("outgroup_specimen_1", "outgroup_specimen_2", "outgroup_specimen_3"), method = "single")
#' GMYC_multiple <- run_GMYC(phylogeny = ultrametric_tree, outgroup = c("outgroup_specimen_1", "outgroup_specimen_2", "outgroup_specimen_3"), method = "multiple")
#'
run_GMYC <- function(phylogeny, outgroup = NULL, method = NULL) {

  # Make sure that the input is a valid phylogenetic tree
  if (!inherits(phylogeny, "phylo")) {
    stop("ERROR: The provided input is not a valid phylogenetic tree ('phylo' object).")
  }

  # Ensure there are no duplicated tip labels
  if (any(duplicated(phylogeny$tip.label))) {
    stop("ERROR: The tree contains duplicate tip labels. Ensure that all taxa have unique names.")
  }

  # Check if the tree has at least 3 tips
  if (length(phylogeny$tip.label) < 3) {
    stop("ERROR: The tree must have at least three taxa for bPTP analysis.")
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

  # Check for negative branch lengths
  if (any(phylogeny$edge.length < 0)) {
    stop("ERROR: The tree contains negative branch lengths.")
  }

  # Make sure that the tree is not ultrametric
  if(ape::is.ultrametric(phylogeny) == FALSE) {
    stop("ERROR: The tree is not ultrametric.")
  }

  # Ensure a valid method is specified
  if (is.null(method) || !method %in% c("single", "multiple")) {
    stop("ERROR: You must specify a valid method ('single' or 'multiple').")
  }

    # Runs GMYC using the splits package, using the single- or multiple-threshold options.
  if (method == "single") {
    GMYC <- splits::gmyc(phylogeny, method = "single")
    GMYC <- splits::spec.list(GMYC)

    # Ensure the directory exists before writing the file
    if (!dir.exists("GMYC/Single")) {
      dir.create("GMYC/Single", recursive = TRUE)
    }

    # Write the file
    write.table(GMYC,
                file = "GMYC/Single/GMYCs_output.txt",
                sep = "\t",
                row.names = FALSE)

    # Tell the user where the new file is
    message(paste0("The file 'GMYCs_output.txt' has been created at '", getwd(), "/GMYC/Single'"))

    return(GMYC)  # Return after writing file
  }

  if (method == "multiple") {
    GMYC <- splits::gmyc(phylogeny, method = "multiple")
    GMYC <- splits::spec.list(GMYC)

    # Ensure the directory exists before writing the file
    if (!dir.exists("GMYC/Multiple")) {
      dir.create("GMYC/Multiple", recursive = TRUE)
    }

    # Write the file
    write.table(GMYC,
                file = "GMYC/Multiple/GMYCm_output.txt",
                sep = "\t",
                row.names = FALSE)

    # Tell the user where the new file is
    message(paste0("The file 'GMYCm_output.txt' has been created at '", getwd(), "/GMYC/Multiple'"))

    return(GMYC)  # Return after writing file
  }
}

