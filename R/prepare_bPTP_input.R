
#' prepare_bPTP_input
#'
#' @description
#' This function processes a non-ultrametric tree and performs several checks to ensure its suitability for Bayesian Poisson Tree Processes (bPTP) analysis.
#' Specifically, it verifies that the input is a valid phylogenetic tree, ensures there are no duplicated tip labels, confirms that the tree has at least three tips, removes any outgroup sequences, checks for negative branch lengths, verifies that the tree is not ultrametric and untroots the tree if it is rooted.
#' Once all checks are successfully passed, the function creates a new directory named "bPTP" and saves the tree as "input_tree.newick" within it.
#' Finally, the function provides the user with instructions on how to run the analysis using the bPTP online tool, including guidance on which files to save and where to store them.
#' If the tree is large (>500 specimens), the function warns the user that bPTP may take a long time to run.
#'
#' @param phylogeny A non-ultrametric phylogenetic tree, as an object of class "phylo" read using ape::read.tree.
#' @param outgroup A string vector containing the name of the outgroup sequences, if they were included in the phylogenetic tree.
#'
#' @returns None. The function creates a newick file ("input_tree.newick") in the "bPTP" directory.
#' @export
#'
#' @examples
#' tree <- read.tree("iqtree.newick")
#' prepare_bPTP_input(phylogeny = tree, outgroup = c("outgroup_specimen_1", "outgroup_specimen_2", "outgroup_specimen_3"))
#'
prepare_bPTP_input <- function(phylogeny, outgroup = NULL) {

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
  if(ape::is.ultrametric(phylogeny) == TRUE) {
    stop("ERROR: The tree is ultrametric. bPTP needs a non-ultrametric tree as input")
  }

  # Unroot the tree if it was rooted
  if (!is.null(phylogeny$root.edge)) {
    message("The tree appears to be rooted. The tree will be unrooted automatically.")
    phylogeny <- unroot(phylogeny)
  }

  # Ensure the directory exists before writing the file
  if (!dir.exists("bPTP")) {
    dir.create("bPTP")
  }

  # Save the filtered alignment
  write.tree(phylogeny, file = "bPTP/input_tree.newick")

  # Tell the user where the new file is
  message(paste0("The file 'input_tree.newick' has been created at '", getwd(), "/bPTP'"))

  # Warn if the tree is too large
  if (length(phylogeny$tip.label) > 500) {
    message("WARNING: The tree contains more than 500 taxa. bPTP may take a considerable amount of time to run with a tree this large..")
  }  else {
    message("All checks have been passed. bPTP should run with no issue.")
  }

  # Give further instructions
  cat("\n")
  cat("To run bPTP:\n")
  cat("\n")
  cat(" - Visit the bPTP online tool at https://species.h-its.org/ptp/\n")
  cat(" - Upload the 'input_tree.newick' file.\n")
  cat(" - Indicate that the tree is unrooted.\n")
  cat(" - Set the number of MCMC generations to 500000\n")
  cat(" - Keep thinning at 100\n")
  cat(" - Set the Burn-in to 0.3.\n")
  cat(" - Keep the seed at 123.\n")
  cat(" - Leave the Outgroup taxa field blank.\n")
  cat(" - (Optional) Enter your email to receive a notification when the analysis is complete.\n")
  cat("\n")
  cat("After the analysis:\n")
  cat("\n")
  cat(" - Go to the 'Highest Bayesian supported solution' section.\n")
  cat(" - Click on 'Download delimitation results'.\n")
  cat(" - Save the file to: ", getwd(), "/bPTP.\n")
  cat(" - Keep the information at the bottom (e.g., Estimated number of species, mean number of species, accptance rate) for future reference.\n")
  cat("\n")

}
