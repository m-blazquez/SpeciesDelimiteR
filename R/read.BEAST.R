
#' read.BEAST
#'
#' @description
#' This function expands on the ape::read.nexus by incorporating the posterior probabilties of a phylogenetic tree generated in BEAST.
#' The posterior probabilities are stored as node labels in the returned "phylo" object.
#'
#' @param file The file path to a nexus tree generated from a BEAST tree population using TreeAnnotator.
#'
#' @returns An object of class "phylo" with the posterior probabilities assigned to the node labels.
#' @export
#'
#' @examples
#' # Example usage of the function
#' tree <- read.BEAST(file = "BEAST.nexus")
#'
read.BEAST <- function(file) {
  # Read the tree from the NEXUS file
  tree <- ape::read.nexus(file)

  # Read file content as text
  nexus_text <- readLines(file)

  # Extract tree line
  tree_line <- grep("tree ", nexus_text, value = TRUE, ignore.case = TRUE)
  tree_string <- sub(".*= \\[&R\\] ", "", tree_line)  # Remove header

  # Extract posterior values from square brackets
  posterior_values <- regmatches(tree_string, gregexpr("posterior=([0-9\\.eE+-]+)", tree_string))[[1]]
  posterior_values <- as.numeric(sub("posterior=", "", posterior_values))

  # Assign posterior values to internal nodes
  if (length(posterior_values) == tree$Nnode) {
    tree$node.label <- as.character(posterior_values)
  } else {
    warning("Mismatch between extracted posterior values and internal nodes")
  }

  return(tree)
}
