
run_gblocks <- function(DNA_alignment_path = paste0(getwd(), "/Alignments/checked_initial_alignment.fasta"), mode = "default", exec = NULL) {

  # Read the FASTA file directly into a DNAbin object
  alignment_matrix <- read.dna(DNA_alignment_path, format = "fasta")

  # Check if the alignment has columns
  if (is.null(ncol(alignment_matrix)) || ncol(alignment_matrix) == 0) {
    stop("Error: Alignment matrix has no columns! Conversion failed.")
  }

  # Define Gblocks parameters based on mode
  params <- if (mode == "default") {
    list(b1 = 0.5, b2 = 0.5, b3 = ncol(alignment_matrix), b4 = 2, b5 = "a")
  } else if (mode == "less_stringent") {
    list(b1 = 0.5, b2 = 0.5, b3 = ncol(alignment_matrix), b4 = 5, b5 = "h")
  } else if (mode == "more_stringent") {
    list(b1 = 0.5, b2 = 0.5, b3 = 50, b4 = 2, b5 = "n")
  } else {
    stop("Invalid mode. Choose from 'default', 'less_stringent', or 'more_stringent'.")
  }

  # Check if the Gblocks executable is set
  if (is.null(exec)) {
    exec <- "Gblocks"
  }

  # Ensure the executable path is valid
  if (!file.exists(exec)) {
    stop("Error: Gblocks executable not found at '", exec, "'. Please provide the correct path.")
  }

  # Run Gblocks
  filtered_alignment <- gblocks(
    x = alignment_matrix,
    b1 = params$b1,
    b2 = params$b2,
    b3 = params$b3,
    b4 = params$b4,
    b5 = params$b5,
    target = "alignment",
    exec = exec)

  # Ensure the Alignments directory exists before writing the file
  if (!dir.exists("Alignments")) {
    dir.create("Alignments")
  }

  # Save the filtered alignment
  write.dna(filtered_alignment,
            file = "Alignments/gblocks_filtered_alignment.fasta",
            format = "fasta",
            colw = 80)

  # Tell the user where the new file is
  message(paste0("The file 'gblocks_filtered_alignment.fasta' has been created at '", getwd(), "/Alignments'"))
}
