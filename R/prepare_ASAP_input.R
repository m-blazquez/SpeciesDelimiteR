
prepare_ASAP_input <- function(DNA_alignment, outgroup = NULL) {

  # Check if DNA_alignment is a list (as read by read.fasta)
  if (!is.list(DNA_alignment)) {
    stop("DNA_alignment must be a list of sequences, like the output of read.fasta()")
  }

  # Remove outgroup sequences if provided
  if (!is.null(outgroup) && length(outgroup) > 0) {

    # Filter the alignment list to exclude outgroup sequences
    DNA_alignment <- DNA_alignment[!(names(DNA_alignment) %in% outgroup)]

    # Check if any outgroup sequence was not found in the alignment
    missing_outgroup <- outgroup[!(outgroup %in% names(DNA_alignment))]

    if (length(missing_outgroup) > 0) {
      warning(paste("WARNING: The following outgroup sequences were not found in the alignment:", paste(missing_outgroup, collapse = ", ")))
    }

  } else {
    message("No outgroup provided or outgroup vector is empty. Proceeding with the full alignment.")
  }

  # Rename remaining sequences by adding "DNA_" as a prefix (if they are only numbers ASAP does not run)
  new_names <- paste0("DNA_", names(DNA_alignment))
  names(DNA_alignment) <- new_names

  # Ensure the ASAP directory exists before writing the file
  if (!dir.exists("ASAP")) {
    dir.create("ASAP")
  }

  # Save the filtered alignment
  write.fasta(sequences = DNA_alignment,
              names = names(DNA_alignment),
              file.out = "ASAP/input_alignment.fasta")

  # Tell the user where the new file is
  message(paste0("The file 'input_alignment.fasta' has been created at '", getwd(), "/ASAP'"))

  # Initialize counters for transitions and transversions
  transitions <- 0
  transversions <- 0

  # Get the number of sequences and alignment length (columns)
  num_sequences <- length(DNA_alignment)
  alignment_length <- nchar(DNA_alignment[[1]])  # Length of the alignment (positions)

  # Loop through each column (position) in the alignment
  for (i in 1:alignment_length) {
    column <- toupper(sapply(DNA_alignment, function(seq) substr(seq, i, i)))  # Convert to uppercase

    # Get all pairwise comparisons of sequences at this position
    for (j in 1:(num_sequences - 1)) {
      for (k in (j + 1):num_sequences) {
        base1 <- column[j]
        base2 <- column[k]

        # Skip gaps (if any)
        if (base1 != "-" && base2 != "-") {
          # Check if it's a transition or transversion
          if ((base1 == "A" && base2 == "G") || (base1 == "G" && base2 == "A") ||
              (base1 == "C" && base2 == "T") || (base1 == "T" && base2 == "C")) {
            transitions <- transitions + 1
          } else if ((base1 == "A" && base2 == "C") || (base1 == "A" && base2 == "T") ||
                     (base1 == "G" && base2 == "C") || (base1 == "G" && base2 == "T") ||
                     (base1 == "C" && base2 == "A") || (base1 == "T" && base2 == "A") ||
                     (base1 == "C" && base2 == "G") || (base1 == "T" && base2 == "G")) {
            transversions <- transversions + 1
          }
        }
      }
    }
  }

  # Calculate Ti/Tv ratio
  TiTv_ratio <- transitions / transversions

  # Output the Ti/Tv ratio
  cat("\n")
  cat("Transition/Transversion ratio:", TiTv_ratio, "\n")
  cat("\n")
  cat("To run ASAP:\n")
  cat("\n")
  cat(" - Visit the ASAP online tool at https://bioinfo.mnhn.fr/abi/public/asap/\n")
  cat(" - Upload the 'input_alignment.fasta' file.\n")
  cat(" - Choose the Kimura (K80) model for distance computation.\n")
  cat(" - Enter the ts/tv ratio: ", TiTv_ratio, "\n")
  cat(" - Click 'Go' to start the analysis.\n")
  cat("\n")
  cat("After the analysis:\n")
  cat("\n")
  cat(" - Identify the partition with the lowest ASAP score.\n")
  cat(" - Download the partition file as CSV under the 'Text' option.\n")
  cat(" - Save the CSV file to: ", getwd(), "/ASAP.\n")
  cat(" - Save the entire webpage for future reference. \n")
  cat("\n")

}
