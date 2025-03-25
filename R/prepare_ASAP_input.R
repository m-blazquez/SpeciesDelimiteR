
#' prepare_ASAP_input
#'
#' @description
#' This function processes a DNA alignment to prepare it as input for ASAP (Assemble Species by Automatic Partitioning).
#' It performs several checks to ensure the data is suitable for analysis: (1) verifies that the input is a valid DNA alignment, (2) removes outgroup sequences that could bias the ASAP analysis, (3) checks that all sequences have the same length, (4) identifies any empty sequences and (5) detects invalid characters in the sequences.
#' Once these checks are successfully passed, the function creates a new "ASAP" directory and saves the alignment as "input_alignment.fasta" within it.
#' The function then calculates the transition/transversion (ts/tv) ratio required to run ASAP using the Kimura (K80) model to compute distances, which we recommend as the default method. Finally, the function directs the user to the ASAP online tool, providing concise instructions on how to run the analysis, what files to save and where to save them.
#'
#' @param DNA_alignment A DNA alignment read using seqinr::read.fasta.
#' @param outgroup A string vector containing the name of the outgroup sequences, if they were included in the alignment.
#'
#' @returns None. The function creates a fasta file ("input_alignment.fasta") in the "ASAP" directory.
#' @export
#'
#' @examples
#' alignment <- seqinr:read.fasta("DNA_alignment.fasta", seqtype = "DNA", as.string = TRUE)
#' prepare_ASAP_input(DNA_alignment = alignment, outgroup = c("outgroup_specimen_1", "outgroup_specimen_2", "outgroup_specimen_3"))
#'
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

  # Ensure all sequences have the same length
  seq_lengths <- sapply(DNA_alignment, nchar)
  if (length(unique(seq_lengths)) > 1) {
    stop("ERROR: Sequences have different lengths. Ensure they are properly aligned.")
  }

  # Check for empty sequences
  empty_seqs <- names(DNA_alignment)[seq_lengths == 0]
  if (length(empty_seqs) > 0) {
    stop(paste("ERROR: The following sequences are empty:", paste(empty_seqs, collapse = ", ")))
  }

  # Define valid nucleotide characters (IUPAC codes + standard bases + gap)
  valid_chars <- c("A", "T", "G", "C", "-", "Y", "R", "W", "S", "M", "D", "V", "H", "N")

  # Check for invalid characters in the alignment (case-insensitive)
  invalid_seqs <- sapply(DNA_alignment, function(seq) {
    seq_upper <- toupper(seq)  # Convert sequence to uppercase
    any(!strsplit(seq_upper, "")[[1]] %in% valid_chars)
  })

  # Stop execution if any sequence contains unexpected characters
  if (any(invalid_seqs)) {
    stop(paste("ERROR: The following sequences contain invalid characters:",
               paste(names(DNA_alignment)[invalid_seqs], collapse = ", ")))
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

  # Calculate Ti/Tv ratio, avoiding division by zero
  if (transversions == 0) {
    TiTv_ratio <- NA
    warning("WARNING: No transversions found, Ts/Tv ratio is undefined.")
  } else {
    TiTv_ratio <- transitions / transversions
  }

  # Give further instructions
  message("All checks have been passed. ASAP should run with no issue.")

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
