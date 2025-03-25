
check_initial_alignment <- function(DNA_alignment, outgroup = NULL) {

  # Check if DNA_alignment is a list (as read by read.fasta)
  if (!is.list(DNA_alignment)) {
    stop("DNA_alignment must be a list of sequences, like the output of read.fasta()")
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

  # Check for long gaps (>20 bp)
  long_gap_seqs <- names(DNA_alignment)[sapply(DNA_alignment, function(seq) {
    any(grepl("-{20,}", seq))  # Regex detects gaps of length 20 or more
  })]

  if (length(long_gap_seqs) > 0) {
    warning(paste("WARNING: The following sequences contain long gaps (>20 bp):",
                  paste(long_gap_seqs, collapse = ", ")))
  }

  # Check for sequences with >25% missing data (Ns)
  high_N_seqs <- names(DNA_alignment)[sapply(DNA_alignment, function(seq) {
    seq_chars <- unlist(strsplit(seq, ""))
    (sum(seq_chars == "N") / length(seq_chars)) > 0.25  # Calculate N proportion
  })]

  if (length(high_N_seqs) > 0) {
    warning(paste("WARNING: The following sequences contain more than 25% missing data (Ns): ",
                  paste(high_N_seqs, collapse = ", ")))
  }

  # Check for highly divergent sequences

  # Convert the list of SeqFastadna objects to a list of character sequences
  seq_strings <- lapply(DNA_alignment, as.character)

  # Convert to DNAbin format and preserve names
  alignment_matrix <- as.DNAbin(seq_strings)

  # Compute pairwise distances
  dist_matrix <- as.matrix(dist.dna(alignment_matrix, model = "raw"))

  # Ensure row and column names are assigned
  rownames(dist_matrix) <- names(DNA_alignment)
  colnames(dist_matrix) <- names(DNA_alignment)

  # Get mean distances per sequence
  mean_distances <- rowMeans(dist_matrix)

  # Ensure names are properly assigned to mean distances
  names(mean_distances) <- names(DNA_alignment)

  # Define threshold for divergence (mean + 2SD)
  threshold <- mean(mean_distances, na.rm = TRUE) + 2 * sd(mean_distances, na.rm = TRUE)

  # Identify divergent sequences (excluding outgroup)
  divergent_seqs <- names(mean_distances)[mean_distances > threshold &
                                            !(names(mean_distances) %in% outgroup)]

  # Check for valid names before issuing warning
  if (length(na.omit(divergent_seqs)) > 0) {
    warning(paste("WARNING: The following sequences are highly divergent:",
                  paste(na.omit(divergent_seqs), collapse = ", ")))
  }

  # Ensure the Alignments directory exists before writing the file
  if (!dir.exists("Alignments")) {
    dir.create("Alignments")
  }

  # Save the filtered alignment
  write.fasta(sequences = DNA_alignment,
              names = names(DNA_alignment),
              file.out = "Alignments/checked_initial_alignment.fasta")

  # Tell the user where the new file is
  message(paste0("The file 'checked_initial_alignment.fasta' has been created at '", getwd(), "/Alignments'"))

}
