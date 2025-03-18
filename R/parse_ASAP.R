
parse_ASAP <- function(ASAP_output, method = NULL) {

  # Make sure a valid method is specified
  if (is.null(method) || !method %in% c("JC69", "K80", "p_distances")) {

    stop("ERROR: Invalid method. Please specify either 'JC69', 'K80' or 'p_distances'")

  }

  # Assign new column names
  colnames(ASAP_output) <- c("Specimen", "ASAP_K80_species")

  # Write new file
  write.table(ASAP_output,
              file = "ASAP/ASAP_parsed_output.txt",
              sep = "\t",
              row.names = FALSE)

  # Tell the user where the new file is
  message(paste0("The file 'ASAP_parsed_output.txt' has been created at '", getwd(), "/ASAP'"))

  # Return the dataframe
  return(ASAP_output)

}
