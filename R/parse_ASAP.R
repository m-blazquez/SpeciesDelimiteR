
parse_ASAP <- function(ASAP_output, prepare_ASAP_input = NULL) {

  # Make sure a valid value for prepare_ASAP_input is specified
  if (is.null(prepare_ASAP_input) || !prepare_ASAP_input %in% c("JC69", "K80", "p_distances")) {

    stop("ERROR: Invalid value for prepare_ASAP_input argument. Please specify either TRUE or FALSE")

  }

  # Return to original specimen names if the prepare_ASAP_input function was used
  if (prepare_ASAP_input) {

    ASAP_output$V1 <- substr(ASAP_output$V1, 5, nchar(ASAP_output$V1))

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
