
parse_ASAP <- function(ASAP_output, method = NULL) {

  # Make sure a valid method is specified
  if (is.null(method) || !method %in% c("JC69", "K80", "p_distances")) {

    stop("ERROR: Invalid method. Please specify either 'JC69', 'K80' or 'p_distances'")

  }

  # Parse the ASAP output based on the method used to produce it
  if(method == "JC69") {

    # Assign new column names
    colnames(ASAP_output) <- c("Specimen", "ASAP_JC69_species")

    # Write new file
    write.table(ASAP_output,
                file = "ASAP_JC69_parsed_output.txt",
                sep = "\t",
                row.names = FALSE)

    # Tell the user where the new file is
    message(paste0("The file 'ASAP_JC69_parsed_output.txt' has been created at '", getwd(), "'"))

    # Return the dataframe
    return(ASAP_output)

  }

  if(method == "K80") {

    # Assign new column names
    colnames(ASAP_output) <- c("Specimen", "ASAP_K80_species")

    # Write new file
    write.table(ASAP_output,
                file = "ASAP_K80_parsed_output.txt",
                sep = "\t",
                row.names = FALSE)

    # Tell the user where the new file is
    message(paste0("The file 'ASAP_K80_parsed_output.txt' has been created at '", getwd(), "'"))

    # Return the dataframe
    return(ASAP_output)

  }

  if(method == "p_distances") {

    # Assign new column names
    colnames(ASAP_output) <- c("Specimen", "ASAP_p_distances_species")

    # Write new file
    write.table(ASAP_output,
                file = "ASAP_p_distances_parsed_output.txt",
                sep = "\t",
                row.names = FALSE)

    # Tell the user where the new file is
    message(paste0("The file 'ASAP_p_distances_parsed_output.txt' has been created at '", getwd(), "'"))

    # Return the dataframe
    return(ASAP_output)

  }

}
