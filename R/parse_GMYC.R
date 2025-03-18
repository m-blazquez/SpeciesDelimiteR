
parse_GMYC <- function(GMYC_output, method = NULL) {

  # Ensure a valid method is specified
  if (is.null(method) || !method %in% c("single", "multiple")) {

    stop("ERROR: You must specify a valid method ('single' or 'multiple').")

  }

  # Parse the GMYC output based on the method used to produce it
  if(method == "single") {

    # Create new dataframe
    parsed_file <- data.frame(GMYC_output$sample_name, GMYC_output$GMYC_spec)

    # Change column names
    colnames(parsed_file) <- c("Specimen", "GMYCs_species")

    # Ensure the directory exists before writing the file
    if (!dir.exists("GMYC/Single")) {
      dir.create("GMYC/Single", recursive = TRUE)
    }

    # Write new file
    write.table(parsed_file,
                file = "GMYC/Single/GMYCs_parsed_output.txt",
                sep = "\t",
                row.names = FALSE)

    # Tell the user where the new file is
    message(paste0("The file 'GMYCs_parsed_output.txt' has been created at '", getwd(), "/GMYC/Single'"))

    # Return the dataframe
    return(parsed_file)

  }

  if(method == "multiple") {

    # Create new dataframe
    parsed_file <- data.frame(GMYC_output$sample_name, GMYC_output$GMYC_spec)

    # Change column names
    colnames(parsed_file) <- c("Specimen", "GMYCm_species")

    # Ensure the directory exists before writing the file
    if (!dir.exists("GMYC/Multiple")) {
      dir.create("GMYC/Multiple", recursive = TRUE)
    }

    # Write new file
    write.table(parsed_file,
                file = "GMYC/Multiple/GMYCm_parsed_output.txt",
                sep = "\t",
                row.names = FALSE)

    # Tell the user where the new file is
    message(paste0("The file 'GMYCm_parsed_output.txt' has been created at '", getwd(), "/GMYC/Multiple'"))

    # Return the dataframe
    return(parsed_file)

  }

}
