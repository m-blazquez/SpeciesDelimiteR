
parse_bPTP <- function(bPTP_output_filename) {

  # Read the file
  lines <- readLines(bPTP_output_filename, warn = FALSE)  # Suppress warning

  # Ensure the last line has a newline (if missing)
  if (length(lines) > 0 && lines[length(lines)] != "") {
    lines <- c(lines, "")  # Append an empty line
  }

  # Remove comment lines
  lines <- lines[!grepl("^#", lines)]

  # Initialize variables
  species_list <- list()
  current_species <- NA

  # Loop through each line
  for (line in lines) {
    line <- trimws(line)  # Remove leading/trailing whitespace

    if (grepl("^Species", line)) {
      # Extract species number
      current_species <- as.numeric(sub("Species (\\d+) .*", "\\1", line))
    } else if (grepl("^\\S", line)) {  # Ensure the line has content (not empty)
      # Extract specimen IDs as character
      specimens <- unlist(strsplit(line, ","))
      specimens <- trimws(specimens)  # Remove extra spaces

      # Append to list
      species_list <- append(species_list, list(data.frame(Specimen = specimens, bPTP_species = current_species, stringsAsFactors = FALSE)))
    }
  }

  # Combine all data frames
  bPTP <- do.call(rbind, species_list)

  # Write new file
  write.table(bPTP,
              file = "bPTP_parsed_output.txt",
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)  # Avoid unnecessary quotes around specimen names

  # Tell the user where the new file is
  message(paste0("The file 'bPTP_parsed_output.txt' has been created at '", getwd(), "'"))

  # Return the dataframe
  return(bPTP)
}
