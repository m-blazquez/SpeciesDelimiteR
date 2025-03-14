
count_putative_species <- function(integrated_results) {

  # Check if input is a dataframe
  if (!is.data.frame(integrated_results)) {

    stop("Error: The input must be a data frame.")

  }

  # Check if the first column is named "Specimen"
  if (names(integrated_results)[1] != "Specimen") {

    stop("Error: The first column must be named 'Specimen'.")

  }

  # Exclude the "Specimen" column
  species_columns <- setdiff(names(integrated_results), "Specimen")

  # Check if all additional columns are factors with at least 1 level
  for (col in species_columns) {

    if (!is.factor(integrated_results[[col]]) || length(levels(integrated_results[[col]])) < 1) {

      stop(paste("Error: Column", col, "must be a factor with at least one level."))

    }

  }

  # Initialize an empty vector to store counts
  species_counts <- c()

  # Loop through each species discovery strategy and count unique species
  for (col in species_columns) {

    species_counts <- c(species_counts, length(unique(na.omit(integrated_results[[col]]))))

  }

  # Create a data frame
  species_counts_df <- data.frame(Strategy = species_columns,
                                  Putative_Species_Count = species_counts,
                                  stringsAsFactors = FALSE)

  # Write new file
  write.table(species_counts_df,
              file = "Putative_species_counts.txt",
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)

  # Tell the user where the new file is
  message(paste0("The file 'Putative_species_counts.txt' has been created at '", getwd(), "'"))

  # Return the dataframe
  return(species_counts_df)
}
