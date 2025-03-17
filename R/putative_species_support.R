
putative_species_support <- function(integrated_results, consensus_species) {

  # Check if the inputs are in the correct format
  if (!is.data.frame(integrated_results) | !is.data.frame(consensus_species)) {

    stop("Both inputs must be dataframes.")

  }

  # Ensure 'Specimen' column is integer or character in both dataframes
  if (!("Specimen" %in% colnames(integrated_results)) | !("Specimen" %in% colnames(consensus_species))) {

    stop("Both dataframes must contain a 'Specimen' column.")

  }

  # Check that consensus_species has two columns
  if (ncol(consensus_species) != 2) {

    stop("consensus_species must have exactly two columns: Specimen and Consensus.")

  }

  # Check that all columns in integrated_results (except Specimen) are factors
  if (any(!sapply(integrated_results[, -1], is.factor))) {

    stop("All columns except 'Specimen' in integrated_results must be factors.")

  }

  # Convert Specimen columns to character for comparison
  integrated_results$Specimen <- as.character(integrated_results$Specimen)
  consensus_species$Specimen <- as.character(consensus_species$Specimen)

  # Make sure Specimen columns match between both dataframes
  if (!all(consensus_species$Specimen %in% integrated_results$Specimen)) {

    stop("Specimen columns do not match between the two dataframes.")

  }

  # Merge dataframes on Specimen
  merged_data <- merge(consensus_species, integrated_results, by = "Specimen")

  # Create an empty result dataframe to store the support information
  result_df <- data.frame(Consensus = levels(consensus_species$Consensus))
  strategies <- colnames(integrated_results)[-1]  # Exclude the Specimen column
  result_df <- cbind(result_df, matrix(NA, nrow = nrow(result_df), ncol = length(strategies)))
  colnames(result_df)[2:(length(strategies) + 1)] <- strategies

  # Evaluate support for each consensus species across strategies
  for (i in 1:nrow(result_df)) {

    consensus_species_value <- result_df$Consensus[i]
    species_support <- rep(NA, length(strategies))

    for (j in 1:length(strategies)) {

      strategy_column <- strategies[j]

      # Get specimens of the current consensus species
      specimens_in_consensus <- merged_data$Specimen[merged_data$Consensus == consensus_species_value]

      # Check if all specimens for the consensus species are ascribed to the same species in the strategy column
      species_in_strategy <- merged_data[[strategy_column]][merged_data$Specimen %in% specimens_in_consensus]
      unique_species_in_strategy <- unique(species_in_strategy)

      if (length(unique_species_in_strategy) == 1) {

        species_support[j] <- "Yes"

      } else {

        species_support[j] <- "No"

      }

    }

    # Assign supported species for each strategy to the result_df
    result_df[i, 2:(length(strategies) + 1)] <- species_support

  }

  # Count the number of strategies supporting each consensus species
  result_df$Supported_strategies <- rowSums(result_df[, 2:(length(strategies) + 1)] == "Yes", na.rm = TRUE)

  # Calculate the percentage of support for each consensus species
  result_df$Support_percentage <- (result_df$Supported_strategies / length(strategies)) * 100

  # Write new file
  write.table(result_df,
              file = "Putative_species_support.txt",
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)

  # Tell the user where the new file is
  message(paste0("The file 'Putative_species_support.txt' has been created at '", getwd(), "'"))

  # Return the result dataframe
  return(result_df)
}

