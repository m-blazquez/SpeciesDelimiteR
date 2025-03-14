
integrate_results <- function(...) {

  # Capture input dataframes
  dfs <- list(...)

  # Make sure some dataframe was provided
  if (length(dfs) == 0) {

    stop("Error: Nothing to integrate: No dataframes were provided")

  }

  # Make sure more than one dataframe was provided
  if (length(dfs) == 1) {

    stop("Error: Nothing to integrate: Just one dataframe was provided")

  }

  # Check if all dataframes have exactly two columns and the first column is "Specimen"
  for (i in seq_along(dfs)) {

    if (ncol(dfs[[i]]) != 2) {

      stop(paste("Error: Dataframe", i, "does not have exactly two columns"))

    }

    if (colnames(dfs[[i]])[1] != "Specimen") {

      stop(paste("Error: The first column of dataframe", i, "is not named 'Specimen'"))

    }
  }

  # Normalize Specimen column: convert to character and trim whitespace
  for (i in seq_along(dfs)) {

    dfs[[i]]$Specimen <- trimws(as.character(dfs[[i]]$Specimen))

  }

  # Extract specimen columns and check if they contain the same values (ignoring order)
  specimen_sets <- lapply(dfs, function(df) sort(df$Specimen))

  if (!all(sapply(specimen_sets, function(x) identical(x, specimen_sets[[1]])))) {

    stop("Error: Specimen columns do not match across all dataframes (check for data type mismatches, extra spaces, or missing values)")

  }

  # Merge all dataframes
  merged_df <- dfs[[1]]  # Use the first dataframe as the base

  for (i in 2:length(dfs)) {

    merged_df <- merge(merged_df, dfs[[i]], by = "Specimen")

  }

  # Convert all columns except "Specimen" to factors
  for (col in names(merged_df)[-1]) {

    merged_df[[col]] <- as.factor(merged_df[[col]])

  }

  # Write new file
  write.table(merged_df,
              file = "Integrated_results.txt",
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)

  # Tell the user where the new file is
  message(paste0("The file 'Integrated_results.txt' has been created at '", getwd(), "'"))

  # Return the dataframe
  return(merged_df)
}
