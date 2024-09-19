library(FourgameteP)
library(vcfR)


vcf <- read.vcfR("isotypes.vcf", verbose = FALSE )
gt <- t(extract.gt(vcf, element = c('GT'), as.numeric = TRUE))
gt <- gt+1
gt[is.na(gt)] <- 0

# Define the window size
window_size <- 24


# Loop through each window
for (i in seq(1, ncol(gt), by = window_size)) {
  
  # Define the subset of columns for the current window
  cols <- i:min(i + window_size - 1, ncol(gt))
  
  # Subset the data
  subset_data <- as.data.frame(gt[, cols])
  
  # Apply FGf function to the subset_data
  FGf(subset_data)
  
  # Store the result in the list
  file_name <- paste("result_", cols[1], "_to_", cols[length(cols)], sep = "")
  
  assign(file_name, FinalResults)
  
}

# List all data frames in the Global Environment that start with "result_"
df_names <- ls(pattern = "^result_")

# Initialize an empty list to store successfully combined data frames
combined_dfs <- list()

# Initialize an empty list to store names of data frames causing errors
error_dfs <- list()

# Loop through each data frame name
for (df_name in df_names) {
  # Try to retrieve the data frame
  df <- get(df_name)
  
  # Check if it's a data frame
  if (is.data.frame(df)) {
    # Ensure all data frames have the same column names as the first one
    if (length(combined_dfs) == 0) {
      combined_dfs[[df_name]] <- df
    } else {
      # Check if column names match with the first data frame
      if (all(names(df) %in% names(combined_dfs[[1]])) && all(names(combined_dfs[[1]]) %in% names(df))) {
        combined_dfs[[df_name]] <- df
      } else {
        # Store the name of the data frame causing the error
        cat("Error binding:", df_name, "- Column names do not match\n")
        error_dfs[[df_name]] <- df_name
      }
    }
  } else {
    # If not a data frame, skip and print a message
    cat("Skipping non-data frame object:", df_name, "\n")
  }
}

# Combine all data frames that were successfully bound into one using rbind
if (length(combined_dfs) > 0) {
  combined_data <- do.call(rbind, combined_dfs)
} else {
  cat("No data frames were successfully combined.\n")
}
