penalized_weighted_median <- function(values, weights, ratios) {
  # Combine values, weights, and ratios into a single dataframe
  data <- data.frame(values, weights, ratios)
  
  # Remove rows with NA in any column
  data <- na.omit(data)
  
  # Create a new weighting scheme that penalizes the values
  # Here we multiply each value by its square and by the number of reads
  data$penalized_weights <- data$ratios^2 * data$weights
  
  # Order by values
  data <- data[order(data$values), ]
  
  # Calculate the cumulative penalized weights
  data$cum_weights <- cumsum(data$penalized_weights)
  
  # Total sum of penalized weights
  total_weights <- sum(data$penalized_weights)
  
  # Find the median based on the cumulative penalized weights
  median_index <- min(which(data$cum_weights >= total_weights / 2))
  return(data$values[median_index])
}
