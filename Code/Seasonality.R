adjust_seasonality_large <- function(data, freq = 12, method = "SEATS", handle_NAs = "omit") {
  # Convert date column to Date format and sort by date
  data$date <- as.Date(data$date)
  data <- data %>% arrange(date)
  
  # Extract the starting year and month from the first date
  start_year <- year(data$date[1])
  start_month <- month(data$date[1])
  start <- c(start_year, start_month)
  
  # Function to adjust a single column
  adjust_column <- function(column, column_name) {
    ts_data <- ts(column, frequency = freq, start = start)
    
    # Try seasonal adjustment
    seas_result <- try(seas(ts_data, na.action = na.exclude, transform.function = "none", outlier = NULL), silent = TRUE)
    
    if (inherits(seas_result, "try-error")) {
      warning(paste("Seasonal adjustment failed for column:", column_name, "- Returning original column"))
      return(column)  # Return original column if adjustment fails
    }
    adjusted_series <- final(seas_result)
    
    # Ensure adjusted series is valid (not all NA, not shorter)
    if (length(adjusted_series) != length(column) || all(is.na(adjusted_series))) {
      warning(paste("Adjusted series is invalid for column:", column_name, "- Returning original column"))
      return(column)
    }
    
    return(adjusted_series)  # Successfully adjusted series
  }
  
  # Apply function to each numeric column (excluding "date") using parallel processing
  plan(multisession, workers = parallel::detectCores() - 1)
  
  adjusted_data <- as.data.frame(future_lapply(names(data)[-1], function(col_name) {
    adjust_column(data[[col_name]], col_name)
  }, future.seed = TRUE))
  
  colnames(adjusted_data) <- names(data)[-1]  # Keep original column names
  adjusted_data <- cbind(date = data$date, adjusted_data)  # Reattach date
  
  return(adjusted_data)
}

