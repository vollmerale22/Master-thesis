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
  old_plan <- future::plan()         # Save the current plan
  on.exit(future::plan(old_plan), add = TRUE)  # Restore the plan on exit
  
  # Use availableCores() and cap the number of workers at 20 to avoid oversubscription
  n_workers <- min(future::availableCores(), 20)
  future::plan(multisession, workers = n_workers)
  
  # Apply the adjustment function to each numeric column (excluding "date")
  adjusted_list <- future.apply::future_lapply(names(data)[-1], function(col_name) {
    adjust_column(data[[col_name]], col_name)
  }, future.seed = TRUE)
  
  # Reassemble the adjusted data
  adjusted_data <- as.data.frame(adjusted_list)
  colnames(adjusted_data) <- names(data)[-1]  # Restore original column names
  adjusted_data <- cbind(date = data$date, adjusted_data)  # Reattach the date column  
  return(adjusted_data)
}
