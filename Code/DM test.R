perform_dm_test <- function(results, h = 1, loss = function(e) e^2) {
  # 'results' must have a numeric 'true_value' column, a numeric 'ar' column, 
  # and other numeric forecast columns (e.g. 'pred_rf', 'pred_mean', etc.).
  #
  # h: forecast horizon (default = 1)
  # loss: a function that takes forecast errors and returns the loss (default = squared error)
  
  # 1) Identify the actual (true) values
  if (!("true_value" %in% names(results))) {
    stop("No 'true_value' column found in 'results'!")
  }
  actual <- results[["true_value"]]
  if (!is.numeric(actual)) {
    stop("'true_value' must be numeric!")
  }
  
  # 2) Identify the benchmark forecast (e.g. 'ar')
  if (!("ar" %in% names(results))) {
    stop("No benchmark forecast column 'ar' found in 'results'!")
  }
  benchmark <- results[["ar"]]
  
  # 3) Internal DM test function comparing two forecasts
  dm_test <- function(actual, f1, f2, h = 1, loss = function(e) e^2) {
    e1 <- actual - f1
    e2 <- actual - f2
    d <- loss(e1) - loss(e2)
    dbar <- mean(d, na.rm = TRUE)
    T_val <- length(na.omit(d))
    var_d <- var(d, na.rm = TRUE)
    # Simple DM statistic (no Newey-West for h>1)
    DM_stat <- dbar / sqrt(var_d / T_val)
    p_value <- 2 * (1 - pnorm(abs(DM_stat)))
    return(list(DM_statistic = DM_stat, p_value = p_value))
  }
  
  # 4) Identify all forecast columns (exclude known columns)
  forecast_names <- setdiff(names(results), c("date", "year", "true_value", "ar"))
  
  # 5) Compare each forecast to the benchmark
  dm_results <- list()
  for (fn in forecast_names) {
    dm_results[[fn]] <- dm_test(actual, results[[fn]], benchmark, h, loss)
  }
  
  return(dm_results)
}


