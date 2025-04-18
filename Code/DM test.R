perform_dm_test <- function(results, 
                            h = 1, 
                            loss = function(e) e^2, 
                            one_sided = TRUE, 
                            harvey_correction = TRUE,
                            benchmark_column) {
  # 1) Identify the actual (true) values
  actual <- results[["true_value"]]
  benchmark <- results[[benchmark_column]]
  
  # 3) Internal DM test function comparing two forecasts
  dm_test <- function(actual, f1, f2, h = 1, loss = function(e) e^2, 
                      one_sided, harvey_correction) {
    e1 <- actual - f1
    e2 <- actual - f2
    d <- loss(e1) - loss(e2)
    dbar <- mean(d, na.rm = TRUE)
    T_val <- sum(!is.na(d))
    var_d <- var(d, na.rm = TRUE)
    
    # Simple DM statistic (no Newey-West for h>1)
    DM_stat <- dbar / sqrt(var_d / T_val)
    
    # Apply Harvey et al. (1997) correction if requested
    if (harvey_correction) {
      CF <- sqrt((T_val + 1 - 2 * h + (h * (T_val - 1)) / T_val) / T_val)
      DM_stat <- DM_stat * CF
    }
    
    # Compute two-sided p-value
    p_value <- 2 * (1 - pnorm(abs(DM_stat)))
    
    # Adjust for one-sided test if required
    if (one_sided) {
      if (dbar < 0) {
        p_value <- pnorm(DM_stat)      # left-tail
      } else {
        p_value <- 1 - pnorm(DM_stat)    # right-tail
      }
    }
    
    return(list(DM_statistic = DM_stat, p_value = p_value, mean_loss_diff = dbar))
  }
  
  # 4) Identify all forecast columns (exclude known columns)
  forecast_names <- setdiff(names(results), c("date", "year", "true_value", benchmark_column))
  
  # 5) Compare each forecast to the benchmark
  dm_results <- list()
  for (fn in forecast_names) {
    dm_results[[fn]] <- dm_test(
      actual,
      f1 = results[[fn]],
      f2 = benchmark,
      h = h,
      loss = loss,
      one_sided = one_sided,
      harvey_correction = harvey_correction
    )
  }
  
  return(dm_results)
}


