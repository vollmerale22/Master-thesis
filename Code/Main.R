# Load packages
# Load packages
library(reticulate)
#reticulate::install_python()
virtualenv_install("r-tf", packages = c("tensorflow==2.15.1", "keras==2.15.0"))
use_virtualenv("r-tf", required = TRUE)
library(data.table)
library(tseries)
library(forecast)
library(purrr)
library(tidyverse)
library(readr)
library(lars)
library(readxl)
library(writexl)
library(haven)
library(magrittr)
library(haven)
library(lubridate)
library(FactoMineR)
library(factoextra)
library(randomForest)
library(glmnet)
library(xgboost)
library(quantreg)
library(MSwM)
library(BMA)
library(zoo)
library(tseries)
library(doBy)
library(devtools)
#install_github("philgoucou/macrorf", force = TRUE)
library(MacroRF)
library(tidyr)
library(ggthemes)
library(sandwich)
library(lmtest)
library(xtable)
library(stargazer)
library(dplyr)
library(scales)
library(data.table)
library(stats)
library(seasonal)
library(lubridate)
library(future.apply)
library(keras)
install_keras()
library(tensorflow)
library(doParallel)
library(foreach)
library(parallel)
library(dfms)


num_cores <- detectCores() - 2  # Reserve a couple of cores for the system
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat("Parallel backend set up with", num_cores, "cores.\n")


rm(list = ls())
# Function to check if a time series is stationary
is_stationary <- function(series) {
  clean_series <- na.omit(series)  # Remove NAs
  
  if (length(clean_series) < 5 || length(unique(clean_series)) == 1) {
    return(TRUE)  # Assume stationary if not enough data or constant
  }
  
  test_result <- tryCatch({
    suppressWarnings(adf.test(clean_series, k = 0))  # Remove ADF test warnings
  }, error = function(e) {
    warning("Stationarity test failed:", e$message)
    return(list(p.value = 1))  # Assume non-stationary if test fails
  })
  
  return(test_result$p.value < 0.05)  # Stationary if p-value < 0.05
}

# Function to standardise a column
standardise <- function(series) {
  return((series - mean(series, na.rm = TRUE)) / sd(series, na.rm = TRUE))
}

# Load the Excel sheet
data <- read_excel("Data/Data combined China Exports.xlsx")

# Load functions
source("Code/Clean.R")         
source("Code/Align.R")        
source("Code/PreSelection.R")  
#source("Code/Regression.R")    
#source("Code/Tuning.R")
source("Code/Seasonality.R")
source("Code/Tuning new.R")
source("Code/Regression new.R")

# ----------------------------------------------------------------------------
# Adjust seasonality and standardise data 
transformed_data <- data %>%
      adjust_seasonality_large()  # Apply seasonal adjustment only to these
# Replace columns 2-11 with the values from data
transformed_data[2:11] <- data[2:11]


final_transformed_data <- transformed_data %>%
  mutate(across(where(is.numeric), ~ {
    series <- .
    # Check for stationarity
    if (!is_stationary(series)) {
      series <- c(NA, diff(series))  # Apply regular differencing
    }
    return(standardise(series))
  }))

# Save the transformed data to a new Excel file
output_file <- final_transformed_data
output_file <- "Data/transformed_data_testEU.xlsx"  
write_xlsx(final_transformed_data, output_file)

# ----------------------------------------------------------------------------
# STEP 0
# ----------------------------------------------------------------------------
target_variable <- "EU imports"
name_input <- output_file  # Name of the input file
min_start_date <- "2001-01-15"          # Minimum start date for variables (otherwise discarded)
start_date_oos <- "2012-01-15"          # Start date for OOS predictions
end_date_oos <- "2024-07-15"            # End date for OOS predictions
list_h <- c(0)                      # List of horizons for back-, now- or fore-cast takes place
# Negative for a back-cast, 0 for a now-cast and positive for a fore-cast
list_methods <- c(1)                  # List of pre-selection methods
# 0 = No pre-selection
# 1 = LARS (Efron et al., 2004)
# 2 = Correlation-based (SIS: Fan and Lv, 2008)
# 3 = t-stat based (Bair et al., 2006)
# 4 = Iterated Bayesian Model Averaging (BMA: Yeung et al., 2005)
list_n <- c(60)                      # List of number of variables kept after pre-selection
list_reg <- c(-2,-1,0,1,2,3,4,5,6,7,8)  # List of regressions techniques
#-2 DFM on whole data
#-1 = DFM on pre-selected variables
# 0 = DFM double PCA
# 1 = OLS
# 2 = Markov-switching regression [requires 1]
# 3 = Quantile regression
# 4 = Random forest
# 5 = XG Boost tree
# 6 = Macroeconomic Random Forest
# 7 = XG Boost linear
# 8 = LSTM
# Optional inputs (calibration)
do_factors <- 1                         # Switch on whether to do factors or not
# 0 = no factors
# 1 = factors
fast_bma <- 0                          # 1 = fast version - i.e. runs with less iterations
# 0 = full number of iterations
n_per <- 12                             # Number of periods (last available ones) on which the optimization of the hyper-parameters for ML is performed
sc_ml <- 1                              # Switch on whether the data should be scaled for ML methods (0 = no, 1 = yes)
fast_MRF <- 0                         # 1 = fast tuning only on number of variables in linear part
# 0 = full tuning on 5 hyper-parameters




# ----------------------------------------------------------------------------
# STEP 2 - READ (AND CLEAN)
# ----------------------------------------------------------------------------

# Read data
data_init <- read_excel(paste0(name_input)) %>%
  mutate(date=ymd(date))

# Clean data-set
data_rmv <- doClean(data_init,min_start_date)


# ----------------------------------------------------------------------------
# STEP 3 - PERFORM THE LOOP
# ----------------------------------------------------------------------------

# Loop over defined horizons
for (hh in 1:length(list_h)){
  
  # Get horizon to test
  horizon = as.numeric(list_h[hh])
  
  # Get re-aligned dataset
  data_real <- doAlign(data_rmv,horizon)
  
  # Initializing the output
  summary_ps_meth <- data.frame(NA)
  
  # Loop over defined methods
  for (mm in 1:length(list_methods)){
    
    # Initialising the output of the results comparison exercise
    select_method = as.numeric(list_methods[mm])
    
    # Loop over defined number of variables
    for (nn in 1:length(list_n)){
      
      # Get values of the test
      n_var = as.numeric(list_n[nn])
      
      # Print
      print('======================================================================================================')
      print(paste0("HORIZON = ",horizon))
      print(paste0("METHOD = ",select_method))
      print(paste0("NUMBER OF VARIABLES = ",n_var))
      print('======================================================================================================')
      
      # Determine starting and ending dates for out-of-sample exercise
      n_start <- which(grepl(start_date_oos, as.character(data_real$date)))
      n_end <- which(grepl(end_date_oos, as.character(data_real$date)))
      
      # Prepare results dataset
      results <- data.frame(matrix(NA,
                                   nrow = n_end - n_start +1,
                                   ncol = length(list_reg)+3))
      
      # Loop over dates for out-of-sample predictions
      for (ii in n_start:n_end){
        
        # Data and dates
        data_all <- head(data_real,ii)
        date_ii <- data_all$date[ii]
        year <- substr(date_ii,1,4)
        month <- substr(date_ii,6,7)
        
        # Print
        print(paste0("Doing out-of-sample predictions for ",year," at month ",month))
        
        
        # --------------------------------------------------------------------
        # STEP A: Pre-selection
        # --------------------------------------------------------------------
        
        # Run pre-selection function
        var_sel <- pre_select(data_real,ii,horizon,select_method,n_var)
        
        # Check that the number of variables is correct
        if(length(var_sel)!=n_var && select_method != 2){
          stop("Pre-selection step did not work, please check")  
        }
        
        # Get LHS and RHS datasets with only the variables from the pre-selection
        rhs_sel <- data_all %>%
          select(date,all_of(var_sel))
        
        lhs_sel <- data_all %>%
          select(date,
                 target,
                 L1st_target,
                 L2nd_target)
        
        # Cleaning
        #rm(var_sel)
        
        # --------------------------------------------------------------------
        # STEP B: Factor extraction (PCA) on pre-selected variables
        # --------------------------------------------------------------------
        
        # Clean NAs in the RHS data
        x_pca <- rhs_sel %>%
          drop_na()
        
        # Clean the NAs also in LHS 
        start_date <- as.Date(x_pca$date[1])
        end_date <- as.Date(x_pca$date[nrow(x_pca)])
        
        lhs_sel %<>%
          filter(date>=start_date) %>%
          filter(date<=end_date)
        
        if(do_factors==1){
          
          # Run PCA
          x_pca %<>% select(-date)
          
          res_pca <- PCA(X = x_pca,
                         scale = FALSE,
                         graph = FALSE,
                         ncp = length(x_pca))
          
          rhs_fct <- res_pca$ind$coord %>%
            as.data.table()
          
          # Keep only the factors with eigen values > 1 = rule of the thumb of Kaiser (1961)
          # If eigen value > 1 then explain more variance than single variable
          n_fct <- dim(rhs_fct)[2]
          rhs_fct_sel <- rhs_fct %>%
            select(1:all_of(n_fct))
          
          # Create final dataset for regression
          don_cb <- cbind(lhs_sel, rhs_fct_sel)
          
          # Cleaning
          rm(res_pca)
          rm(rhs_fct)
          rm(rhs_fct_sel)
          rm(n_fct)
          
        }else{
          
          # Create alternative dataset without factors and CSV output
          don_cb <- cbind(lhs_sel,select(x_pca,-date))
          
        }
        
        # Cleaning
        #rm(rhs_sel)
        #rm(end_date)
        #rm(start_date)
        #rm(lhs_sel)
        #rm(x_pca)
        
        
        # --------------------------------------------------------------------
        # STEP C: Regression on factors
        # --------------------------------------------------------------------
        
        # Get current date
        don_reg <- don_cb %>%
          drop_na()
        n_date_ii <- which(grepl(date_ii, as.character(don_reg$date)))
        
        # Define datasets in- and out-of-sample
        smpl_in <- head(don_reg,n_date_ii-horizon-3)
        smpl_out <- as.data.frame(dplyr::slice(don_reg,n_date_ii), drop=FALSE)
        results[ii-n_start+1,1] <- date_ii
        # Run regressions
        temp_res <- run_regressions_new(data_all,
                                        don_cb,
                                        lhs_sel, 
                                        x_pca,
                                        smpl_in,
                                        smpl_out,
                                        list_reg,
                                        sc_ml,
                                        fast_MRF)
        
        # Write results  
        results[ii-n_start+1,2:length(results)] <- temp_res
        colnames(results) <- c('date',colnames(temp_res))
        
      }
      
      # Create mean of predictions (excluding AR)
      #results %<>%
       ## mutate(pred_mean = rowMeans(select(results,starts_with("pred_"))),
         #      date = as_date(as.Date(as.POSIXct(date*24*60*60, origin="1970-01-01"))),
          #     year = year(date))
      
      
      # Summarising and writing aggregate results
      err <- function(X,Y){sqrt(sum((X-Y)^2)/length(Y))}
      
      total <- results %>%
        select(-date,-year) %>%
        apply(2,err,Y=results[,2])
      
      crisis <- results %>%
        filter((year %in% c(2008,2009,2020,2021))) %>%
        select(-date,-year)
      
      crisis <- apply(crisis,2,err,Y=crisis[,1])
      
      normal <- results %>%
        filter(!(year %in% c(2008,2009,2020,2021))) %>%
        select(-date,-year)
      
      normal <- apply(normal,2,err,Y=normal[,1])
      
      summary_all <- cbind(cbind(total,crisis),normal) %>%
        as.data.frame()
      summary_all <- summary_all[!(row.names(summary_all) %in% c('true_value')),]
      
      # Writing results
      # Results are written for RMSE on crisis sample (2008-2009 and 2020-2021) and non-crisis sample
      output_rmse <- paste0("./Output/",
                            paste0("h", horizon, "/"))
      if (!dir.exists(output_rmse)) {
        dir.create(output_rmse, recursive = TRUE)
      }
      write.csv(summary_all, file = paste0("./Output/",
                                           paste0("h",horizon,"/"),
                                           "rmse_",
                                           target_variable,
                                           "_reg_",
                                           paste(list_reg, collapse='_'),
                                           "_h_",
                                           horizon,
                                           "_",
                                           start_date_oos,
                                           "_",
                                           end_date_oos,
                                           ".csv"))
    } # End of loop on number of variables
  } # End of loop on pre-selection method
} # End of loop on horizon
stopCluster(cl)

#Plot the first 12 rows  of data with date being on the x axis
library(ggplot2)
library(reshape2)
library(ggplot2)


results <- read.csv("./Final results/China/Exports/results_China_Exports.csv", row.names = 1)
summary_all <- read.csv("./Final results/China/Exports/summaryall_China_Exports.csv" ,row.names = 1)
resultsDFM <- read.csv("./Final results/China/Exports/LSTM_results_China_Exports.csv", row.names = 1)
summary_allDFM <- read.csv("./Final results/China/Exports/LSTM_summaryall_China_Exports.csv", row.names = 1)

#add results_LSTM_custom from resultsDFM to results by date
results <- results %>%
  left_join(resultsDFM %>% select(date, pred_lstm_custom), by = "date")
#add summary_allDFM to summary_all by model
summary_all <- bind_rows(summary_all, summary_allDFM["pred_lstm_custom",])

# ----------------------------------------------------------------------------
# STEP 4 - PERFORM TEST
# ----------------------------------------------------------------------------
source("Code/DM test.R")
benchmark_AR <- "ar"                  # AR benchmark column 
benchmark_DFM <- "pred_dfm_alldata"   # DFM benchmark column
# --- Run the DM tests for both benchmarks using the 'results' object:
dm_results_AR <- perform_dm_test(results, h = 1, loss = function(e) e^2, one_sided = TRUE, 
                                 harvey_correction = TRUE, benchmark_column = benchmark_AR)
dm_results_DFM <- perform_dm_test(results, h = 1, loss = function(e) e^2, one_sided = TRUE, 
                                  harvey_correction = TRUE, benchmark_column = benchmark_DFM)

# Save results to CSV
write.csv(dm_results_AR, file = paste0("./Output/",paste0("DM Tests AR2"), "dm_test_results_", target_variable, ".csv"))
write.csv(dm_results_DFM, file = paste0("./Output/",paste0("DM Tests DFM2"), "dm_test_results_", target_variable, ".csv"))


# ----------------------------------------------------------------------------
# STEP 5 - CREATE FINAL TABLE
# ----------------------------------------------------------------------------

display_names <- c("ar" = "AR",
                   "pred_dfm_alldata" = "DFM",
                   "pred_qr" = "QR",
                   "pred_ms" = "MS",
                   "pred_rf" = "RF",
                   "pred_mrf" = "MRF",
                   "pred_xgbt" = "XGBoost (Tree)",
                   "pred_lstm_custom" = "LSTM")

#--- Define the desired order and grouping ---
stat_models_keys <- c("ar", "pred_dfm_alldata", "pred_qr", "pred_ms")
ml_models_keys <- c("pred_rf", "pred_mrf", "pred_xgbt", "pred_lstm_custom")
ordered_keys <- c(stat_models_keys, ml_models_keys)

#--- Build the final summary table ---
final_table <- data.frame(Model = character(),
                          RMSE_Total = numeric(),
                          RMSE_Crisis = numeric(),
                          RMSE_Normal = numeric(),
                          Benchmark_AR = character(),
                          Benchmark_DFM = character(),
                          stringsAsFactors = FALSE)

# For benchmarks, extract RMSE from summary_all
rmse_AR_total <- summary_all["ar", "total"]
rmse_DFM_total <- summary_all["pred_dfm_alldata", "total"]

for (key in ordered_keys) {
  # Get display name for current model
  model_disp <- display_names[key]
  
  # Extract RMSE values from summary_all 
  rmse_model_total <- summary_all[key, "total"]
  rmse_model_crisis <- summary_all[key, "crisis"]
  rmse_model_normal <- summary_all[key, "normal"]
  
  # For benchmarks AR and DFM, no DM test is applied
  if (key == "ar") {
    bench_AR <- "--"
    # For AR row, you might display the relative improvement versus DFM
    bench_DFM <- paste0(round(((rmse_AR_total - rmse_DFM_total) / rmse_DFM_total) * 100,3), "%")
  } else if (key == "pred_dfm_alldata") {
    bench_AR <- paste0(round(((rmse_DFM_total - rmse_AR_total) / rmse_AR_total) * 100, 3), "%")
    bench_DFM <- "--"
  } else {
    # Relative improvement vs. AR benchmark:
    improvement_AR <- (rmse_model_total - rmse_AR_total) / rmse_AR_total * 100
    p_val_AR <- dm_results_AR[[key]]$p_value
    stars_AR <- ifelse(p_val_AR < 0.05, "**", ifelse(p_val_AR < 0.1, "*", ""))
    bench_AR <- paste0(round(improvement_AR, 1), "%", stars_AR)
    
    # Relative improvement vs. DFM benchmark:
    improvement_DFM <- (rmse_model_total - rmse_DFM_total) / rmse_DFM_total * 100
    p_val_DFM <- dm_results_DFM[[key]]$p_value
    stars_DFM <- ifelse(p_val_DFM < 0.05, "**", ifelse(p_val_DFM < 0.1, "*", ""))
    bench_DFM <- paste0(round(improvement_DFM, 1), "%", stars_DFM)
  }
  
  final_table <- rbind(final_table,
                       data.frame(Model = model_disp,
                                  RMSE_Total = rmse_model_total,
                                  RMSE_Crisis = rmse_model_crisis,
                                  RMSE_Normal = rmse_model_normal,
                                  Benchmark_AR = bench_AR,
                                  Benchmark_DFM = bench_DFM,
                                  stringsAsFactors = FALSE))
}

# --- Create xtable with dynamic caption and label ---
xt <- xtable(final_table, 
             caption = paste("Results Overview", target_variable),
             label = paste0("tab:", gsub(" ", "_", tolower(target_variable))))

# --- Insert custom header rows for grouping ---
# Create a list with 'pos' and 'command' elements
addtorow <- list(pos = list(), command = c())

# Before the first row, add a header for Statistical models
addtorow$pos[[1]] <- -1
addtorow$command[1] <- "\\hline\n\\multicolumn{6}{l}{\\textbf{Statistical models}} \\\\\n\\hline\n"

# Find the position where ML models start (after the last statistical model row)
stat_count <- length(stat_models_keys)
addtorow$pos[[2]] <- stat_count  # insert after the last statistical model row
addtorow$command[2] <- "\\hline\n\\multicolumn{6}{l}{\\textbf{ML Models}} \\\\\n\\hline\n"

# --- Print the table with xtable ---
print(xt, 
      digits = c(0, 0, 3, 3, 3, 0, 0),
      add.to.row = addtorow, 
      include.colnames = TRUE, 
      include.rownames = FALSE,
      floating = TRUE, 
      caption.placement = "top",
      sanitize.text.function = identity)


# ----------------------------------------------------------------------------
# STEP 6 - PERFORMANCE METRICS & DIRECTIONAL TESTS
# ----------------------------------------------------------------------------

library(rugarch)

# Calculate MAE and MAPE for each forecasting model
actual <- results[, 2]  
forecast_names <- names(results)[-(1:2)]  # Exclude 'date' and 'true_value' columns

mae <- sapply(forecast_names, function(model) {
  mean(abs(results[[model]] - actual), na.rm = TRUE)
})

mape <- sapply(forecast_names, function(model) {
  mean(abs((results[[model]] - actual) / actual), na.rm = TRUE) * 100
})

mae_mape_df <- data.frame(Model = forecast_names, MAE = mae, MAPE = mape)
print("MAE and MAPE for each model:")
print(mae_mape_df)

da_test_results <- list()
ep_test_results <- list()

# Loop through each forecasting model to perform the tests
for (model in forecast_names) {
  forecast_vals <- results[[model]]
  
  # Directional Accuracy Test (Pesaran and Timmerman)
  da_test_results[[model]] <- DACTest(forecast = forecast_vals, 
                                      actual = actual, 
                                      test = "PT", 
                                      conf.level = 0.95)
  
  # Excess Profitability Test (Anatolyev and Gerko)
  ep_test_results[[model]] <- DACTest(forecast = forecast_vals, 
                                      actual = actual, 
                                      test = "AG", 
                                      conf.level = 0.95)
}

# Print the test results
print("Directional Accuracy Test (Pesaran and Timmerman) Results:")
print(da_test_results)

print("Excess Profitability Test (Anatolyev and Gerko) Results:")
print(ep_test_results)


