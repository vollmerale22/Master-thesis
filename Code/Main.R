# Load packages
library(data.table)
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
library(keras)
library(tensorflow)
library(quantreg)
library(MSwM)
library(BMA)
library(zoo)
library(tseries)
library(doBy)
library(MacroRF)
library(tidyr)
library(ggthemes)
library(devtools)
library(sandwich)
library(lmtest)
library(xtable)
library(stargazer)
library(dplyr)
library(scales)
library(data.table)
library(stats)

rm(list = ls())
# Function to check stationarity using the Augmented Dickey-Fuller (ADF) test
is_stationary <- function(series, significance_level = 0.05) {
  series <- na.omit(series)  # Remove NA values for the test
  adf_test <- adf.test(series)
  return(adf_test$p.value <= significance_level)  # Stationary if p-value <= significance level
}

# Function to check for seasonality (using ACF)
has_seasonality <- function(series, seasonal_lag = 12, acf_threshold = 0.2) {
  series <- na.omit(series)  # Remove NA values for ACF calculation
  acf_values <- acf(series, plot = FALSE)  # Compute ACF
  significant_lags <- which(abs(acf_values$acf) > acf_threshold)  # Significant lags
  return(seasonal_lag %in% significant_lags)  # Check if seasonal lag is significant
}
# Function to standardise a column
standardise <- function(series) {
  return((series - mean(series, na.rm = TRUE)) / sd(series, na.rm = TRUE))
}

# Load the Excel sheet
data <- read_excel("Data/EU DATA.xlsx")
# Data transformation: handle seasonality, stationarity, and standardisation
transformed_data <- data %>%
  mutate(across(where(is.numeric), ~ {
    series <- .  # Current column data
    
    # Check for seasonality
    if (has_seasonality(series)) {
      print(paste(cur_column(), "has seasonality. Applying seasonal differencing."))
      series <- c(rep(NA, 12), diff(series, lag = 12))  # Apply seasonal differencing
    }
    
    # Check for stationarity after handling seasonality
    if (!is_stationary(series)) {
      print(paste(cur_column(), "is non-stationary after seasonal adjustment. Applying differencing."))
      series <- c(NA, diff(series))  # Apply regular differencing
    }
    
    # Standardise the series
    print(paste(cur_column(), "is now stationary. Standardising."))
    return(standardise(series))
  }))

# Save the transformed data to a new Excel file
output_file <- "Data/transformed_data_EU.xlsx"  # Replace with your desired output file name
write_xlsx(transformed_data, output_file)

# ----------------------------------------------------------------------------
# STEP 0 - SET USER PARAMETERS FOR THE HORSERACE
# ----------------------------------------------------------------------------
name_input <- "Data/transformed_data_EU.xlsx"  # Name of the input file
min_start_date <- "2015-11-15"          # Minimum start date for variables (otherwise discarded)
start_date_oos <- "2024-01-15"          # Start date for OOS predictions
end_date_oos <- "2024-07-15"            # End date for OOS predictions

list_h <- c(0)                      # List of horizons for back-, now- or fore-cast takes place
# Negative for a back-cast, 0 for a now-cast and positive for a fore-cast
list_methods <- c(2)                  # List of pre-selection methods
# 0 = No pre-selection
# 1 = LARS (Efron et al., 2004)
# 2 = Correlation-based (SIS: Fan and Lv, 2008)
# 3 = t-stat based (Bair et al., 2006)
# 4 = Iterated Bayesian Model Averaging (BMA: Yeung et al., 2005)
list_n <- c(40,60)                      # List of number of variables kept after pre-selection
list_reg <- c(1,2,3,4,5,6,7)                      # List of regressions techniques
# The AR benchmark is always performed - regardless of selection
# 1 = OLS
# 2 = Markov-switching regression [requires 1]
# 3 = Quantile regression
# 4 = Random forest
# 5 = XG Boost tree
# 6 = Macroeconomic Random Forest
# 7 = XG Boost linear

# User parameters - 3
# Optional inputs (calibration)
do_factors <- 1                         # Switch on whether to do factors or not
# 0 = no factors
# 1 = factors
fast_bma <- 1                           # 1 = fast version - i.e. runs with less iterations
# 0 = full number of iterations
n_per <- 12                             # Number of periods (last available ones) on which the optimization of the hyper-parameters for ML is performed
sc_ml <- 1                              # Switch on whether the data should be scaled for ML methods (0 = no, 1 = yes)
fast_MRF <- 1                           # 1 = fast tuning only on number of variables in linear part
# 0 = full tuning on 5 hyper-parameters (see code)


# Load functions
source("Code/Clean.R")         # Program cleaning the dataset (interpolation of missing data, elimination of variables starting after min_start_date)
source("Code/Align.R")         # Program doing the re-alignment - see section 2.1 of Chinn et al. (2023)
source("Code/PreSelection.R")  # Program doing the pre-selection - see section 1.3 of Chinn et al. (2023)
source("Code/Regression.R")    # Program doing the regressions - see section 1.5 of Chinn et al. (2023)
source("COde/Tuning.R")        # Program doing the tuning of hyper-parameters for machine learning methods


# ----------------------------------------------------------------------------
# STEP 2 - READ (AND CLEAN) THE USER-PROVIDED DATASET
# ----------------------------------------------------------------------------

# Read data
data_init <- read_excel(paste0(name_input)) %>%
  mutate(date=ymd(date))

# Clean data-set
data_rmv <- doClean(data_init,min_start_date)


# ----------------------------------------------------------------------------
# STEP 3 - PERFORM THE HORSERACE (LOOP OVER USER PARAMETERS)
# ----------------------------------------------------------------------------

# Loop over user-defined horizons
for (hh in 1:length(list_h)){
  
  # Get horizon to test
  horizon = as.numeric(list_h[hh])
  
  # Get re-aligned dataset
  data_real <- doAlign(data_rmv,horizon)
  
  # Initializing the output
  summary_ps_meth <- data.frame(NA)
  
  # Loop over user-defined methods
  for (mm in 1:length(list_methods)){
    
    # Initializing the output of the results comparison exercise
    select_method = as.numeric(list_methods[mm])
    
    # Loop over user-defined number of variables
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
        rm(var_sel)
        
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
                         scale = TRUE,
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
          #rm(eig_val)
          rm(res_pca)
          rm(rhs_fct)
          rm(rhs_fct_sel)
          rm(n_fct)
          
        }else{
          
          # Create alternative dataset without factors and CSV output
          don_cb <- cbind(lhs_sel,select(x_pca,-date))
          
        }
        
        # Cleaning
        rm(rhs_sel)
        rm(end_date)
        rm(start_date)
        rm(lhs_sel)
        rm(x_pca)
        
        
        # --------------------------------------------------------------------
        # STEP C: Regression on factors
        # --------------------------------------------------------------------
        
        # Get current date
        don_reg <- don_cb %>%
          drop_na()
        n_date_ii <- which(grepl(date_ii, as.character(don_reg$date)))
        
        # Define datasets in- and out-of-sample
        smpl_in <- head(don_reg,n_date_ii-horizon-3)
        smpl_out <- dplyr::slice(don_reg,n_date_ii)
        results[ii-n_start+1,1] <- date_ii
        
        # Run regressions
        temp_res <- run_regressions(smpl_in,
                                    smpl_out,
                                    list_reg,
                                    n_sel,
                                    sc_ml,
                                    fast_MRF)
        
        # Write results  
        results[ii-n_start+1,2:length(results)] <- temp_res
        colnames(results) <- c('date',colnames(temp_res))
        
      }
      
      # Create mean of predictions (excluding AR)
      results %<>%
        mutate(pred_mean = rowMeans(select(results,starts_with("pred_"))),
               date = as_date(as.Date(as.POSIXct(date*24*60*60, origin="1970-01-01"))),
               year = year(date))
      
      # Create directory
      dir.create(file.path(""), showWarnings = FALSE)
      dir.create(file.path(paste0("",horizon)), showWarnings = FALSE)
      
      # Write results (predictions)
      write.csv(select(results,-year), file = {output_dir <- paste0("./2-Output/",
                                                                    paste0("h", horizon, "/"))
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }
        paste0("./2-Output/",
                                                     paste0("h",horizon,"/"),
                                                     "pred_sel_",
                                                     select_method,
                                                     "_n_",
                                                     n_var,
                                                     "_reg_",
                                                     #paste(list_reg, collapse='_'),
                                                     "_h_",
                                                     horizon,
                                                     "_",
                                                     start_date_oos,
                                                     "_",
                                                     end_date_oos,
                                                     ".csv")},
                row.names=FALSE)
      
      # Summarizing and writing aggregate results
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
      write.csv(summary_all, file = paste0("./2-Output/",
                                           paste0("h",horizon,"/"),
                                           "rmse_sel_",
                                           select_method,
                                           "_n_",
                                           n_var,
                                           "_reg_",
                                           paste(list_reg, collapse='_'),
                                           "_h_",
                                           horizon,
                                           "_",
                                           start_date_oos,
                                           "_",
                                           end_date_oos,
                                           ".csv"))
      
      # Write results to the summary
      tot_meth <- length(list_methods)
      num_col <- 1+(mm-1)*tot_meth+nn
      summary_ps_meth[1,num_col] <- select_method
      summary_ps_meth[2,num_col] <- n_var
      summary_ps_meth[3:(nrow(summary_all)+2),num_col] <- summary_all$total
      summary_ps_meth[3:(nrow(summary_all)+2),1] <- rownames(summary_all)
      
    } # End of loop on number of variables
  } # End of loop on pre-selection method
  
  # Write results
  # This is a summary of all the tested methods
  summary_ps_meth[1,1] <- "pre-selection"
  summary_ps_meth[2,1] <- "nb variables"
  write.csv(summary_ps_meth, file = paste0("./2-Output/",
                                           paste0("h",horizon,"/"),
                                           "summaryALL_sel_",
                                           paste(list_methods,collapse='_'),
                                           "_n_",
                                           paste(list_n,collapse='_'),
                                           "_reg_",
                                           paste(list_reg, collapse='_'),
                                           "_h_",
                                           horizon,
                                           "_",
                                           start_date_oos,
                                           "_",
                                           end_date_oos,
                                           ".csv"),
            row.names=FALSE)
  
} # End of loop on horizon
