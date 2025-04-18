doClean <- function(data_init,min_start_date){
  
  # Inputs:
  # data_init [data.frame] = initial data-set read from Stata
  # min_start_date [string] = minimum start date for variables
  #                           NB: should be in "YYYY-MM-15" format (pay attention to put "15" for the day)
  # Output: 
  # data_rmv = cleaned data-set 
  
  # -----------------------------------
  # 1 - Remove variables starting after the date set by the user
  # -----------------------------------
  
  data_rmv <- data_init
  
  # Remove variables starting after min_start_date
  ctrl_min <- data_rmv %>%
    filter(date==as.Date(min_start_date)) %>%
    reshape2::melt(id.vars = "date", 
                   variable.name = "variable",
                   value.name = "value") %>%
    filter(!is.na(value)) %>%
    mutate(variable=as.character(variable))
  
  var_to_keep <- c(unique(ctrl_min$variable),"date")
  
  data_rmv %<>%
    dplyr::select(all_of(var_to_keep))
  
  
  # -----------------------------------
  # 2 - Clean data (interpolate missing observations and remove non-stationary variables)
  # -----------------------------------
  
  # Interpolate data missing between two points (but not leading or trailing NAs)
  data_rmv %<>%
    mutate(across(-date, ~na.approx(.,na.rm = FALSE)))
  
  # Remove non-stationary variables - by ADF test, removing if p-value > 10%
  temp <- data_rmv %>%
    dplyr::select(-date,-target)
  var_nonstat <- c()
  
  for (jj in 1:ncol(temp)){
    
    # Select individual variable and drop NAs
    temp_col <- temp[jj]
    temp_col <- drop_na(as.data.frame(temp_col))
    temp_name <- colnames(temp_col)
    
    # Run ADF test and add variable to removal list if p-value > 10%
    temp_adf <- adf.test(temp_col[,1])
    if(temp_adf$p.value > 0.10){
      var_nonstat <- c(var_nonstat,temp_name)
    }
  }
  
  data_rmv %<>%
    dplyr::select(-all_of(var_nonstat))
  
  
  # -----------------------------------
  # 3 - Return data-set
  # -----------------------------------
  
  return(data_rmv)
  
}


x <- arima.sim(list(order = c(1,1,0),ar = 0.2),n = 100)
adf.test(x)
