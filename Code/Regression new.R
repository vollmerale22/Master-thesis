run_regressions_new <- function(data_all,lhs_sel, x_pca, don_cb, smpl_in, smpl_out, list_methods, sc_ml, fast_MRF) {
  #n_sel
  
  # Prepare results dataset: true value and predictions for each method (plus AR)
  results <- data.frame(matrix(NA,
                               nrow = 1,
                               ncol = length(list_methods) + 2))
  results[1, 1] <- smpl_out[1, 2]
  names_col <- c("true_value", "ar")
  
  # Number of factors (used in MS regression)
  n_fct <- length(smpl_in) - 4
  
  # Prepare training and test samples for ML methods
  y_train <- as.matrix(select(smpl_in, target))[, 1]
  y_train_sc <- scale(y_train)[, 1]
  x_train <- as.matrix(select(smpl_in, -date, -target, -L1st_target, -L2nd_target))
  x_train_sc <- scale(x_train)
  x_test <- as.matrix(select(smpl_out, -date, -target, -L1st_target, -L2nd_target))
  
  # For scaling consistency, combine all and then extract test row
  x_all <- rbind(x_train, x_test)
  x_all_sc <- scale(x_all)
  x_test_sc <- t(as.matrix(x_all_sc[nrow(x_all_sc), ]))
  
  if (sc_ml == 1) {
    y_train_ml <- y_train_sc
    x_train_ml <- x_train_sc
    x_test_ml <- x_test_sc
    x_all_ml <- x_all_sc
  } else {
    y_train_ml <- y_train
    x_train_ml <- x_train
    x_test_ml <- x_test
    x_all_ml <- x_all
  }
  
  # AR model 
  max_lag <- 2
  
  # Initialize vectors to store models and their AIC values
  models <- vector("list", max_lag)
  bic_values <- numeric(max_lag)
  
  # Loop over possible lag lengths
  for (lag in 1:max_lag) {
    if (lag == 1) {
      data_model <- select(smpl_in, target, L1st_target)
      new_data <- as.data.frame(select(smpl_out, L1st_target))
    } else if (lag == 2) {
      data_model <- select(smpl_in, target, L1st_target, L2nd_target)
      new_data <- as.data.frame(select(smpl_out, L1st_target, L2nd_target))
    }
    
    # Fit the AR model with the specified number of lags
    model_temp <- lm(target ~ ., data = data_model, na.action = "na.omit")
    models[[lag]] <- model_temp
    bic_values[lag] <- BIC(model_temp)
  }
  
  # Select the model with the lowest AIC
  best_lag <- which.min(bic_values)
  best_model <- models[[best_lag]]
  
  # Use the selected model for prediction
  results[1, 2] <- predict(best_model,
                           newdata = new_data,
                           na.action = "na.omit")
  
  
  count_col <- 3  # starting column for ML methods
  
  if (-2 %in% list_methods) {
    print("Forecasting with DFM whole data")
    dfm_data1 <- data_all %>% select(-date) %>% drop_na()
    
    ic_out1 <- ICr(dfm_data1, max.r = 5) # Bai-Ng criterion for factors
    r_selected1 <- which.min(ic_out1$IC[,1])
    
    var_sel1 <- vars::VARselect(dfm_data1, lag.max = 5, type = "const") # Select the lag order for the VAR model
    p_optimal1 <- as.numeric(var_sel1$selection["SC(n)"])
    # Estimate the final DFM with the chosen r and p.
    dfm_model1 <- DFM(X = dfm_data1, r = r_selected1, p = p_optimal1, em.method = "none")
    
    # Forecast the target variable h steps ahead.
    # (Assume that predict() from dfms returns a forecast for the common factors and/or directly for the target.)
    dfm_pred1 <- predict(dfm_model1, h = 1)
    
    results[1, count_col] <- dfm_pred1$X_fcst[1, "target"]
    count_col <- count_col + 1
    names_col <- c(names_col, "pred_dfm_alldata")
  }
  
  if (-1 %in% list_methods) {
    print("Forecasting with DFM after preselection")
    dfm_data <- cbind(lhs_sel, 
                      select(x_pca, setdiff(colnames(x_pca), colnames(lhs_sel)))) %>% 
      select(-date) %>% 
      drop_na()
    
    
    ic_out <- ICr(dfm_data, max.r = 5) # Bai-Ng criterion for factors
    r_selected <- which.min(ic_out$IC[,1])
    
    var_sel <- vars::VARselect(dfm_data, lag.max = 5, type = "const") # Select the lag order for the VAR model
    p_optimal <- as.numeric(var_sel$selection["SC(n)"])
    # Estimate the final DFM with the chosen r and p.
    dfm_model <- DFM(X = dfm_data, r = r_selected, p = p_optimal, em.method = "none")
    
    # Forecast the target variable h steps ahead.
    # (Assume that predict() from dfms returns a forecast for the common factors and/or directly for the target.)
    dfm_pred <- predict(dfm_model, h = 1)
    
    results[1, count_col] <- dfm_pred$X_fcst[1, "target"]
    count_col <- count_col + 1
    names_col <- c(names_col, "pred_dfm_after_pre")
  }
  
  #DFM
  if (0 %in% list_methods) {
    print("Forecasting with DFM double pca")
    dfm_data3 <- don_reg %>%drop_na() %>% select(-date)
  
    ic_out3 <- ICr(dfm_data3, max.r = 5) # Bai-Ng criterion for factors
    r_selected3 <- which.min(ic_out3$IC[,1])

    var_sel3 <- vars::VARselect(dfm_data3, lag.max = 5, type = "const") # Select the lag order for the VAR model
    p_optimal3 <- as.numeric(var_sel3$selection["SC(n)"])
    # Estimate the final DFM with the chosen r and p.
    dfm_model3 <- DFM(X = dfm_data3, r = r_selected3, p = p_optimal3, em.method = "none")
    
    # Forecast the target variable h steps ahead.
    # (Assume that predict() from dfms returns a forecast for the common factors and/or directly for the target.)
    dfm_pred3 <- predict(dfm_model3, h = 1)
    
    results[1, count_col] <- dfm_pred3$X_fcst[1, "target"]
    count_col <- count_col + 1
    names_col <- c(names_col, "pred_dfm_doublepca")
  }
  
  
  
  # OLS
  if (1 %in% list_methods) {
    print("Forecasting with OLS")
    eq_lm_allf <- lm(target ~ .,
                     data = select(smpl_in, target, 5:ncol(smpl_in)),
                     na.action = "na.omit")
    results[1, count_col] <- predict(eq_lm_allf,
                                     newdata = select(smpl_out, 5:ncol(smpl_in)),
                                     na.action = "na.omit")
    count_col <- count_col + 1
    names_col <- c(names_col, "pred_ols")
  }
  
  # Markov-switching
  if (2 %in% list_methods) {
    print("Forecasting with Markov-switching")
    set.seed(1234)
    out <- tryCatch({
      eq_ms_allf <- msmFit(eq_lm_allf,
                           k = 2,
                           sw = rep(TRUE, n_fct + 2))
      temp1 <- as.matrix(eq_ms_allf@Coef[1, ])
      temp2 <- as.matrix(eq_ms_allf@Coef[2, ])
      temp_var <- as.matrix(select(smpl_out, 5:ncol(smpl_in)))
      temp_var <- cbind(1, temp_var)
      pred_1 <- sum(temp1 * temp_var)
      pred_2 <- sum(temp2 * temp_var)
      preds <- c(pred_1, pred_2)
      proba <- eq_ms_allf@Fit@smoProb[nrow(eq_ms_allf@Fit@smoProb), ]
      proba_trans <- eq_ms_allf@transMat
      if (proba[1] > proba[2]) {
        temp_prob <- proba_trans[, 1]
        pred <- sum(temp_prob * preds)
      } else {
        temp_prob <- proba_trans[, 2]
        pred <- sum(temp_prob * preds)
      }
    }, error = function(e) {
      return(NA)
    })
    
    if (is.na(out)) {
      pred <- predict(eq_lm_allf,
                      newdata = select(smpl_out, 5:ncol(smpl_in)),
                      na.action = "na.omit")
    }
    
    results[1, count_col] <- pred
    count_col <- count_col + 1
    names_col <- c(names_col, "pred_ms")
  }
  
  # Quantile regression
  if (3 %in% list_methods) {
    print("Forecasting with Quantile Regression")
    eq_qr_allf <- rq(target ~ .,
                     data = select(smpl_in, target, 5:ncol(smpl_in)),
                     na.action = "na.omit")
    results[1, count_col] <- predict(eq_qr_allf,
                                     newdata = select(smpl_out, 5:ncol(smpl_in)),
                                     na.action = "na.omit")
    count_col <- count_col + 1
    names_col <- c(names_col, "pred_qr")
  }
  
  # Random Forest
  if (4 %in% list_methods) {
    set.seed(1234)
    # Here we pass n_per (global) as the initial window, use a fixed horizon (e.g., 12), and n_folds = 5.
    param <- tune_RF_new(x_train_ml, y_train_ml, 
                     initial_window = n_per, 
                     horizon = 12, 
                     n_folds = 5, 
                     seed = 1234)
    print("Forecasting with Random Forest")
    eq_rf <- randomForest(y = y_train_ml,
                          x = x_train_ml,
                          na.action = "na.omit",
                          do.trace = FALSE,
                          ntree = param$ntree,
                          mtry = param$mtry,
                          nodesize = param$nodesize,
                          maxnodes = param$maxnodes,
                          corr.bias = TRUE)
    results[1, count_col] <- predict(eq_rf,
                                     newdata = x_test_ml,
                                     na.action = "na.omit")
    if (sc_ml == 1) {
      results[1, count_col] <- results[1, count_col] * sd(y_train) + mean(y_train)
    }
    count_col <- count_col + 1
    names_col <- c(names_col, "pred_rf")
  }
  
  # XGBoost tree
  if (5 %in% list_methods) {
    set.seed(1234)
    param <- tune_XGBT_new(x_train_ml, y_train_ml, 
                       initial_window = n_per, 
                       horizon = 12, 
                       n_folds = 5, 
                       seed = 1234)
    temp_train <- xgb.DMatrix(data = x_train_ml, label = y_train_ml)
    print("Forecasting with gradient boosting (trees)")
    eq_boost <- xgboost(data = temp_train,
                        nrounds = param$nrounds,
                        eta = param$eta,
                        max_depth = param$max_depth,
                        min_child_weight = param$min_child_weight,
                        gamma = param$gamma,
                        objective = "reg:squarederror",
                        verbose = 0)
    results[1, count_col] <- predict(eq_boost,
                                     newdata = x_test_ml,
                                     na.action = "na.omit")
    if (sc_ml == 1) {
      results[1, count_col] <- results[1, count_col] * sd(y_train) + mean(y_train)
    }
    count_col <- count_col + 1
    names_col <- c(names_col, "pred_xgbt")
  }
  
  # Macroeconomic Random forest
  if(6 %in% list_methods){
    
    # Prepare data
    data.in <- rbind(cbind(y_train_ml,x_train_ml),cbind(NA,x_test_ml))
    
    # Set seed for reproducibility
    # NB: has to be the same as in DoTuning.R
    set.seed(1234)
    
    if(fast_MRF==0){
      
      # Calling a tuning function from DoTuning.R
      # The tuning is done on the last n_per periods
      param <- tune_MRF(data.in,
                        n_per)
      
      print("Forecasting with macroeconomic random forest")
      eq_mrf <- MRF(data.in,
                    B = param[1],
                    x.pos = c(2:param[2]),
                    oos.pos = nrow(data.in),
                    cheap.look.at.GTVPs = FALSE,
                    ridge.lambda = param[3],
                    resampling.opt = param[4],
                    block.size = param[5],
                    printb = FALSE)
      
    }else if(fast_MRF==1){
      
      # Calling a tuning function from DoTuning.R
      # The tuning is done on the last n_per periods
      param <- tune_MRF_fast(data.in,
                             n_per)
      
      print("Forecasting with macroeconomic random forest")
      eq_mrf <- MRF(data.in,
                    x.pos = c(2:param),
                    oos.pos = nrow(data.in),
                    cheap.look.at.GTVPs = FALSE,
                    printb = FALSE)
      
    }
    
    results[1,count_col] <- eq_mrf$pred
    if(sc_ml==1){results[1,count_col] <- results[1,count_col]*sd(y_train) + mean(y_train)}
    
    count_col <- count_col + 1
    names_col <- c(names_col,"pred_mrf") 
    
  }
  
  
  
  # XGBoost linear
  if (7 %in% list_methods) {
    param <- tune_XGBL_new(x_train_ml, y_train_ml, 
                       initial_window = n_per, 
                       horizon = 12, 
                       n_folds = 5, 
                       seed = 1234)
    print("Forecasting with gradient linear boosting")
    temp_train <- xgb.DMatrix(data = x_train_ml, label = y_train_ml)
    eq_boost <- xgboost(data = temp_train,
                        nrounds = as.numeric(param$nrounds),
                        eta = as.numeric(param$eta),
                        alpha = as.numeric(param$alpha),
                        booster = "gblinear",
                        objective = "reg:squarederror",
                        verbose = 0)
    results[1, count_col] <- predict(eq_boost,
                                     newdata = x_test_ml,
                                     na.action = "na.omit")
    if (sc_ml == 1) {
      results[1, count_col] <- results[1, count_col] * sd(y_train) + mean(y_train)
    }
    count_col <- count_col + 1
    names_col <- c(names_col, "pred_xgbl")
  }
  
  
  # LSTM Nowcasting
  if (8 %in% list_methods) {
    set.seed(1234)
    lag <- 2
    
    # Tune LSTM with a 2-lag structure.
    param <- tune_LSTM(x_train_ml, y_train_ml, 
                       initial_window = n_per, 
                       horizon = 12, 
                       n_folds = 5, 
                       seed = 1234, 
                       lag = lag)
    print("Forecasting with LSTM")
    
    # Create training sequences for the final model.
    num_train <- nrow(x_train_ml) - lag + 1
    x_train_seq <- array(NA, dim = c(num_train, lag, ncol(x_train_ml)))
    for (j in 1:num_train) {
      x_train_seq[j, , ] <- as.matrix(x_train_ml[j:(j+lag-1), ])
    }
    y_train_seq <- y_train_ml[lag:length(y_train_ml)]
    
    # For the test sample, we need a sequence of length 'lag'.
    # Assuming x_test_ml is a single new observation, combine it with the last (lag - 1) rows of x_train_ml.
    x_test_seq <- array(NA, dim = c(1, lag, ncol(x_test_ml)))
    if (lag > 1) {
      x_test_seq[1, 1:(lag-1), ] <- as.matrix(x_train_ml[(nrow(x_train_ml) - lag + 2):nrow(x_train_ml), ])
    }
    x_test_seq[1, lag, ] <- as.matrix(x_test_ml)
    
    # Build the final LSTM model.
    model <- keras_model_sequential() %>%
      layer_lstm(
        units = param$units,
        input_shape = c(lag, ncol(x_train_ml)),
        dropout = param$dropout,
        recurrent_dropout = param$recurrent_dropout,
        return_sequences = TRUE
      ) %>%
      layer_lstm(
        units = param$units,
        dropout = param$dropout,
        recurrent_dropout = param$recurrent_dropout
      ) %>%
      layer_dense(units = 1)
    
    model %>% compile(
      loss = "mean_squared_error",
      optimizer = optimizer_adam(lr = 0.001)
    )
    
    early_stop <- callback_early_stopping(monitor = "val_loss", patience = 5)
    lr_reduce <- callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.5, patience = 3)
    
    history <- model %>% fit(
      x = x_train_seq,
      y = y_train_seq,
      epochs = param$epochs,
      batch_size = param$batch_size,
      verbose = 1,
      validation_split = 0.2,
      callbacks = list(early_stop, lr_reduce)
    )
    
    pred <- model %>% predict(x_test_seq)
    pred <- as.numeric(pred)
    
    if (sc_ml == 1) {
      pred <- pred * sd(y_train_ml) + mean(y_train_ml)
    }
    
    results[1, count_col] <- pred
    count_col <- count_col + 1
    names_col <- c(names_col, "pred_lstm_custom")
    
    rm(model, history, x_train_seq, x_test_seq, pred)
    keras::k_clear_session()
    gc()
  }
  
  colnames(results) <- names_col
  return(results)
}
