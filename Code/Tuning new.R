# =======================================================================================================================================
# Random Forest

tune_RF_new <- function(x,y,n,initial_window, horizon, n_folds = 5, ntree = 300, seed = 1234) {

  print("Tuning random forest")
  
  set.seed(seed)
  
  # Define grid of hyperparameters
  grid <- expand.grid(
    mtry = seq(from = floor(ncol(x) / 3), to = ncol(x), length.out = 5),
    nodesize = seq(from = 3, to = 12, by = 3)
  )
  grid$RMSE <- NA
  
  n_total <- nrow(x)
  total_splits <- n_total - initial_window - horizon + 1
  if(total_splits < n_folds) {
    n_folds <- total_splits
  }
  splits <- sort(sample(1:total_splits, n_folds))
  
  for (i in 1:nrow(grid)) {
    cv_errors <- c()
    for (fold in splits) {
      train_end <- initial_window + fold - 1
      train_indices <- 1:train_end
      valid_indices <- (train_end + 1):(train_end + horizon)
      
      x_train <- x[train_indices, , drop = FALSE]
      y_train <- y[train_indices]
      x_valid <- x[valid_indices, , drop = FALSE]
      y_valid <- y[valid_indices]
      
      set.seed(seed + fold)
      model <- randomForest(x = x_train, y = y_train,
                            mtry = grid$mtry[i],
                            nodesize = grid$nodesize[i],
                            ntree = ntree)
      preds <- predict(model, newdata = x_valid)
      rmse <- sqrt(mean((preds - y_valid)^2))
      cv_errors <- c(cv_errors, rmse)
    }
    grid$RMSE[i] <- mean(cv_errors)
  }
  
  best <- grid[which.min(grid$RMSE), ]
  message("Best RF parameters using rolling origin CV:")
  print(best)
  return(best)
}


# =======================================================================================================================================
# XG Boost - tree

tune_XGBT_new <- function(x, y, initial_window, horizon, n_folds = 5, seed = 1234) {

  print("Tuning XG Boost tree")
  set.seed(seed)
  
  grid <- expand.grid(
    nrounds = seq(10, 100, by = 30),
    eta = c(0.05, 0.1, 0.2, 0.3),
    max_depth = c(3, 5, 6, 7, 9),
    min_child_weight = 1:4,
    gamma = c(0, 0.1, 0.3, 0.5)
  )
  grid$RMSE <- NA
  
  n_total <- nrow(x)
  total_splits <- n_total - initial_window - horizon + 1
  if(total_splits < n_folds) {
    n_folds <- total_splits
  }
  splits <- sort(sample(1:total_splits, n_folds))
  
  for (i in 1:nrow(grid)) {
    cv_errors <- c()
    for (fold in splits) {
      train_end <- initial_window + fold - 1
      train_indices <- 1:train_end
      valid_indices <- (train_end + 1):(train_end + horizon)
      
      x_train <- x[train_indices, , drop = FALSE]
      y_train <- y[train_indices]
      x_valid <- x[valid_indices, , drop = FALSE]
      y_valid <- y[valid_indices]
      
      dtrain <- xgb.DMatrix(data = x_train, label = y_train)
      
      set.seed(seed + fold)
      model <- xgboost(
        data = dtrain,
        nrounds = grid$nrounds[i],
        eta = grid$eta[i],
        max_depth = grid$max_depth[i],
        min_child_weight = grid$min_child_weight[i],
        gamma = grid$gamma[i],
        objective = "reg:squarederror",
        verbose = 0
      )
      preds <- predict(model, newdata = x_valid)
      rmse <- sqrt(mean((preds - y_valid)^2))
      cv_errors <- c(cv_errors, rmse)
    }
    grid$RMSE[i] <- mean(cv_errors)
  }
  
  best <- grid[which.min(grid$RMSE), ]
  message("Best XGBoost Tree parameters using rolling origin CV:")
  print(best)
  return(best)
}

# =======================================================================================================================================
# Macroeconomic Random Forest

# Full version
tune_MRF_new <- function(data.in, initial_window, horizon = 12, n_folds = 5, seed = 1234) {
  # data.in: Matrix (or data.frame) with target in column 1 and predictors in remaining columns
  # initial_window: Number of initial observations used for training in the first fold
  # horizon: Number of periods to forecast (e.g., 12)
  # n_folds: Number of rolling CV folds
  # seed: For reproducibility
  
  set.seed(seed)
  print("Tuning macroeconomic random forest")
  # Determine candidate values for n_var (number of variables in the linear part)
  total_preds <- ncol(data.in) - 1  # predictors count
  m_min <- min(5, total_preds)
  m_max <- max(total_preds, 9)
  # Ensure we don’t exceed available predictors (x.pos are columns 2 to n_var)
  n_var_candidates <- seq(from = m_min, to = min(m_max, total_preds + 1), by = 2)
  
  # Create grid for all hyper-parameters
  grid <- expand.grid(
    n_trees = c(25, 50),
    n_var   = n_var_candidates,
    lambda  = c(0.001, 0.01, 0.1),
    re_meth = c(2, 4),
    bl_size = c(6, 12)
  )
  grid$RMSE <- NA
  
  n <- nrow(data.in)
  total_splits <- n - initial_window - horizon + 1
  if (total_splits < n_folds) {
    n_folds <- total_splits
  }
  # Select n_folds evenly spaced splits (or randomly sample splits)
  splits <- sort(sample(1:total_splits, n_folds))
  
  # Loop over each candidate parameter combination
  for (i in 1:nrow(grid)) {
    params <- grid[i, ]
    cv_errors <- c()
    
    for (fold in splits) {
      train_end <- initial_window + fold - 1
      data_tr <- data.in[1:train_end, , drop = FALSE]
      n_train <- nrow(data_tr)
      test_indices <- (n_train - horizon + 1):n_train
      y_te_cv <- data_tr[test_indices, 1]  # Actual target values for test period
      
      # Remove target values in test set to mimic out‑of‑sample forecast
      data_tr[test_indices, 1] <- NA
      
      # Set seed for this fold (using an offset)
      set.seed(seed + fold)
      model <- MRF(data_tr,
                   B = params$n_trees,
                   x.pos = 2:params$n_var,
                   oos.pos = nrow(data_tr),
                   cheap.look.at.GTVPs = FALSE,
                   printb = FALSE,
                   ridge.lambda = params$lambda,
                   resampling.opt = params$re_meth,
                   block.size = params$bl_size)
      pred <- model$pred
      error <- sqrt(mean((pred - y_te_cv)^2))
      cv_errors <- c(cv_errors, error)
    }
    grid$RMSE[i] <- mean(cv_errors)
  }
  
  best <- grid[which.min(grid$RMSE), ]
  message("Best full MRF parameters using rolling origin CV:")
  print(best)
  return(best)
}
# Faster version
tune_MRF_fast_new <- function(data.in, initial_window, horizon = 12, n_folds = 5, seed = 1234) {
  # data.in: Matrix (or data.frame) with target in column 1 and predictors in remaining columns
  # initial_window: Number of initial observations for training in the first fold
  # horizon: Forecast horizon (e.g., 12)
  # n_folds: Number of CV folds
  # seed: For reproducibility
  
  set.seed(seed)
  print("Tuning macroeconomic random forest (fast version)")
  
  total_preds <- ncol(data.in) - 1
  m_min <- min(5, total_preds)
  m_max <- max(total_preds, 9)
  n_var_candidates <- seq(from = m_min, to = min(m_max, total_preds + 1), by = 2)
  
  results <- data.frame(n_var = n_var_candidates, RMSE = NA)
  
  n <- nrow(data.in)
  total_splits <- n - initial_window - horizon + 1
  if (total_splits < n_folds) {
    n_folds <- total_splits
  }
  splits <- sort(sample(1:total_splits, n_folds))
  
  for (i in 1:length(n_var_candidates)) {
    cv_errors <- c()
    candidate <- n_var_candidates[i]
    
    for (fold in splits) {
      train_end <- initial_window + fold - 1
      data_tr <- data.in[1:train_end, , drop = FALSE]
      n_train <- nrow(data_tr)
      test_indices <- (n_train - horizon + 1):n_train
      y_te_cv <- data_tr[test_indices, 1]
      
      # Remove test targets to simulate forecasting
      data_tr[test_indices, 1] <- NA
      
      set.seed(seed + fold)
      model <- MRF(data_tr,
                   x.pos = 2:candidate,
                   oos.pos = nrow(data_tr),
                   cheap.look.at.GTVPs = FALSE,
                   printb = FALSE)
      pred <- model$pred
      error <- sqrt(mean((pred - y_te_cv)^2))
      cv_errors <- c(cv_errors, error)
    }
    results$RMSE[i] <- mean(cv_errors)
  }
  
  best <- results[which.min(results$RMSE), ]
  message("Best fast MRF parameter (n_var) using rolling origin CV:")
  print(best)
  return(best)
}


# =======================================================================================================================================
# XG Boost - linear

tune_XGBL_new <- function(x, y, initial_window, horizon, n_folds = 5, seed = 1234) {
  set.seed(seed)

  print("Tuning linear gradient boosting")
  grid <- expand.grid(
    nrounds = seq(10, 130, by = 30),
    eta = c(0.05, 0.1, 0.2, 0.3),
    alpha = c(0.001, 0.01, 0.1, 0.2, 0.3)
  )
  grid$RMSE <- NA
  
  n_total <- nrow(x)
  total_splits <- n_total - initial_window - horizon + 1
  if(total_splits < n_folds) {
    n_folds <- total_splits
  }
  splits <- sort(sample(1:total_splits, n_folds))
  
  for (i in 1:nrow(grid)) {
    cv_errors <- c()
    for (fold in splits) {
      train_end <- initial_window + fold - 1
      train_indices <- 1:train_end
      valid_indices <- (train_end + 1):(train_end + horizon)
      
      x_train <- x[train_indices, , drop = FALSE]
      y_train <- y[train_indices]
      x_valid <- x[valid_indices, , drop = FALSE]
      y_valid <- y[valid_indices]
      
      dtrain <- xgb.DMatrix(data = x_train, label = y_train)
      
      set.seed(seed + fold)
      model <- xgb.train(
        data = dtrain,
        nrounds = grid$nrounds[i],
        eta = grid$eta[i],
        alpha = grid$alpha[i],
        booster = "gblinear",
        objective = "reg:squarederror",
        verbose = 0
      )
      preds <- predict(model, newdata = x_valid)
      rmse <- sqrt(mean((preds - y_valid)^2))
      cv_errors <- c(cv_errors, rmse)
    }
    grid$RMSE[i] <- mean(cv_errors)
  }
  
  best <- grid[which.min(grid$RMSE), ]
  message("Best XGBoost Linear parameters using rolling origin CV:")
  print(best)
  return(best)
}
tune_LSTM <- function(x, y, initial_window, horizon, n_folds = 5, seed = 1234) {
  # x: predictor matrix (ordered in time)
  # y: target vector (ordered in time)
  # initial_window: number of initial observations for training in the first fold
  # horizon: forecast horizon for each fold
  # n_folds: number of rolling CV folds
  # seed: reproducibility seed
  
  set.seed(seed)
  
  # Define grid for LSTM hyperparameters.
  # Here we tune: units (number of LSTM cells), dropout rate, epochs, and batch_size.
  grid <- expand.grid(
    units = c(10, 20, 50),
    dropout = c(0, 0.2),
    epochs = c(50, 100),
    batch_size = c(32, 64)
  )
  grid$RMSE <- NA
  
  n_total <- nrow(x)
  total_splits <- n_total - initial_window - horizon + 1
  if(total_splits < n_folds) {
    n_folds <- total_splits
  }
  splits <- sort(sample(1:total_splits, n_folds))
  
  for(i in 1:nrow(grid)){
    cv_errors <- c()
    for(fold in splits){
      train_end <- initial_window + fold - 1
      train_indices <- 1:train_end
      valid_indices <- (train_end + 1):(train_end + horizon)
      
      x_train <- x[train_indices, , drop = FALSE]
      y_train <- y[train_indices]
      x_valid <- x[valid_indices, , drop = FALSE]
      y_valid <- y[valid_indices]
      
      # Reshape predictors to 3D array: [samples, timesteps, features]
      x_train_array <- array(x_train, dim = c(nrow(x_train), 1, ncol(x_train)))
      x_valid_array <- array(x_valid, dim = c(nrow(x_valid), 1, ncol(x_valid)))
      
      # Build a simple LSTM model
      model <- keras_model_sequential() %>%
        layer_lstm(units = grid$units[i],
                   input_shape = c(1, ncol(x_train)),
                   dropout = grid$dropout[i]) %>%
        layer_dense(units = 1)
      
      model %>% compile(
        loss = "mean_squared_error",
        optimizer = "adam"
      )
      
      set.seed(seed + fold)
      history <- model %>% fit(
        x = x_train_array,
        y = y_train,
        epochs = grid$epochs[i],
        batch_size = grid$batch_size[i],
        verbose = 0
      )
      
      preds <- model %>% predict(x_valid_array)
      rmse <- sqrt(mean((preds - y_valid)^2))
      cv_errors <- c(cv_errors, rmse)
      
      # Clear session to free up memory
      k_clear_session()
    }
    grid$RMSE[i] <- mean(cv_errors)
  }
  
  best <- grid[which.min(grid$RMSE), ]
  message("Best LSTM parameters using rolling origin CV:")
  print(best)
  return(best)
}