# =======================================================================================================================================
# Random Forest

tune_RF_new <- function(x, y, n, initial_window, horizon, n_folds = 5, seed = 1234) {
  cat("Tuning Random Forest\n")
  set.seed(seed)
  
  grid <- expand.grid(
    ntree = c(100, 200, 300, 500, 1000),
    mtry = seq(from = floor(ncol(x) / 3), to = ncol(x), length.out = 5),
    nodesize = seq(from = 3, to = 12, by = 3),
    maxnodes = c(5, 10, 15, 20)
  )
  grid$RMSE <- NA_real_
  
  n_total <- nrow(x)
  total_splits <- n_total - initial_window - horizon + 1
  if(total_splits < n_folds) n_folds <- total_splits
  splits <- floor(seq(1, total_splits, length.out = n_folds))
  
  res <- foreach(i = 1:nrow(grid), .combine = rbind, .packages = "randomForest") %dopar% {
    cv_errors <- numeric(n_folds)
    for (fold in seq_along(splits)) {
      train_end <- initial_window + splits[fold] - 1
      train_indices <- 1:train_end
      valid_indices <- (train_end + 1):(train_end + horizon)
      
      x_train <- x[train_indices, , drop = FALSE]
      y_train <- y[train_indices]
      x_valid <- x[valid_indices, , drop = FALSE]
      y_valid <- y[valid_indices]
      
      set.seed(seed + splits[fold])
      model <- randomForest(x = x_train, y = y_train,
                            mtry = grid$mtry[i],
                            nodesize = grid$nodesize[i],
                            ntree = grid$ntree[i],
                            maxnodes = grid$maxnodes[i])
      preds <- predict(model, newdata = x_valid)
      cv_errors[fold] <- sqrt(mean((preds - y_valid)^2))
    }
    data.frame(index = i, RMSE = mean(cv_errors))
  }
  grid$RMSE <- res$RMSE
  best <- grid[which.min(grid$RMSE), ]
  message("Best RF parameters using parallel rolling origin CV:")
  print(best)
  return(best)
}


# =======================================================================================================================================
# XG Boost - tree

tune_XGBT_new <- function(x, y, initial_window, horizon, n_folds = 5, seed = 1234, n_iter = 70) {
  print("Tuning XGBoost tree with rolling CV")
  set.seed(seed)
  
  # Reduced grid example – adjust as needed
  nrounds_options <- seq(10, 100, by = 30)
  eta_options <- c(0.05, 0.1, 0.2, 0.3)
  max_depth_options <- c(3, 5, 6, 7, 9)
  min_child_weight_options <- 1:4
  gamma_options <- c(0, 0.1, 0.3, 0.5)
  
  full_grid <- expand.grid(nrounds = nrounds_options,
                           eta = eta_options,
                           max_depth = max_depth_options,
                           min_child_weight = min_child_weight_options,
                           gamma = gamma_options)
  sampled_grid <- full_grid[sample(nrow(full_grid), n_iter), ]
  sampled_grid$RMSE <- NA_real_
  
  n_total <- nrow(x)
  total_splits <- n_total - initial_window - horizon + 1
  if(total_splits < n_folds) n_folds <- total_splits
  splits <- floor(seq(1, total_splits, length.out = n_folds))
  
  res <- foreach(i = 1:nrow(sampled_grid), .combine = rbind, .packages = "xgboost") %dopar% {
    cv_errors <- numeric(n_folds)
    for (fold in seq_along(splits)) {
      train_end <- initial_window + splits[fold] - 1
      train_indices <- 1:train_end
      valid_indices <- (train_end + 1):(train_end + horizon)
      
      if(max(valid_indices) > n_total) next
      
      x_train <- x[train_indices, , drop = FALSE]
      y_train <- y[train_indices]
      x_valid <- x[valid_indices, , drop = FALSE]
      y_valid <- y[valid_indices]
      
      dtrain <- xgb.DMatrix(data = x_train, label = y_train)
      set.seed(seed + splits[fold])
      model <- xgboost(
        data = dtrain,
        nrounds = sampled_grid$nrounds[i],
        eta = sampled_grid$eta[i],
        max_depth = sampled_grid$max_depth[i],
        min_child_weight = sampled_grid$min_child_weight[i],
        gamma = sampled_grid$gamma[i],
        objective = "reg:squarederror",
        verbose = 0,
        nthread = 1
      )
      preds <- predict(model, newdata = x_valid)
      cv_errors[fold] <- sqrt(mean((preds - y_valid)^2))
    }
    data.frame(index = i, RMSE = mean(cv_errors, na.rm = TRUE))
  }
  
  sampled_grid$RMSE <- res$RMSE
  best <- sampled_grid[which.min(sampled_grid$RMSE), ]
  message("Best XGBoost Tree parameters using parallel random search and rolling origin CV:")
  print(best)
  return(best)
}


# =======================================================================================================================================
# XG Boost - linear
tune_XGBL_new <- function(x, y, initial_window, horizon, n_folds = 5, seed = 1234, n_iter = 50) {
  cat("Tuning linear gradient boosting (gblinear) with rolling window CV\n")
  set.seed(seed)
  
  full_grid <- expand.grid(
    nrounds = c(10,40,70,100),
    eta = c(0.1, 0.2, 0.3),
    alpha = c(0.001, 0.01, 0.1)
  )
  n_iter_actual <- min(n_iter, nrow(full_grid))
  sampled_grid <- full_grid[sample(nrow(full_grid), n_iter_actual), ]
  sampled_grid$RMSE <- NA_real_
  
  n_total <- nrow(x)
  total_splits <- n_total - initial_window - horizon + 1
  if(total_splits < n_folds) n_folds <- total_splits
  splits <- floor(seq(1, total_splits, length.out = n_folds))
  
  # Create a reusable function to compute splits inside the foreach loop
  # This avoids serializing large objects across nodes
  compute_split <- function(fold_idx, n_total, initial_window, horizon, splits) {
    if (fold_idx > length(splits)) return(NULL)
    
    train_end <- initial_window + splits[fold_idx] - 1
    train_indices <- 1:train_end
    valid_indices <- (train_end + 1):(train_end + horizon)
    
    if(max(valid_indices) > n_total) return(NULL)
    
    return(list(
      train_indices = train_indices,
      valid_indices = valid_indices
    ))
  }
  
  res <- foreach(i = 1:nrow(sampled_grid), .combine = rbind, .packages = "xgboost") %dopar% {
    cv_errors <- numeric(length(splits))
    
    for (fold in 1:length(splits)) {
      # Compute split indices inside the parallel loop
      split_info <- compute_split(fold, n_total, initial_window, horizon, splits)
      if(is.null(split_info)) next
      
      # Extract training and validation data for this fold
      x_train <- x[split_info$train_indices, , drop = FALSE]
      y_train <- y[split_info$train_indices]
      x_valid <- x[split_info$valid_indices, , drop = FALSE]
      y_valid <- y[split_info$valid_indices]
      
      # Create DMatrix on demand
      dtrain <- xgb.DMatrix(data = x_train, label = y_train)
      
      set.seed(seed + splits[fold])
      model <- xgb.train(
        params = list(
          booster = "gblinear",
          objective = "reg:squarederror",
          eta = sampled_grid$eta[i],
          alpha = sampled_grid$alpha[i],
          feature_selector = "shuffle",  # Use 'shuffle' for faster feature selection
          nthread = 1  # Explicitly set to 1 thread for better parallel control
        ),
        data = dtrain,
        nrounds = sampled_grid$nrounds[i],
        verbose = 0
      )
      
      preds <- predict(model, newdata = x_valid)
      cv_errors[fold] <- sqrt(mean((preds - y_valid)^2))
    }
    data.frame(index = i, RMSE = mean(cv_errors, na.rm = TRUE))
  }
  
  sampled_grid$RMSE <- res$RMSE
  best <- sampled_grid[which.min(sampled_grid$RMSE), ]
  message("Best XGBoost Linear parameters using parallel random search and rolling origin CV:")
  print(best)
  return(best)
}

# =======================================================================================================================================
# Macroeconomic Random Forest

# Full version
tune_MRF <- function(data.in,n){
  
  # Inputs:
  # data.in = dataset with target in first column and variables in other columns - should be a matrix
  # n = number of periods on which the hyper-parameter should be optimized
  # Output:
  # optim_param = vector of best-performing hyper-parameters in the following order:
  #           1 = n_trees (number of trees - default value = 50)
  #           2 = n_var (number of variables in the linear part)
  #           3 = lambda (ridge shrinkage parameter in the linear part - default value = 0.01)
  #           4 = re_meth (resampling method - default value = 2, goes from 0 to 4)
  #           5 = bl_size (block size for resampling - default value = 12)
  # NB: other variables might be optimized but are not looked at here for the sake of time
  
  print("Tuning macroeconomic random forest")
  
  # Start and end values
  n_tot <- nrow(data.in)
  n_test <- n_tot - n + 1
  
  # Define train and test samples
  data.tr <- data.in[-n_tot,]
  y_te <- data.tr[(n_test-1):(n_tot-1),1]
  data.tr[(n_test-1):(n_tot-1),1] <- NA
  
  #
  # Step 1: n_trees and n_var
  #
  
  # Summary dataframe
  summary <- data.frame(matrix(NA,
                               nrow = 1,
                               ncol = 3))
  colnames(summary) <- c("n_trees","n_var","RMSE")
  count_mod <- 1 
  
  m_min <- min(5,ncol(x))
  m_max <- max(ncol(x),9)
  
  for (n_trees in c(25,50)){
    for (n_var in seq(m_min,m_max,by=2)){
      
      # Set seed for reproducibility
      # NB: has to be the same as in Main.R
      set.seed(1234) 
      
      # Fit model on train
      eq_mrf <- MRF(data.tr,
                    B = n_trees,
                    x.pos = c(2:n_var),
                    oos.pos = c((nrow(data.tr)-11):nrow(data.tr)),
                    cheap.look.at.GTVPs = FALSE,
                    printb = FALSE)
      
      # Predict on test
      fit <- eq_mrf$pred
      
      # Write summary and advance counter
      summary[count_mod,1] <- n_trees
      summary[count_mod,2] <- n_var
      summary[count_mod,3] <- sqrt(sum((fit - y_te)^2)/length(y_te))
      count_mod <-  count_mod + 1
      
    }
  }
  
  # Select best hyper-parameter and return it
  summary %<>%
    arrange(RMSE)
  optim_trees <- summary[1,1]
  optim_var <- summary[1,2]
  
  
  #
  # Step 2: L2 regularization parameter
  #
  
  # Summary dataframe
  summary <- data.frame(matrix(NA,
                               nrow = 1,
                               ncol = 2))
  colnames(summary) <- c("lambda","RMSE")
  count_mod <- 1 
  
  for (lambda in c(0.001,0.01,0.1)){
    
    # Set seed for reproducibility
    # NB: has to be the same as in Main.R
    set.seed(1234) 
    
    # Fit model on train
    eq_mrf <- MRF(data.tr,
                  B = optim_trees,
                  x.pos = c(2:optim_var),
                  oos.pos = c((nrow(data.tr)-11):nrow(data.tr)),
                  cheap.look.at.GTVPs = FALSE,
                  printb = FALSE,
                  ridge.lambda = lambda)
    
    # Predict on test
    fit <- eq_mrf$pred
    
    # Write summary and advance counter
    summary[count_mod,1] <- lambda
    summary[count_mod,2] <- sqrt(sum((fit - y_te)^2)/length(y_te))
    count_mod <-  count_mod + 1
    
  }
  
  # Select best hyper-parameter and return it
  summary %<>%
    arrange(RMSE)
  optim_lambda <- summary[1,1]
  
  
  #
  # Step 3: re-sampling
  #
  
  # Summary dataframe
  summary <- data.frame(matrix(NA,
                               nrow = 1,
                               ncol = 3))
  colnames(summary) <- c("re_meth","bl_size","RMSE")
  count_mod <- 1 
  
  for (re_meth in c(2,4)){
    for (bl_size in c(6,12)){
      
      # Set seed for reproducibility
      # NB: has to be the same as in Main.R
      set.seed(1234) 
      
      # Fit model on train
      eq_mrf <- MRF(data.tr,
                    B = optim_trees,
                    x.pos = c(2:optim_var),
                    oos.pos = c((nrow(data.tr)-11):nrow(data.tr)),
                    cheap.look.at.GTVPs = FALSE,
                    printb = FALSE,
                    ridge.lambda = optim_lambda,
                    resampling.opt = re_meth,
                    block.size = bl_size)
      
      # Predict on test
      fit <- eq_mrf$pred
      
      # Write summary and advance counter
      summary[count_mod,1] <- re_meth
      summary[count_mod,2] <- bl_size
      summary[count_mod,3] <- sqrt(sum((fit - y_te)^2)/length(y_te))
      count_mod <-  count_mod + 1
      
    }
  }
  
  # Select best hyper-parameter and return it
  summary %<>%
    arrange(RMSE)
  optim_meth <- summary[1,1]
  optim_size <- summary[1,2]
  
  
  #
  # Step 4: return list of optimal parameters
  #
  
  optim_param <- c(optim_trees,
                   optim_var,
                   optim_lambda,
                   optim_meth,
                   optim_size)
  
  names(optim_param) <- c("ntrees",
                          "nvar",
                          "ridge",
                          "resampling_method",
                          "block_size")
  
  return(optim_param) 
  
}

# Faster version
tune_MRF_fast <- function(data.in,n){
  
  # Inputs:
  # data.in = dataset with target in first column and variables in other columns - should be a matrix
  # n = number of periods on which the hyper-parameter should be optimized
  # Output:
  # optim_param = vector of best-performing hyper-parameters in the following order:
  #           1 = n_trees (number of trees - default value = 50)
  # NB: other variables can be optimized in non-fast version (and even more in principle)
  
  print("Tuning macroeconomic random forest (fast version)")
  
  # Start and end values
  n_tot <- nrow(data.in)
  n_test <- n_tot - n + 1
  
  # Define train and test samples
  data.tr <- data.in[-n_tot,]
  y_te <- data.tr[(n_test-1):(n_tot-1),1]
  data.tr[(n_test-1):(n_tot-1),1] <- NA
  
  # Summary dataframe
  summary <- data.frame(matrix(NA,
                               nrow = 1,
                               ncol = 2))
  colnames(summary) <- c("n_var","RMSE")
  count_mod <- 1 
  
  m_min <- min(5,ncol(x))
  m_max <- max(ncol(x),9)
  
  for (n_var in seq(m_min,m_max,by=2)){
    
    # Set seed for reproducibility
    # NB: has to be the same as in Main.R
    set.seed(1234) 
    
    # Fit model on train
    eq_mrf <- MRF(data.tr,
                  x.pos = c(2:n_var),
                  oos.pos = c((nrow(data.tr)-11):nrow(data.tr)),
                  cheap.look.at.GTVPs = FALSE,
                  printb = FALSE)
    
    # Predict on test
    fit <- eq_mrf$pred
    
    # Write summary and advance counter
    summary[count_mod,1] <- n_var
    summary[count_mod,2] <- sqrt(sum((fit - y_te)^2)/length(y_te))
    count_mod <-  count_mod + 1
    
  }
  
  # Select best hyper-parameter and return it
  summary %<>%
    arrange(RMSE)
  optim_param <- summary[1,1]
  
  names(optim_param) <- c("nvar")
  return(optim_param) 
  
}
# =======================================================================================================================================
# LSTM

tune_LSTM <- function(x, y, initial_window, horizon, n_folds = 3, seed = 1234, lag = 2) {
  print("Tuning LSTM")
  set.seed(seed)
  
  # Define a grid of hyperparameters to search over.
  # Here we vary only the number of hidden units in each LSTM layer over the values 10, 20, and 40.
  grid <- expand.grid(
    units = c(10, 20, 40),
    dropout = 0.3,
    recurrent_dropout = 0.2,
    epochs = 30,
    batch_size = 30
  )
  grid$RMSE <- NA_real_
  
  n_total <- nrow(x)
  total_splits <- n_total - initial_window - horizon + 1
  if(total_splits < n_folds) n_folds <- total_splits
  
  # Select a set of rolling CV splits.
  splits <- sort(sample(1:total_splits, n_folds))
  
  for (i in 1:nrow(grid)) {
    cv_errors <- c()
    
    # For each rolling CV fold:
    for (fold in splits) {
      train_end <- initial_window + fold - 1
      train_indices <- 1:train_end
      valid_indices <- (train_end + 1):(train_end + horizon)
      
      # Extract training and validation sets.
      x_train <- x[train_indices, , drop = FALSE]
      y_train <- y[train_indices]
      x_valid <- x[valid_indices, , drop = FALSE]
      y_valid <- y[valid_indices]
      
      # Create training sequences: each sample has 'lag' time steps.
      num_train <- nrow(x_train) - lag + 1
      if(num_train < 1) next  # skip if not enough data
      x_train_seq <- array(NA, dim = c(num_train, lag, ncol(x_train)))
      for (j in 1:num_train) {
        x_train_seq[j, , ] <- as.matrix(x_train[j:(j+lag-1), ])
      }
      # The label for each sequence is taken from the last row of the sequence.
      y_train_seq <- y_train[lag:length(y_train)]
      
      # Create validation sequences similarly.
      num_valid <- nrow(x_valid) - lag + 1
      if(num_valid < 1) next
      x_valid_seq <- array(NA, dim = c(num_valid, lag, ncol(x_valid)))
      for (j in 1:num_valid) {
        x_valid_seq[j, , ] <- as.matrix(x_valid[j:(j+lag-1), ])
      }
      y_valid_seq <- y_valid[lag:length(y_valid)]
      
      # Build the LSTM model with 2 LSTM layers.
      model <- keras_model_sequential() %>%
        layer_lstm(
          units = grid$units[i],
          input_shape = c(lag, ncol(x_train)),
          dropout = grid$dropout[i],
          recurrent_dropout = grid$recurrent_dropout[i],
          return_sequences = TRUE
        ) %>%
        layer_lstm(
          units = grid$units[i],
          dropout = grid$dropout[i],
          recurrent_dropout = grid$recurrent_dropout[i]
        ) %>%
        layer_dense(units = 1)
      
      model %>% compile(
        loss = "mean_squared_error",
        optimizer = optimizer_adam(lr = 0.001)
      )
      
      # Callbacks for early stopping and reducing learning rate.
      early_stop <- callback_early_stopping(monitor = "val_loss", patience = 5)
      lr_reduce <- callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.5, patience = 3)
      
      history <- model %>% fit(
        x = x_train_seq,
        y = y_train_seq,
        epochs = grid$epochs[i],
        batch_size = grid$batch_size[i],
        verbose = 0,
        validation_data = list(x_valid_seq, y_valid_seq),
        callbacks = list(early_stop, lr_reduce)
      )
      
      preds <- model %>% predict(x_valid_seq)
      preds <- as.numeric(preds)
      rmse <- sqrt(mean((preds - y_valid_seq)^2))
      cv_errors <- c(cv_errors, rmse)
      
      rm(model, history, preds, x_train_seq, x_valid_seq)
      keras::k_clear_session()
      gc(full = TRUE)
    }
    
    grid$RMSE[i] <- mean(cv_errors, na.rm = TRUE)
  }
  
  best <- grid[which.min(grid$RMSE), ]
  message("Best LSTM parameters using rolling CV:")
  print(best)
  return(best)
}
