pre_select <- function(data_real,ii,horizon,select_method,n_var){
  smpl <- head(data_real,ii-horizon-3) %>%
    select(-date,
           -L1st_target,
           -L2nd_target) %>%
    drop_na()
    
    if(select_method==2){
      order <- smpl %>%
        mutate(across(-target,~cor(.,
                                   target,
                                   use = "pairwise.complete.obs",
                                   method  = "pearson"))) %>%
        select(-target) %>%
        distinct() %>%
        t() %>%
        as.data.frame() %>%
        rename(corr=V1) %>%
        mutate(variable = rownames(.))
      rownames(order) <- NULL # Reset row names to avoid conflicts
      
      out_sel <- order %>%
        mutate(corr = abs(corr)) %>%
        arrange(desc(corr)) %>%
        head(n_var)
      
      var_sel <- out_sel$variable
      
      
    }else if(select_method==1){
      
      # Creating batches for LARS loop
      batch <- smpl
      count <- 0
      order <- data.frame(V1=c(0))
      
      # LARS loop (have to run in different batchs if number of variables is too high)
      while (count < n_var) {
        
        # Running LARS equation
        x <- as.matrix(select(batch,-target))
        y <- as.matrix(select(batch,target))
        eq <- lars(x = x,
                   y = y,
                   type="lar")
        
        # Ordering the variables
        out <- as.data.frame(coef(eq)) %>%
          summarise(across(everything(),~sum(.==0))) %>%
          t() %>%
          as.data.frame() %>%
          rownames_to_column('name') %>%
          arrange(V1) %>%
          group_by(V1) %>%
          filter(n()==1) %>%
          column_to_rownames('name')
        
        order <- rbind(order,out)
        
        # Deleting variables already ordered from the sample
        var <- row.names(out)
        batch %<>%
          select(!all_of(var))
        
        # Checking if nrow(out) = 0 to avoid infinite loops
        if (nrow(out)==0){
          # Putting all remaining variables at the end
          x <- as.matrix(select(batch,-target))
          n_end <- ncol(x)
          out <- data.frame(V1=rep(1,n_end))
          rownames(out) <- colnames(x)
          order <- rbind(order,out)
          count <- count + n_end
          print("Warning: Using special procedure for LARS")
        }else{
          # Just updating the count      
          count <- count + nrow(out)
        } 
      }
      
      # Getting the list of selected variables
      out_sel <- out %>%
        head(n_var)
      var_sel <- row.names(out_sel)
      
    }else if(select_method==0){
      
      order <- smpl %>%
        select(-target)
      
      out_sel <- order
      
      var_sel <- colnames(out_sel)
      
    }else if(select_method==3){
      
      init <- head(data_real,ii-horizon-3) %>%
        select(-date) %>%
        drop_na()
      
      list_var <- colnames(select(smpl,-target))
      order <- data.frame(name_var=NA,
                          tstat=NA)
      
      for (v in seq_along(list_var)){
        
        data_eq <- init %>%
          select(target,
                 L1st_target,
                 L2nd_target,
                 list_var[v])
        
        eq <- lm(target ~ .,
                 data = data_eq,
                 na.action = "na.omit")
        
        order[v,1] <- list_var[v]
        order[v,2] <- coef(summary(eq))[4, "t value"]
      }
      
      # Order by t-stat (higher to lower)
      order %<>%
        arrange(desc(tstat))
      
      # Getting the list of selected variables
      out_sel <- order %>%
        head(n_var)
      var_sel <- out_sel$name_var
      
    }else if(select_method==4){
      
      # Hyper-parameters of BMA
      # Max number of variables (maxNvar) should not be > 30 or code will crash
      # NB: for fast BMA, the threshold of posterior probe and number of iterations have been determined empirically
      # This corresponds to value balancing time consumption with accuracy 
      # I.e. after 20 iterations, there are already 20 variables with posterior probability above 20
      # Continuing would just mean removing the variable with lowest probability and testing individually remaining ones (which can be inefficient)
      max_bma <- 20
      th_bma <- 20
      niter_fast_bma <- 20
      
      # Inputs to BMA
      x_all <- select(smpl,-target)
      y <- as.matrix(select(smpl,target))
      var_sel <- c()
      n_var_sel <- 0
      
      while(n_var_sel < n_var){
        
        # Get number of variables in this batch
        var_left <- n_var - n_var_sel
        max_batch <- min(max_bma,var_left)
        
        # Get variables for the batch
        x <- x_all %>%
          select(-all_of(var_sel)) %>%
          as.matrix()
        
        # If max_batch = 1, the do the selection based on sortedX (which is itself based on BMA)
        if(max_batch==1){
          bma <- iBMA.bicreg(x, 
                             y, 
                             thresProbne0 = th_bma,
                             verbose = TRUE,
                             maxNvar = max_batch+1,
                             nIter = 1)
          # ---- SHORT-CIRCUIT FOR TOO FEW VARIABLES ----
          if (length(bma$currentSet) == 0) {
            message("No more variables to add from BMA (currentSet is empty). Stopping early.")
            break
          }
          
          sel_fast <- bma$sortedX
          var_sel <- c(var_sel,colnames(sel_fast)[1])
          
        }else{
          if(fast_bma==1){
            
            # In fast BMA, the number of iterations (nIter) is set at a low value
            # For large datasets, this means some variables might NOT be considered by BMA but:
            # (1) variables are sorted prior to the BMA so after some iterations, variables are likely uninformative
            # (2) empirically this is confirmed: no variables enter the model after some iterations
            # (3) a large nIter results in a highly time-consuming procedure (for little to no gains given 1 and 2)
            bma <- iBMA.bicreg(x, 
                               y, 
                               thresProbne0 = th_bma,
                               verbose = TRUE,
                               maxNvar = max_batch,
                               nIter = niter_fast_bma)
            # ---- SHORT-CIRCUIT FOR TOO FEW VARIABLES ----
            if (length(bma$currentSet) == 0) {
              message("No more variables to add from BMA (currentSet is empty). Stopping early.")
              break
            }
            # Select intermediary selection of variables
            sel_fast <- bma$sortedX %>%
              as.data.frame() %>%
              select(all_of(bma$currentSet))
            var_sel <- c(var_sel,colnames(sel_fast))
            
            
          }else if(fast_bma==0){
            
            # Run BMA
            bma <- iBMA.bicreg(x, 
                               y, 
                               thresProbne0 = th_bma,
                               verbose = TRUE,
                               maxNvar = max_batch,
                               nIter = ncol(x))
            # ---- SHORT-CIRCUIT FOR TOO FEW VARIABLES ----
            if (length(bma$currentSet) == 0) {
              message("No more variables to add from BMA (currentSet is empty). Stopping early.")
              break
            }
            
            # Add results to var_sel
            var_sel <- c(var_sel,bma$bma$namesx)
            
          }
        }
        
        # Count number of var selected
        n_var_sel <- length(var_sel)
        
      }
      
      # Outputs
      order <- NA
      out_sel <- NA
      
    }
    
    return(var_sel)
}
