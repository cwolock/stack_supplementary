do_one <- function(n_train, n_test = 1000, estimator, dgp){

  dimension <- 10

  # training data
  data_gen <- generate_data(n = n_train*6,
                            truncation = "covariate",
                            direction = "prospective",
                            dgp = dgp)
  train <- data_gen$data
  indices <- sample(1:nrow(train), n_train)
  train <- train[indices,] # b/c of truncation,
  # generate more samples than needed, then randomly select n_train
  true_S_T <- data_gen$true_S_T

  # test data generated without truncation
  data_gen <- generate_data(n = n_test,
                            truncation = "none",
                            direction = "prospective",
                            dgp = dgp)
  test <- data_gen$data
  theo_quant <- round(quantile(test$Y[test$Delta == 1], probs = c(0.5, 0.75, 0.9)),
                      digits = 0)
  # benchmarks
  approx_times <- sort(unique(train$Y))#[train$Delta == 1]))
  benchmark_times <- seq(0.1, 100, by = 0.1)

  # calculate true survival function values
  true_df_uni <- matrix(NA, nrow = nrow(test), ncol = length(benchmark_times))
  for (i in 1:length(benchmark_times)){
    vals <- unlist(lapply(X = 1:nrow(test), FUN = function(x){
      true_S_T(as.numeric(test[x,1:dimension]), benchmark_times[i])
    }))
    true_df_uni[,i] <- vals
  }

  # tuning parameters
  tune <- list(ntrees = c(250, 500, 1000),
               max_depth = c(1,2),
               minobspernode = 1,
               shrinkage = 0.01)
  xgb_grid <- create.SL.xgboost(tune = tune)
  SL.library <- c("SL.mean", "SL.glm.interaction", "SL.earth",
                  "SL.gam", "SL.ranger", xgb_grid$names)

  start_time <- Sys.time()
  if (estimator == "stackG_fine_W"){ # global stacking, all times grid
    F_Y_1_grid <- sort(unique(c(0,train$Y[train$Delta == 1])))
    F_Y_0_grid <- sort(unique(c(0,train$Y[train$Delta == 0])))
    G_W_1_grid <- sort(unique(c(0,train$W[train$Delta == 1])))
    G_W_0_grid <- sort(unique(c(0,train$W[train$Delta == 0])))
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          time_grid_fit = list(F_Y_1_grid = F_Y_1_grid,
                                               F_Y_0_grid = F_Y_0_grid,
                                               G_W_1_grid = G_W_1_grid,
                                               G_W_0_grid = G_W_0_grid),
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
    est_df <- out$S_T_preds
  } else if (estimator == "stackG_medium_W"){
    F_Y_1_grid <- c(0,
                    sort(unique(quantile(train$Y[train$Delta == 1],
                                         probs = seq(0, 1,length.out = 40)))))
    F_Y_0_grid <- c(0,
                    sort(unique(quantile(train$Y[train$Delta == 0],
                                         probs = seq(0, 1,length.out = 40)))))
    G_W_1_grid <- c(0,
                    sort(unique(quantile(train$W[train$Delta == 1],
                                         probs = seq(0, 1,length.out = 40)))))
    G_W_0_grid <- c(0,
                    sort(unique(quantile(train$W[train$Delta == 0],
                                         probs = seq(0, 1,length.out = 40)))))
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          time_grid_fit = list(F_Y_1_grid = F_Y_1_grid,
                                               F_Y_0_grid = F_Y_0_grid,
                                               G_W_1_grid = G_W_1_grid,
                                               G_W_0_grid = G_W_0_grid),
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
    est_df <- out$S_T_preds
  } else if (estimator == "stackG_coarse_W"){
    F_Y_1_grid <- c(0,
                    sort(unique(quantile(train$Y[train$Delta == 1],
                                         probs = seq(0, 1,length.out = 10)))))
    F_Y_0_grid <- c(0,
                    sort(unique(quantile(train$Y[train$Delta == 0],
                                         probs = seq(0, 1,length.out = 10)))))
    G_W_1_grid <- c(0,
                    sort(unique(quantile(train$W[train$Delta == 1],
                                         probs = seq(0, 1,length.out = 10)))))
    G_W_0_grid <- c(0,
                    sort(unique(quantile(train$W[train$Delta == 0],
                                         probs = seq(0, 1,length.out = 10)))))
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          time_grid_fit = list(F_Y_1_grid = F_Y_1_grid,
                                               F_Y_0_grid = F_Y_0_grid,
                                               G_W_1_grid = G_W_1_grid,
                                               G_W_0_grid = G_W_0_grid),
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
    est_df <- out$S_T_preds
  } else if (estimator == "stackG_fine_Y"){ # global stacking, all times grid
    F_Y_1_grid <- sort(unique(c(0,train$Y[train$Delta == 1])))
    F_Y_0_grid <- sort(unique(c(0,train$Y[train$Delta == 0])))
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          time_grid_fit = list(F_Y_1_grid = F_Y_1_grid,
                                               F_Y_0_grid = F_Y_0_grid,
                                               G_W_1_grid = F_Y_1_grid,
                                               G_W_0_grid = F_Y_0_grid),
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
    est_df <- out$S_T_preds
  } else if (estimator == "stackG_medium_Y"){
    F_Y_1_grid <- c(0,
                    sort(unique(quantile(train$Y[train$Delta == 1],
                                         probs = seq(0, 1,length.out = 40)))))
    F_Y_0_grid <- c(0,
                    sort(unique(quantile(train$Y[train$Delta == 0],
                                         probs = seq(0, 1,length.out = 40)))))
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          time_grid_fit = list(F_Y_1_grid = F_Y_1_grid,
                                               F_Y_0_grid = F_Y_0_grid,
                                               G_W_1_grid = F_Y_1_grid,
                                               G_W_0_grid = F_Y_0_grid),
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
    est_df <- out$S_T_preds
  } else if (estimator == "stackG_coarse_Y"){
    F_Y_1_grid <- c(0,
                    sort(unique(quantile(train$Y[train$Delta == 1],
                                         probs = seq(0, 1,length.out = 10)))))
    F_Y_0_grid <- c(0,
                    sort(unique(quantile(train$Y[train$Delta == 0],
                                         probs = seq(0, 1,length.out = 10)))))
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          time_grid_fit = list(F_Y_1_grid = F_Y_1_grid,
                                               F_Y_0_grid = F_Y_0_grid,
                                               G_W_1_grid = F_Y_1_grid,
                                               G_W_0_grid = F_Y_0_grid),
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
    est_df <- out$S_T_preds
  } 
  end_time <- Sys.time()

  ### wrap things up
  squared_errors <- (est_df - true_df_uni)^2
  MSE_uni <- mean(squared_errors)
  landmark_indices <- which(round(benchmark_times, digits = 2) %in% theo_quant)
  landmark_estimates <- est_df[, landmark_indices]
  landmark_truth <- true_df_uni[, landmark_indices]
  landmark_sq_error <- (landmark_estimates - landmark_truth)^2
  landmark_MSE <- colSums(landmark_sq_error)/n_test
  output <- data.frame(MSE_uni = MSE_uni,
                       landmark_MSE_50 = landmark_MSE[1],
                       landmark_MSE_75 = landmark_MSE[2],
                       landmark_MSE_90 = landmark_MSE[3])
  output$dgp <- rep(dgp, nrow(output))
  output$n_train <- rep(n_train, nrow(output))
  output$estimator <- rep(estimator, nrow(output))
  runtime <- difftime(end_time, start_time, units = "secs")
  output$runtime <- rep(runtime, nrow(output))
  return(output)
}
