do_one <- function(n_train, n_test = 1000, estimator, dgp){

  dimension <- 10

  # training data
  data_gen <- generate_data(n = n_train,
                            truncation = "none",
                            direction = "prospective",
                            dgp = dgp)
  train <- data_gen$data
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
  if (estimator == "stackG_all"){ # global stacking, all times grid
    approx_times <- sort(unique(train$Y))#[train$Delta == 1]))
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          bin_size = NULL,
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
    est_df <- out$S_T_preds
  } else if (estimator == "stackG_250"){
    approx_times <- quantile(train$Y, probs = seq(0, 1, length.out = 250))
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          bin_size = NULL,
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
    est_df <- out$S_T_preds
  } else if (estimator == "stackG_100"){
    approx_times <- quantile(train$Y, probs = seq(0, 1, length.out = 100))
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          bin_size = NULL,
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
    est_df <- out$S_T_preds
  } else if (estimator == "stackG_50"){ # global stacking, all times grid
    approx_times <- quantile(train$Y, probs = seq(0, 1, length.out = 50))
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          bin_size = NULL,
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
