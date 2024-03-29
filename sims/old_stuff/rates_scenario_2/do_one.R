do_one <- function(n_train, n_test=1000, rate, dgp){

  dimension <- 10

  # training data
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

  # test data
  data_gen <- generate_data(n = n_test,
                            truncation = "none",
                            direction = "prospective",
                            dgp = dgp)
  test <- data_gen$data
  theo_quant <- round(quantile(test$Y[test$Delta == 1],
                               probs = c(0.5, 0.75, 0.9)),
                      digits = 0)
  # benchmarks
  approx_times <- sort(unique(train$Y[train$Delta == 1])) # just a placeholder
  benchmark_times <- seq(0.1, 100, by = 0.1)

  # calculate true survival function values
  true_df_uni <- matrix(NA, nrow = nrow(test), ncol = length(benchmark_times))
  for (i in 1:length(benchmark_times)){
    vals <- unlist(lapply(X = 1:nrow(test), FUN = function(x){
      true_S_T(as.numeric(test[x,1:dimension]), benchmark_times[i])
    }))
    true_df_uni[,i] <- vals
  }

  # set up tuning parameters
  tune <- list(ntrees = c(250, 500, 1000),
               max_depth = c(1,2),
               minobspernode = 1,
               shrinkage = 0.01)
  xgb_grid <- create.SL.xgboost(tune = tune)
  SL.library <- c("SL.mean", "SL.glm.interaction", "SL.earth",
                  "SL.gam", "SL.ranger", xgb_grid$names)
  start_time <- Sys.time()
  approx_rates <- c("0.25", "0.333", "0.5", "0.667", "0.75", "1")
  if (rate == "1/4"){
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          bin_size = n_train^(-1/4),
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
  } else if (rate == "1/3"){ # global stacking, all times grid
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          bin_size = n_train^(-1/3),
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
  } else if (rate == "1/2"){ # global stacking 0.025 grid
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          bin_size = n_train^(-1/2),
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
  } else if (rate == "2/3"){ # global stacking, 0.1 grid
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          bin_size = n_train^(-2/3),
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
  } else if (rate == "3/4"){ # global stacking, 0.1 grid
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          bin_size = n_train^(-3/4),
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
  } else if (rate == "1"){ # global stacking, 0.1 grid
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          bin_size = n_train^(-1),
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
  }
  for (i in 1:length(approx_rates)){
    approx_rate_numeric <- as.numeric(approx_rates[i])
    approx_times <- quantile(train$Y, probs = seq(0, 1, length.out = n_train^approx_rate_numeric))
    preds <- predict(out,
                     time_grid_approx = approx_times,
                     newX = test[,1:dimension],
                     newtimes = benchmark_times)
    est_df <- preds$S_T_preds

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
    output$approx_rate <- approx_rates[i]
    if (i == 1){
      pooled_output <- output
    } else{
      pooled_output <- rbind(pooled_output, output)
    }
  }
  output <- pooled_output
  end_time <- Sys.time()
  output$dgp <- rep(dgp, nrow(output))
  output$n_train <- rep(n_train, nrow(output))
  output$rate <- rep(rate, nrow(output))
  runtime <- difftime(end_time, start_time, units = "secs")
  output$runtime <- rep(runtime, nrow(output))
  return(output)
}
