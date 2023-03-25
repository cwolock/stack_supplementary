do_one <- function(n_train, n_test = 1000, estimator, dgp, cens){

  dimension <- 10

  # training data
  data_gen <- generate_data(n = n_train*6,
                            truncation = "covariate",
                            direction = "prospective",
                            dgp = dgp, 
                            cens = cens)
  train <- data_gen$data
  censoring_rate <- mean(train$Delta)
  truncation_rate <- 1 - nrow(train)/(6*n_train)
  indices <- sample(1:nrow(train), n_train)
  train <- train[indices,] # b/c of truncation,
  # generate more samples than needed, then randomly select n_train
  true_S_T <- data_gen$true_S_T

  # test data generated without truncation
  data_gen <- generate_data(n = n_test,
                            truncation = "none",
                            direction = "prospective",
                            dgp = dgp,
                            cens = cens)
  test <- data_gen$data
  theo_quant <- round(quantile(test$Y[test$Delta == 1], 
                               probs = c(0.25, 0.5, 0.75, 0.9, 0.95)),
                      digits = 0)
  # benchmarks
  approx_times <- sort(unique(train$Y))
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
  if (estimator == "stackG_fine"){ # global stacking, all times grid
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
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
    fits <- names(out$fits[!sapply(out$fits, is.null)])
    algos <- list()
    weights <- list()
    for (i in 1:length(fits)){
      curr_fit <- fits[i]
      curr_algos <- names(out$fits[[curr_fit]]$reg.object$coef)
      curr_weights <- out$fits[[curr_fit]]$reg.object$coef
      algos[[i]] <- curr_algos
      weights[[i]] <- curr_weights
    }
    algos <- unlist(algos)
    weights <- unlist(weights)
    fits <- rep(fits, each = length(algos)/length(fits))
  } else if (estimator == "stackG_medium"){ # global stacking, 0.025 grid
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          bin_size = 0.025,
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
    est_df <- out$S_T_preds
    fits <- names(out$fits[!sapply(out$fits, is.null)])
    algos <- list()
    weights <- list()
    for (i in 1:length(fits)){
      curr_fit <- fits[i]
      curr_algos <- names(out$fits[[curr_fit]]$reg.object$coef)
      curr_weights <- out$fits[[curr_fit]]$reg.object$coef
      algos[[i]] <- curr_algos
      weights[[i]] <- curr_weights
    }
    algos <- unlist(algos)
    weights <- unlist(weights)
    fits <- rep(fits, each = length(algos)/length(fits))
  } else if (estimator == "stackG_coarse"){ # global stacking, 0.1 grid
    out <- survML::stackG(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          time_grid_approx = approx_times,
                          bin_size = 0.1,
                          time_basis = "continuous",
                          surv_form = "exp",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
    est_df <- out$S_T_preds
    fits <- names(out$fits[!sapply(out$fits, is.null)])
    algos <- list()
    weights <- list()
    for (i in 1:length(fits)){
      curr_fit <- fits[i]
      curr_algos <- names(out$fits[[curr_fit]]$reg.object$coef)
      curr_weights <- out$fits[[curr_fit]]$reg.object$coef
      algos[[i]] <- curr_algos
      weights[[i]] <- curr_weights
    }
    algos <- unlist(algos)
    weights <- unlist(weights)
    fits <- rep(fits, each = length(algos)/length(fits))
  } else if (estimator ==  "stackL_fine"){ # local stacking, all times grid
    out <- survML::stackL(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          bin_size = NULL,
                          time_basis = "continuous",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
    est_df <- out$S_T_preds
    fits <- rep("stack_fit", length(out$fit$coef))
    algos <- names(out$fit$coef)
    weights <- out$fit$coef
  } else if(estimator == "stackL_medium"){ # local stacking, 0.025 grid
    out <- survML::stackL(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          bin_size = 0.025,
                          time_basis = "continuous",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
    est_df <- out$S_T_preds
    fits <- rep("stack_fit", length(out$fit$coef))
    algos <- names(out$fit$coef)
    weights <- out$fit$coef
  } else if(estimator == "stackL_coarse"){ # local stacking, 0.1 grid
    out <- survML::stackL(time = train$Y,
                          event = train$Delta,
                          entry = train$W,
                          X = train[,1:dimension],
                          newX = test[,1:dimension],
                          newtimes = benchmark_times,
                          bin_size = 0.1,
                          time_basis = "continuous",
                          SL_control = list(SL.library = SL.library,
                                            V = 5))
    est_df <- out$S_T_preds
    fits <- rep("stack_fit", length(out$fit$coef))
    algos <- names(out$fit$coef)
    weights <- out$fit$coef
  } else if (estimator == "coxph"){
    fit <- survival::coxph(
      survival::Surv(entry, time, event) ~ .,
      data = as.data.frame(cbind(time=train$Y, 
                                 entry = train$W, 
                                 event=train$Delta, 
                                 train[1:dimension]))
    )
    pred <- t(summary(survival::survfit(fit,
                                        newdata=test[,1:dimension],
                                        se.fit = FALSE,
                                        conf.int = FALSE),
                      times=benchmark_times)$surv)
    pred <- pred[,-ncol(pred)]
    if(ncol(pred) < length(benchmark_times)) {
      pred <- cbind(pred, matrix(pred[,ncol(pred)], 
                                 nrow=nrow(pred), 
                                 ncol=length(benchmark_times) - ncol(pred)))
    }
    est_df <- pred
    fits <- NA
    algos <- NA
    weights <- NA
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
                       landmark_MSE_25 = landmark_MSE[1],
                       landmark_MSE_50 = landmark_MSE[2],
                       landmark_MSE_75 = landmark_MSE[3],
                       landmark_MSE_90 = landmark_MSE[4],
                       landmark_MSE_95 = landmark_MSE[5])
  output$dgp <- rep(dgp, nrow(output))
  output$cens <- rep(cens, nrow(output))
  output$n_train <- rep(n_train, nrow(output))
  output$estimator <- rep(estimator, nrow(output))
  output$censoring_rate <- rep(censoring_rate, nrow(output))
  output$truncation_rate <- rep(truncation_rate, nrow(output))
  output$fits <- paste(fits, collapse = ",")
  output$algos <- paste(algos, collapse = ",")
  output$weights <- paste(weights, collapse = ",")
  runtime <- difftime(end_time, start_time, units = "secs")
  output$runtime <- rep(runtime, nrow(output))
  return(output)

}
