.libPaths(c(
  "/home/cwolock/R_lib",
  .libPaths()
))

set.seed(123)

library(survival)
library(tidyverse)
library(survML)

metabric <- readRDS("/home/cwolock/stack_supplementary/public_data_analysis/metabric.rds")
metabric <- metabric %>% na.omit()

X <- metabric[,-c(1,2)]
Y <- metabric[,1]
Delta <- metabric[,2]

landmark_t <- quantile(Y[Delta == 1], probs = c(0.5, 0.75, 0.9))

folds <- sample(rep(seq_len(5), length = length(Y)))

stackG_briers <- matrix(NA, nrow = length(landmark_t), ncol = 5)
stackL_briers <- matrix(NA, nrow = length(landmark_t), ncol = 5)
cox_briers <- matrix(NA, nrow = length(landmark_t), ncol = 5)
km_briers <- matrix(NA, nrow = length(landmark_t), ncol = 5)
naive_briers <- matrix(NA, nrow = length(landmark_t), ncol = 5)

for (k in 1:5){
  test_folds <- k
  train_folds <- which(1:5 != k)
  
  train_indices <- which(folds %in% train_folds)
  test_indices <- which(folds %in% test_folds)
  
  train <- data.frame(Y = Y[train_indices],
                      Delta = Delta[train_indices],
                      X[train_indices,])
  
  test <- data.frame(Y = Y[test_indices],
                     Delta = Delta[test_indices],
                     X[test_indices,])
  
  tune <- list(ntrees = c(250, 500, 1000),
               max_depth = c(1,2),
               minobspernode = 1,
               shrinkage = 0.01)
  xgb_grid <- create.SL.xgboost(tune = tune)
  SL.library <- c("SL.mean", "SL.glm.interaction", "SL.earth",
                  "SL.gam", "SL.ranger", xgb_grid$names)
  
  approx_times <- quantile(sort(unique(Y)), probs = seq(0, 1, by = 0.01))
  
  stackG_out <- survML::stackG(time = train$Y,
                               event = train$Delta,
                               X = train[,-c(1,2)],
                               newX = test[,-c(1,2)],
                               newtimes = approx_times,
                               time_grid_approx = approx_times,
                               bin_size = 0.025,
                               time_basis = "continuous",
                               surv_form = "PI",
                               SL_control = list(SL.library = SL.library,
                                                 V = 5))
  
  stackL_out <- survML::stackL(time = train$Y,
                               event = train$Delta,
                               X = train[,-c(1,2)],
                               newX = test[,-c(1,2)],
                               newtimes = approx_times,
                               bin_size = 0.025,
                               time_basis = "continuous",
                               SL_control = list(SL.library = SL.library,
                                                 V = 5))
  
  stackG_pred <- stackG_out$S_T_preds
  stackL_pred <- stackL_out$S_T_preds
  
  cox_out <-  survival::coxph(
    survival::Surv(time, event) ~ .,
    data = as.data.frame(cbind(time=train$Y, event=train$Delta, train[,-c(1:2)]))
  )
  cox_pred <- t(summary(survival::survfit(cox_out,
                                          newdata=test[,-c(1,2)],
                                          se.fit = FALSE,
                                          conf.int = FALSE),
                        times=approx_times)$surv)
  cox_pred <- cox_pred[,-ncol(cox_pred)]
  if(ncol(cox_pred) < length(approx_times)) {
    cox_pred <- cbind(cox_pred, matrix(cox_pred[,ncol(cox_pred)],
                                       nrow=nrow(cox_pred),
                                       ncol=length(approx_times) - ncol(cox_pred)))
  }
  
  cens_km <- survival::survfit(
    survival::Surv(time, event)~1,
    data = as.data.frame(cbind(time=test$Y, event=1-test$Delta, test[,-c(1:2)])))
  cens_pred <- matrix(stats::stepfun(cens_km$time, c(1,cens_km$surv), 
                                     right = FALSE)(approx_times),
                      nrow=nrow(test), ncol = length(approx_times), byrow=TRUE)
  
  event_km <- survival::survfit(
    survival::Surv(time, event)~1,
    data = as.data.frame(cbind(time=train$Y, event=train$Delta, train[,-c(1:2)])))
  event_km_pred <- matrix(stats::stepfun(event_km$time, c(1,event_km$surv), 
                                         right = FALSE)(approx_times),
                          nrow=nrow(test), ncol = length(approx_times), byrow=TRUE)
  
  compute_brier <- function(i, t, S_T_preds){
    t_index <- which.min(abs(approx_times - t))
    Y_index <- which.min(abs(approx_times - test$Y[i]))
    pred_i <- S_T_preds[i, t_index]
    if (test$Y[i] <= t & test$Delta[i] == 1){
      sq_err <- (0 - pred_i)^2 / cens_pred[i, Y_index]
    } else if (test$Y[i] > t){
      sq_err <- (1 - pred_i)^2 / cens_pred[i, t_index]
    } else{
      sq_err <- 0
    }
    return(sq_err)
  }
  
  for (j in 1:length(landmark_t)){
    outcome <- as.numeric(train$Y > landmark_t[j])
    design <- train[,-c(1,2)]
    naive <- SuperLearner::SuperLearner(Y = outcome,
                                        X = design,
                                        family = stats::binomial(),
                                        SL.library = SL.library,
                                        verbose = FALSE,
                                        cvControl = list(V = 5))
    naive_preds <- stats::predict(naive, newdata = test[,-c(1,2)])$pred
    
    naive_preds <- matrix(rep(naive_preds, length(approx_times)), 
                          nrow = nrow(test), ncol = length(approx_times))
    
    stackG_brier <- apply(X = matrix(1:nrow(test)),
                          MARGIN = 1,
                          FUN = compute_brier,
                          t = landmark_t[j],
                          S_T_preds = stackG_pred)
    
    stackL_brier <- apply(X = matrix(1:nrow(test)),
                          MARGIN = 1,
                          FUN = compute_brier,
                          t = landmark_t[j],
                          S_T_preds = stackL_pred)
    cox_brier <- apply(X = matrix(1:nrow(test)),
                       MARGIN = 1,
                       FUN = compute_brier,
                       t = landmark_t[j],
                       S_T_preds = cox_pred)
    
    km_brier <- apply(X = matrix(1:nrow(test)), 
                      MARGIN = 1, 
                      FUN = compute_brier, 
                      t = landmark_t[j], 
                      S_T_preds = event_km_pred)
    
    naive_brier <- apply(X = matrix(1:nrow(test)), 
                         MARGIN = 1, 
                         FUN = compute_brier, 
                         t = landmark_t[j], 
                         S_T_preds = naive_preds)
    
    stackG_briers[j,k] <- mean(stackG_brier)
    stackL_briers[j,k] <- mean(stackL_brier)
    cox_briers[j,k] <- mean(cox_brier)
    km_briers[j,k] <- mean(km_brier)
    naive_briers[j,k] <- mean(naive_brier)
  }
}

results <- data.frame(stackG = rowMeans(stackG_briers),
                      stackL = rowMeans(stackL_briers),
                      cox = rowMeans(cox_briers),
                      naive = rowMeans(naive_briers),
                      km = rowMeans(km_briers))

saveRDS(results, "metabric_results.rds")

