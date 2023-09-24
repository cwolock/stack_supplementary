# setup
library(survival)
library(tidyverse)
library(survML)
set.seed(123)

# data are available in survival package
data(nwtco)

### variables of interest
# edrel - time to event outcome
# rel - event indicator
# instit - local histology (binary)
# histol - central histology (binary)
# age - age in months
# stage - disease stage (4 stages, use dummy variables)
# study - which study (binary)
## not including id number or included in subcohort indicator

# wrangle data
nwtco$stage2 <- ifelse(nwtco$stage == 2, 1, 0)
nwtco$stage3 <- ifelse(nwtco$stage == 3, 1, 0)
nwtco$stage4 <- ifelse(nwtco$stage == 4, 1, 0)

# organize relevant variables
X <- nwtco %>% select(instit, histol, age, stage2, stage3, stage4, study)
Y <- nwtco$edrel
Delta <- nwtco$rel

# landmark times
landmark_t <- quantile(Y[Delta == 1], probs = c(0.5, 0.75, 0.9))

# set up folds for CV estimation of prediction error
nfolds <- 5
folds <- sample(rep(seq_len(nfolds), length = length(Y)))

# store results here
stackG_briers <- matrix(NA, nrow = length(landmark_t), ncol = nfolds)
stackL_briers <- matrix(NA, nrow = length(landmark_t), ncol = nfolds)
cox_briers <- matrix(NA, nrow = length(landmark_t), ncol = nfolds)
km_briers <- matrix(NA, nrow = length(landmark_t), ncol = nfolds)
naive_briers <- matrix(NA, nrow = length(landmark_t), ncol = nfolds)
rf_briers <- matrix(NA, nrow = length(landmark_t), ncol = nfolds)

for (k in 1:nfolds){
  # organize test and train data
  test_folds <- k
  train_folds <- which(1:nfolds != k)
  train_indices <- which(folds %in% train_folds)
  test_indices <- which(folds %in% test_folds)
  train <- data.frame(Y = Y[train_indices],
                      Delta = Delta[train_indices],
                      X[train_indices,])
  test <- data.frame(Y = Y[test_indices],
                     Delta = Delta[test_indices],
                     X[test_indices,])

  # set up super learner stuff
  tune <- list(ntrees = c(250, 500, 1000),
               max_depth = c(1,2),
               minobspernode = 1,
               shrinkage = 0.01)
  xgb_grid <- create.SL.xgboost(tune = tune)
  SL.library <- c("SL.mean", "SL.glm.interaction", "SL.earth",
                  "SL.gam", "SL.ranger", xgb_grid$names)

  approx_times <- quantile(sort(unique(Y)), probs = seq(0, 1, by = 0.01))

  # fit random forest
  formula <- paste0("survival::Surv(entry, time, event) ~ ",
                    paste0( names(train[,-c(1,2)]), collapse = "+"))
  rf_fit <- LTRCforests::ltrccif(as.formula(formula),
                                 data = data.frame(time=train$Y,
                                                   event=train$Delta,
                                                   entry = 0,
                                                   train[,-c(1,2)]),
                                 mtry = ceiling(sqrt(ncol(X))))
  rf_pred <- t(LTRCforests::predictProb(rf_fit,
                                        newdata = data.frame(entry = 0,
                                                             time = test$Y,
                                                             event = test$Delta,
                                                             test[,-c(1,2)]),
                                        time.eval = approx_times)$survival.probs)

  # fit global stacking
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
  stackG_pred <- stackG_out$S_T_preds

  # fit local stacking
  stackL_out <- survML::stackL(time = train$Y,
                               event = train$Delta,
                               X = train[,-c(1,2)],
                               newX = test[,-c(1,2)],
                               newtimes = approx_times,
                               bin_size = 0.025,
                               time_basis = "continuous",
                               SL_control = list(SL.library = SL.library,
                                                 V = 5))
  stackL_pred <- stackL_out$S_T_preds

  # fit cox model
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

  # fit Kaplan Meier for censoring weights
  cens_km <- survival::survfit(
    survival::Surv(time, event)~1,
    data = as.data.frame(cbind(time=test$Y, event=1-test$Delta, test[,-c(1:2)])))
  cens_pred <- matrix(stats::stepfun(cens_km$time, c(1,cens_km$surv), right = FALSE)(approx_times),
                      nrow=nrow(test), ncol = length(approx_times), byrow=TRUE)

  # fit Kaplan-Meier for estimates as a comparator method
  event_km <- survival::survfit(
    survival::Surv(time, event)~1,
    data = as.data.frame(cbind(time=train$Y, event=train$Delta, train[,-c(1:2)])))
  event_km_pred <- matrix(stats::stepfun(event_km$time, c(1,event_km$surv), right = FALSE)(approx_times),
                          nrow=nrow(test), ncol = length(approx_times), byrow=TRUE)

    # function to compute Brier score
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

  # loop through landmark times
  for (j in 1:length(landmark_t)){
    # fit naive estimator ignoring censoring
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

    # compute brier scores
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
    rf_brier <- apply(X = matrix(1:nrow(test)),
                      MARGIN = 1,
                      FUN = compute_brier,
                      t = landmark_t[j],
                      S_T_preds = rf_pred)

    stackG_briers[j,k] <- mean(stackG_brier)
    stackL_briers[j,k] <- mean(stackL_brier)
    cox_briers[j,k] <- mean(cox_brier)
    km_briers[j,k] <- mean(km_brier)
    naive_briers[j,k] <- mean(naive_brier)
    rf_briers[j,k] <- mean(rf_brier)
  }
}

results <- data.frame(stackG = rowMeans(stackG_briers),
                      stackL = rowMeans(stackL_briers),
                      cox = rowMeans(cox_briers),
                      naive = rowMeans(naive_briers),
                      rf = rowMeans(rf_briers),
                      km = rowMeans(km_briers))

# normalize by KM brier score
results <- results/results$km
saveRDS(results, "nwtco_results.rds")
