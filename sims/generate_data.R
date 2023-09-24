generate_data <- function(n = 500,
                          truncation = "none",
                          dgp = "leftskew",
                          direction = "prospective",
                          ph = FALSE,
			  cens = 0.25,
                          discretize = FALSE,
                          nbin = 10){
  if (!ph){
    if (dgp == "leftskew" & truncation == "none"){
      c_intercept <- 5.1
    } else if(dgp == "leftskew" & truncation != "none"){
      c_intercept <- 4.8
    } else if(dgp == "rightskew" & truncation == "none"){
      c_intercept <- 4.7
    } else if(dgp == "rightskew" & truncation != "none"){
      c_intercept <- 4.8
    }
  } else{
    if (dgp == "leftskew" & truncation == "none"){
      c_intercept <- 5.1
    } else if(dgp == "leftskew" & truncation != "none"){
      c_intercept <- 5.0
    } else if(dgp == "rightskew" & truncation == "none"){
      c_intercept <- 4.8
    } else if(dgp == "rightskew" & truncation != "none"){
      c_intercept <- 4.9
    }
  }
  
  train <- data.frame(X1 = runif(n, min = -1, max = 1),
                      X2 = runif(n, min = -1, max = 1),
                      X3 = (2*rbinom(n, 1, 0.5)-1),
                      X4 = (2*rbinom(n, 1, 0.5)-1),
                      X5 = rnorm(n, 0, 1),
                      X6 = rnorm(n, 0, 1),
                      X7 = rnorm(n, 0, 1),
                      X8 = rnorm(n, 0, 1),
                      X9 = rnorm(n, 0, 1),
                      X10 = rnorm(n, 0, 1))
  betaC <- c(c_intercept, 0.5, 0.5, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0) 
  hazC <- exp(betaC[1] + betaC[2]*train$X1 + betaC[3]*train$X2 + betaC[4]*train$X3 +
                betaC[5]*train$X4 + betaC[6]*train$X5)
  train$C <- rweibull(n = n, shape = 1.5, scale = hazC)
  
  if (!ph){
    betaT1 <- 2*c(0, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0,0.5, 0.5, 0.5)
    a <- exp(betaT1[1] + betaT1[2]*train$X1 + betaT1[3]*train$X2 + betaT1[4]*train$X3 +
               betaT1[5]*train$X4 + betaT1[6]*train$X5 + betaT1[12]*train$X1*train$X2 +
               betaT1[13]*train$X3*train$X4 + betaT1[14]*train$X5*train$X1) + 2
    b <- 2
    if (dgp == "leftskew"){
      train$T <- 100*rbeta(n = n, shape1 = a, shape2 = b)
    } else if (dgp == "rightskew"){
      train$T <- 100*rbeta(n = n, shape1 = b, shape2 = a)
    }
    
    true_S_T <- function(X, t){
      X <- c(1,X, X[1]*X[2], X[3]*X[4], X[5]*X[1]) # intercept
      a <- exp(betaT1 %*% X) + 2
      b <- 2
      if (dgp == "leftskew"){
        S_T <- pbeta(t/100, shape1 = a, shape2 = 2, lower.tail = FALSE)
      } else if (dgp == "rightskew"){
        S_T <- pbeta(t/100, shape1 = 2, shape2 = a, lower.tail = FALSE)
      }
      return(S_T)
    }
  } else if (ph){
    baseline_haz = function(t){
      a <- 3
      b <- 2
      if (dgp == "leftskew"){
        S_T <- pbeta(t/100, shape1 = a, shape2 = 2, lower.tail = FALSE)
      } else if (dgp == "rightskew"){
        S_T <- pbeta(t/100, shape1 = 2, shape2 = a, lower.tail = FALSE)
      }
      return(-log(S_T))
    }
    
    approx_grid <- seq(0, 99.9, by = 0.01)
    cuhaz_grid <- apply(X = matrix(approx_grid), MARGIN = 1, FUN = baseline_haz)
    
    cuhaz_mat <- do.call("rbind", replicate(n, cuhaz_grid, simplify = FALSE))
    
    betaT1 <- 1*c(0, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0)
    a <- exp(betaT1[1] + betaT1[2]*train$X1 + betaT1[3]*train$X2 + betaT1[4]*train$X3 +
               betaT1[5]*train$X4 + betaT1[6]*train$X5)
    
    a <- t(do.call("rbind", replicate(length(approx_grid), a, simplify = FALSE)))
    cuhaz_mat <- a * cuhaz_mat
    
    exp_sample <- rexp(n = n, rate = 1)
    
    inverse_baseline_haz <- function(i){
      if (exp_sample[i] > max(cuhaz_mat[i,])){
        inv <- max(cuhaz_mat[i,])
      } else{
        inv <- approx_grid[min(which(cuhaz_mat[i,] >= exp_sample[i]))]
      }
      return(inv)
    }
    
    beta_sample <- apply(X = matrix(1:n), MARGIN = 1, FUN = inverse_baseline_haz)
    
    train$T <- beta_sample
    
    true_S_T <- function(X, t){
      a <- 3
      b <- 2
      mult <- exp(betaT1[-1] %*% X)
      if (dgp == "leftskew"){
        H_T <- -log(pbeta(t/100, shape1 = a, shape2 = 2, lower.tail = FALSE))
      } else if (dgp == "rightskew"){
        H_T <- -log(pbeta(t/100, shape1 = 2, shape2 = a, lower.tail = FALSE))
      }
      HX_T <- mult * H_T
      S_T <- exp(-HX_T)
      return(S_T)
    }
  }
  
  train$Delta <- ifelse(train$T <= train$C, 1, 0)
  train$Y <- ifelse(train$T <= train$C, train$T, train$C)
  
  # truncation
  if (truncation == "uniform"){
    train$W <- runif(n, 0, 100)
  } else if (truncation == "none"){
    if (direction == "prospective"){
      train$W <- 0
    } else if (direction == "retrospective"){
      train$W <- 100
    }
  } else if (truncation == "covariate"){
    a <- ifelse(train$X1 < 0, 1, 1.5)
    b <- ifelse(train$X1 > 0, 1, 1.5)
    train$W <- 100*rbeta(n = n, shape1 = a, shape2 = b)
  }
  
  if (direction == "prospective"){
    train <- train[train$Y >= train$W,]
  } else if (direction == "retrospective"){
    train <- train[train$T <= train$W,]
  }
  
  if (discretize){
    time_grid <- c(0, seq((100/nbin), 100, length.out = nbin))
    discretize_time <- function(t, forward = TRUE){
      if (t <= max(time_grid) & forward){
        time <- time_grid[min(which(time_grid >= t))]
      } else if (t <= max(time_grid) & !forward){
        time <- time_grid[max(which(time_grid <= t))]
      } else{
        time <- max(time_grid)
      }
      return(time)
    }
    train$Y <- apply(X = matrix(train$Y), MARGIN = 1, FUN = discretize_time)
    train$T <- apply(X = matrix(train$T), MARGIN = 1, FUN = discretize_time)
    train$C <- apply(X = matrix(train$C), MARGIN = 1, FUN = discretize_time)
    train$W <- apply(X = matrix(train$W), MARGIN = 1, FUN = discretize_time,
                     forward = FALSE)
  }
  
  return(list(data = train, true_S_T = true_S_T))
}
