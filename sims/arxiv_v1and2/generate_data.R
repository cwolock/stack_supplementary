generate_data <- function(n = 500,
                          truncation = "none",
                          dgp = "leftskew",
                          direction = "prospective"){

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

  betaC <- c(-5.5, 0.5, 0.5, 0.2, 0.2, 0.2,
             0, 0, 0, 0, 0) #X1 and X2 have equal effects on censoring
  hazC <- exp(betaC[1] + betaC[2]*train$X1 + betaC[3]*train$X2 + betaC[4]*train$X3 +
                betaC[5]*train$X4 + betaC[6]*train$X5)
  train$C <- rexp(n = n, rate = hazC)
  betaT1 <- 2*c(0, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0,0.5, 0.5, 0.5)
  a <- exp(betaT1[1] + betaT1[2]*train$X1 + betaT1[3]*train$X2 + betaT1[4]*train$X3 +
             betaT1[5]*train$X4 + betaT1[6]*train$X5 + betaT1[7]*train$X1*train$X2 +
             betaT1[8]*train$X3*train$X4 + betaT1[9]*train$X5*train$X1) + 2
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

  }

  if (direction == "prospective"){
    train <- train[train$Y >= train$W,]
  } else if (direction == "retrospective"){
    train <- train[train$T <= train$W,]
  }
  return(list(data = train, true_S_T = true_S_T))
}
