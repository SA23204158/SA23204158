#' @title Based on EM algorithm and boot algorithm, the survival analysis of censored data was carried out
#' @description Based on EM algorithm and boot algorithm, the survival analysis of censored data was carried out
#' @param x the results of the first type of experiment
#' @param s cut-off values in the first type of experiment
#' @param y the results of the second type of experiment
#' @param u the results of the third type of experiment
#' @param v the results of the third type of experiment
#' @param lambda0 the initial value of the iteration
#' @param max.it the maximum number of iterations
#' @param lambda parameters of the survival model
#' @param B number of bootstrap samples
#' @param alpha the level of significance of the confidence interval
#' @return  a list with maximum likelihood estimates for parameters, full values for censored data, and estimated expectations, variances, and confidence intervals
#' @export
EM <- function(x, s, y, u, v, lambda0 = 0.1, max.it = 1e4, tol = 1e-6){
  m <- length(x)
  r <- m - sum(x)
  n <- length(y)
  k <- length(u)
  lambda <- numeric(length = max.it)
  lambda[1] <- lambda0
  
  for (i in 2:max.it) {
    n1 <- m + n + k - 1
    n2 <- s * exp(-s * lambda[i - 1])
    d2 <- 1 - exp(-s * lambda[i - 1])
    n3 <- (u * exp(-lambda[i - 1] * u) - v * exp(-lambda[i - 1] * v))
    d3 <- exp(-lambda[i - 1] * u) - exp(-lambda[i - 1] * v)
    d1 <- n * mean(y) + (m - r) * (s + 1 / lambda[i - 1]) + r * (1 / lambda[i - 1] - n2 / d2) + k / lambda[i - 1] + sum(n3 / d3)
    lambda[i] <- n1 / d1
    
    # Check convergence
    if (abs(lambda[i] - lambda[i - 1]) < tol) {
      return(list(lambda = lambda[i], iteration = i))
    }
  }
  
  # If convergence criteria are not met, return the result at the maximum iteration
  return(list(lambda = lambda[max.it], iteration = max.it))
}

Supvalue <- function(x, s, y, u, v, lambda){
  m <- length(x)
  r <- m - sum(x)
  n <- length(y)
  k <- length(u)
  Y <- numeric(m + n + k)
  u1 <- runif(r ,0 , 1)
  u2 <- runif(m-r,0 , 1)
  for (i in 1:r) {
    Y[i] <- -1/lambda * log(1 - u1[i]*(1 - exp(-lambda*s)))
  }
  for (i in 1:(m-r)) {
    Y[r + i] <- s - log(u2[i]) / lambda
  }
  for (i in 1:n) {
    Y[m + i] <- y[i]
  }
  for (i in 1:k) {
    Y[m + n + i] <- ((u[i] + 1/lambda)*exp(-lambda*u[i]) - (v[i] + 1/lambda)*exp(-lambda*v[i])) / (exp(-lambda*u[i]) - exp(-lambda*v[i]))
  }
  return(Y)
}

Bootstrap <- function(x,B,alpha=0.05){
  lambdahat <-  1 / mean(x)
  lambdastar <- numeric(B)
  for(b in 1:B){
    xstar <- sample(x,replace=TRUE)
    lambdastar[b] <- 1 / mean(xstar)
  }
  mean.boot <- mean(lambdastar)
  bias.boot <- mean(lambdastar) - lambdahat
  se.boot   <- sd(lambdastar)
  
  normal.CI1 <- lambdahat - qnorm(1-0.5*alpha)*se.boot
  normal.CI2 <- lambdahat + qnorm(1-0.5*alpha)*se.boot
  basic.CI1 <- 2*lambdahat - unname(quantile(lambdastar, 1-0.5*alpha))
  basic.CI2 <- 2*lambdahat - unname(quantile(lambdastar, 0.5*alpha))
  percentile.CI1 <- unname(quantile(lambdastar, 0.5*alpha))
  percentile.CI2 <- unname(quantile(lambdastar, 1-0.5*alpha))
  
  
  return(c(list(mean.boot = mean.boot, se.boot = se.boot, bias.boot=bias.boot), 
           list(normal.CI1=normal.CI1,normal.CI2=normal.CI2), 
           list(basic.CI1=basic.CI1,basic.CI2=basic.CI2),
           list(percentile.CI1=percentile.CI1,percentile.CI2=percentile.CI2)))
}