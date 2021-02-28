# European basket call option with 3 assets with equal weights (qi = 1/3)
S01 = 1.3     # Initial price of asset 1
S02 = 1.5     # Initial price of asset 2
S03 = 1.7     # Initial price of asset 3
S0 = c(S01, S02, S03) # Initial values

K = 1.1   # Strike Price
T = 1     # Maturity of the option
M = 10000   # Number of simulations

r = 0.05 # Interest rate on bank account

# Standard deviations of each assets are 0.2, 0.35, and 0.15 respectively
# Correlation between asset 1, 2, 3 are 0.1, 0.2, and 0.3 respectively
sigma = matrix(c(0.2, 0.35, 0.15), ncol = 3)   
corr = matrix(c(1,0.1,0.3, 0.1,1,0.2, 0.3,0.2,1), nrow = 3, 
              byrow = TRUE)  # Correlation matrix


## Price European basket call option using naive Monte Carlo
BasketPrice = function(InAssVal, r, sigma, corr, M_sim, K, T) {
  
  dim = nrow(corr)              # Number of assets
  b = t(chol(corr))             # Cholesky decomposition of Correlation matrix   
  Y = matrix(rnorm(M_sim*dim, 0, diag(dim)), 
             nrow = M_sim, byrow = TRUE)    # generate Y ~ N(0,I)
  Z = Y %*% b    # ~ N(0,p) where p is the correlation matrix
  
  spot.price = InAssVal * exp(matrix(rep(r - 0.5*(sigma^2)*T, M_sim), nrow = M_sim, 
                                     byrow = TRUE) 
                              + matrix(rep(sigma, M_sim), nrow = M_sim, 
                                       byrow = TRUE)*sqrt(T)*Z) 
  fz = exp(-r*T) * pmax(rowSums(spot.price)/3 - K, 0)  # discounted payoff
  
  option.price = mean(fz)
  variance = var(fz)
  
  return(c(option.price, variance))
}

Naive.MonteCarlo = BasketPrice(S0, r, sigma, corr, M, K, T)

## Apply importance sampling technique
# Calculate Option Price at time T
library(mvtnorm)
OptionPrice = function(InAssVal, r, sigma, corr, M_sim, lambda, K, T) {
  # InAssVal: Initial Asset Values (S0)
  # r: Interest rate on bank account
  # corr: Correlation matrix
  
  dim = nrow(corr)              # Number of assets
  b = t(chol(corr))             # Cholesky decomposition of Correlation matrix
  B = lambda * b %*% t(sigma)   # Beta
  
  x = matrix(rnorm(M_sim*dim, B, diag(dim)), nrow = M_sim, 
             byrow = TRUE)  # generate x from h ~ N(B,I) where I is identity matrix
  
  spot.price = InAssVal * exp(matrix(rep(r - 0.5*(sigma^2)*T, M_sim), 
                                     nrow = M_sim,byrow = TRUE)+ 
                                matrix(rep(sigma, M_sim), nrow = M_sim, 
                                       byrow = TRUE)*sqrt(T)*x) 
  fx = exp(-r*T) * pmax(rowSums(spot.price)/3 - K, 0)  # discounted payoff
  gx = dmvnorm(x, rep(0,dim), corr)      # pdf of g ~ MVN(0,p) where p is correlation matrix
  hx = dmvnorm(x, B, diag(dim))          # importance sampling distribution N(B,I)
  
  wx = gx/hx   # importance weight
  option.price = mean(fx * wx)
  variance = var(fx * wx)
  return(c(option.price, variance))
}

# Calculate lambda using uniroot(searching the interval for a root)
L = function(lambda){
  lambda- sqrt(T) * exp(lambda * sqrt(T)* sigma %*% corr %*% t(sigma))/
    exp(lambda * sqrt(T) * sigma %*% corr %*% t(sigma) - K/(sum(S0)/3))
}

lowerbound = log(K/(sum(S0)/3))/sqrt(T)* sigma %*% corr %*% t(sigma)
upperbound = 3  

lambda = uniroot(L, c(lowerbound, upperbound))
lambdastar = lambda$root

Importance.Sampling = OptionPrice(S0, r, sigma, corr, M, lambdastar, K, T)

### Compare the results
result = rbind(Naive.MonteCarlo, Importance.Sampling) 
colnames(result) = c("Price", "Variance")
result
