# constant for studentization
tau <- function(x, y) {
  xbar <- mean(x)
  ybar <- mean(y)
  mu_20_vec <- (x - xbar)^2
  mu_02_vec <- (y - ybar)^2
  mu_22_vec <- mu_20_vec * mu_02_vec
  mu_20 <- mean(mu_20_vec)
  mu_02 <- mean(mu_02_vec)
  mu_22 <- mean(mu_22_vec)
  tau_hat <- sqrt(mu_22/(mu_20*mu_02))
  return(tau_hat)
}

# studentized pearson correlation
rho_pearson_s_func <- function(x, y) {
  rho_c <- cor(x, y, method = "pearson")
  tau_hat <- tau(x, y)
  rho_s <- rho_c/tau_hat
  return(rho_s)
}

# CCC
rho_ccc_func <- function(x, y) {
  xbar=mean(x)
  ybar=mean(y)
  Sigma <- cov(data.frame(x, y), method = "pearson") 
  rho_c <- 2*Sigma[1,2]/(Sigma[1,1] + Sigma[2,2] + (xbar - ybar)^2)
  return(rho_c)
}

# studentized CCC
rho_ccc_s_func <- function(x, y) {
  rho_c <- rho_ccc_func(x, y)
  tau_hat <- tau(x, y)
  rho_s <- rho_c/tau_hat
  return(rho_s)
}

