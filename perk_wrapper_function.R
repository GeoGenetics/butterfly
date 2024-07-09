# permutation test for pearson or spearman correlation coefficients
perm_cor <- function(x, y, B, method = c("pearson", "spearman"), alternative = c("two.sided", "less", "greater")) {
  n <- length(x)
  if (method == "spearman") {
    x <- rank(x)
    y <- rank(y)
  }
  rho_c <- cor(x, y, method = "pearson")
  rho_s <- rho_pearson_s_func(x, y) # studentized correlation
  rho_s_star <- replicate(B, rho_pearson_s_func(sample(x), y), simplify = TRUE)
  if (alternative == "two.sided") {
    p_value <- mean(rho_s_star > abs(rho_s) | rho_s_star < -abs(rho_s))
  } else if (alternative == "less") {
    p_value <- mean(rho_s_star < abs(rho_s))
  } else if (alternative == "greater") {
    p_value <- mean(rho_s_star > rho_s)
  }
  res <- list(estimate = rho_c, p.value = p_value, method = method, alternative = alternative)
  return(res)
}

# permutation test for ccc
perm_ccc <- function(x, y, B) {
  n <- length(x)
  rho_c <- rho_ccc_func(x, y)
  rho_s <- rho_ccc_s_func(x, y)
  rho_s_star <- replicate(B, rho_ccc_s_func(sample(x), y), simplify = TRUE)
  p_value <- mean(rho_s_star > rho_s)
  res <- list(estimate = rho_c, p.value = p_value, method = "ccc", alternative = "greater")
  return(res)
}

test_permute_stu_2 <- function(x, y, rc0, B=500) {
  
  n <- length(x)
  u <- (x - mean(x)) / sqrt(var(x)) 
  v <- (y - mean(y)) / sqrt(var(y))
  stat_sample <- stat_func(x, y, rc0)
  
  ###########
  rho_sample <- cor(u, v)
  C_0 <- 2*sqrt(var(x)*var(y))/(var(x) + var(y) + (mean(x) - mean(y))^2)
  rho_hat_0 <- rc0/C_0
  if (rho_hat_0 >= 1) return(res <- list(estimate = rho_ccc_func(x, y), p.value = 1, method = "ccc", alternative = "greater"))
  
  # de-correlate
  u_1 <- u
  v_1 <- (u - v/rho_sample)/sqrt(1 + 1/rho_sample^2)
  # re-correlate
  u_2 <- u_1
  v_2 <- rho_hat_0*u_1 + sqrt(1 - rho_hat_0^2)*v_1
  # re-scale and shift
  x_2 <- u_2*sqrt(var(x)) + mean(x)
  y_2 <- v_2*sqrt(var(y)) + mean(y) 
  #######
  
  rho_star_vec <- c()
  for (b in 1:n) {
    ind_x <- setdiff(1:n, b)
    x_star <- x_2[ind_x]
    y_star <- y_2[ind_x]
    rho_star_vec[b] <- stat_func(x_star, y_star, rc0)
  }
  
  tau_hat <- sqrt(var(rho_star_vec))
  stat_star_vec <- c()
  
  for (b in 1:B) {
    ind_x <- sample(1:n, n, replace = TRUE)
    ind_y <- sample(1:n, n, replace = TRUE)
    u_star <- u[ind_x]
    v_star <- v[ind_y]
    stat_star_vec[b] <- cor(u_star, v_star)
  }
  
  p_value <- mean(stat_star_vec*tau_hat*sqrt(n*(n-1)) > stat_sample)
  res <- list(estimate = rho_ccc_func(x, y), p.value = p_value, method = "ccc", alternative = "greater")
  return(res)
}


stat_func <- function(x, y, rc0) {
  C_hat <- 2*sqrt(var(x)*var(y))/(var(x) + var(y) + (mean(x) - mean(y))^2)
  rho_hat_0 <- rc0/C_hat
  u <- (x - mean(x)) / sqrt(var(x)) 
  v <- (y - mean(y)) / sqrt(var(y))
  v2 <- (v/rho_hat_0-u)/sqrt(1 + 1/rho_hat_0^2)
  rho_hat <- cor(u, v2)
}


