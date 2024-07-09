library(DescTools)
# https://github.com/hyu-ub/perk
# library(perk)
library(tidyverse)
library(furrr)
dmg_fwd_CCC <- function(df, smp, ci = "z-transform", nperm = NULL, nproc = 1) {
  # Register a parallel backend to use
  options(future.globals.maxSize = 2291289600)
  plan(multisession, workers = nproc) # Choose the number of workers based on your CPU
  
  # Parallelized code
  results <- future_map_dfr(smp, function(X) {
    dist_data <- df %>%
      filter(label == X) %>%
      select(tax_name, label, starts_with("fwf")) %>%
      pivot_longer(names_to = "type", values_to = "f_fwd", c(-tax_name, -label)) %>%
      mutate(x = gsub("fwf", "", type)) %>%
      select(-type) %>%
      mutate(x = as.numeric(x))
    dist_fit <- df %>%
      filter(label == X) %>%
      select(tax_name, label, matches("^fwdx\\d+")) %>%
      pivot_longer(names_to = "type", values_to = "Dx_fwd", c(-tax_name, -label)) %>%
      mutate(x = gsub("fwdx", "", type)) %>%
      select(-type) %>%
      mutate(x = as.numeric(x))
    data <- dist_data %>%
      inner_join(dist_fit, by = join_by(tax_name, label, x)) %>%
      group_by(tax_name) %>%
      arrange(x, .by_group = TRUE)
    if (nrow(dist_data) == 0) {
      return(NULL)
    } else {
      fits <- data %>%
        do({
          ccc <- DescTools::CCC(.$f_fwd, .$Dx_fwd, ci = ci)
          if (!is.null(nperm)) {
            ccc_perm <- perk_test(.$f_fwd, .$Dx_fwd, B = nperm, method = "ccc", alternative = "greater")
            tibble(rho_c = ccc$rho.c$est, rho_lwr_ci = ccc$rho.c$lwr.ci, rho_upr_ci = ccc$rho.c$upr.ci, C_b = ccc$C.b, l_shift = ccc$l.shift, s_shift = ccc$s.shift, label = X, rho_c_perm = ccc_perm$estimate, rho_c_perm_pval = ccc_perm$p.value)
          } else {
            tibble(rho_c = ccc$rho.c$est, rho_lwr_ci = ccc$rho.c$lwr.ci, rho_upr_ci = ccc$rho.c$upr.ci, C_b = ccc$C.b, l_shift = ccc$l.shift, s_shift = ccc$s.shift, label = X)
          }
        })
    }
  }, .progress = TRUE) # Optional progress bar
  
  # Reset the default plan if needed
  plan(sequential)
  return(results)
}
