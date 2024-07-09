calculate_plot_grid <- function(n) {
  cols <- ceiling(sqrt(n))
  rows <- ceiling(n / cols)
  
  list(rows = rows, cols = cols)
}
