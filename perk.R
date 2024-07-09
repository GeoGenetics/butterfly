#' Robust permutation tests of correlation coefficients
#' Studentized permutation test of Lin's concordance correlation coefficient.
#' @details Robust permutation test of Pearson correlation coefficient, Spearman correlation coefficient, and Lin's concordance correlation coefficient.
#' The test robustly controls type I error under general scenarios
#' when sample size is small and two variables are dependent but uncorrelated.
#' Currently, for Pearson and Spearman correlations, the test only supports a zero null hypothesis, but the alternative hypothesis can be 
#' either one-sided or two-sided. For CCC, the test supports a more general null hypothesis. Currently, the package only supports a one-sided 
#' alternative hypothesis (greater).
#' @param x a \code{numeric} vector
#' @param y a \code{numeric} vector
#' @param B an \code{interger} number of permutations
#' @param r0 a \code{numeric} denoting the CCC under the null hypothesis. It should be in the range between -1 and 1. 
#' This parameter will be ignored for the tests of Pearson or Spearman's correlation coefficient.
#' @param method the correlation coefficient to be tested, options include Pearson's correlation coefficient (\code{pearson}), 
#' Spearman's correlation coefficient (\code{spearman}), and Lin's concordance correlation coefficient (\code{ccc}).
#' @param alternative the alternative hypothesis, can be \code{two.sided}, \code{less}, or \code{greater}
#' @return 
#' \describe{
#' \item{\code{estimate}}{the estimated correlation coefficient.}
#' \item{\code{p.value}}{the p-value from the studentized permutation test.}
#' \item{\code{method}}{the method for measuring correlation coefficient.}
#' \item{\code{alternative}}{the alternative hypothesis.}
#' }
#' 
#' @author Han Yu, Alan Hutson
#' @references Lawrence, I., & Lin, K. (1989). A concordance correlation coefficient to evaluate reproducibility. Biometrics, 255-268. \cr
#' \cr
#' DiCiccio, C. J., & Romano, J. P. (2017). Robust permutation tests for correlation and regression coefficients. Journal of the American Statistical Association, 112(519), 1211-1220.
#' @importFrom stats cov cor var
#' @examples
#' set.seed(123)
#' x <- rnorm(20)
#' y <- rnorm(20)
#' perk_test(x, y, B = 500, method = "pearson", alternative = "greater")
#' @export

perk_test <- function(x, y, B = 1000, r0=0, method = c("pearson", "spearman", "ccc"), 
                      alternative = c("two.sided", "less", "greater")) {
  if (length(method)> 1) method <- method[1]
  if (!method %in% c("pearson", "spearman", "ccc")) {
    warning("Method not defined, using pearson correlation coefficient as default.")
    method <- "pearson"
  }
  if (!alternative %in% c("two.sided", "less", "greater")) {
    warning("Alternative hypothesis not defined, using default two-sided test.")
    alternative <- "two.sided"
  }
  
  if (length(alternative)> 1) alternative <- alternative[1]
  if (method == "pearson" || method == "spearman") {
    if (r0 != 0) {
      warning("Non-zero null hypothesis is not supported for Pearson/Spearman correlation, using default H0: r0 = 0.")
    }
    res <- perm_cor(x, y, B, method, alternative = alternative) 
  } else if (method == "ccc") {
    if (r0 == 0) {
      res <- perm_ccc(x, y, B) 
    } else if (r0 >= -1 && r0 <= 1) {
      res <- test_permute_stu_2(x, y, rc0=r0, B=500)
    } else {
      stop("r0 must be in the rage between -1 and 1.")
    }
    
    if (alternative != "greater") {
      warning("Only one-sided (greater) hypothesis is supported for CCC, setting alternative as 'greater'.")
    }
    
  }
  return(res)
}
