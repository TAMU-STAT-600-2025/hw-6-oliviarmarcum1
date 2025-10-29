
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

# Source your LASSO functions from HW4 (make sure to move the corresponding .R file in the current project folder)
source("LassoFunctions.R")

# Do at least 2 tests for soft-thresholding function below. You are checking output agreements on at least 2 separate inputs
#################################################
set.seed(1)
stopifnot(all.equal(soft_c(3.0, 0.5), soft(3.0, 0.5)))
stopifnot(all.equal(soft_c(-0.12, 0.2), soft(-0.12, 0.2)))
# (|a| < lambda and |a| > lambda)
stopifnot(all.equal(soft_c(0.05, 0.2), soft(0.05, 0.2)))
stopifnot(all.equal(soft_c(-2.3, 0.1), soft(-2.3, 0.1)))

# Do at least 2 tests for lasso objective function below. You are checking output agreements on at least 2 separate inputs
#################################################
n <- 120; p <- 50
Xraw <- matrix(rnorm(n * p), n, p)
beta_true <- c(rep(2, 5), rep(0, p - 5))
Yraw <- as.numeric(Xraw %*% beta_true + rnorm(n))

std <- standardizeXY(Xraw, Yraw)
Xtilde <- std$Xtilde
Ytilde <- std$Ytilde

lambda1 <- 0.15
lambda2 <- 0.6
b_try1 <- rnorm(p)
b_try2 <- rnorm(p)

stopifnot(all.equal(lasso_c(Xtilde, Ytilde, b_try1, lambda1),
                    lasso(Xtilde, Ytilde, b_try1, lambda1),
                    tolerance = 1e-10))
stopifnot(all.equal(lasso_c(Xtilde, Ytilde, b_try2, lambda2),
                    lasso(Xtilde, Ytilde, b_try2, lambda2),
                    tolerance = 1e-10))

# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################
lam_seq1 <- sort(exp(seq(log(1.2), log(0.05), length.out = 25)), decreasing = TRUE)
b0 <- rep(0, p)

fit_R_1 <- fitLASSOstandardized(Xtilde, Ytilde, lambda1, beta_start = b0, eps = 1e-8)
fit_C_1 <- fitLASSOstandardized_c(Xtilde, Ytilde, lambda1, beta_start = b0, eps = 1e-8)
stopifnot(all.equal(as.numeric(fit_R_1$beta), as.numeric(fit_C_1), tolerance = 1e-6))

# warm start from previous solution, different lambda
fit_R_2 <- fitLASSOstandardized(Xtilde, Ytilde, lambda2, beta_start = fit_R_1$beta, eps = 1e-8)
fit_C_2 <- fitLASSOstandardized_c(Xtilde, Ytilde, lambda2, beta_start = fit_R_1$beta, eps = 1e-8)
stopifnot(all.equal(as.numeric(fit_R_2$beta), as.numeric(fit_C_2), tolerance = 1e-6))


# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################
library(microbenchmark)
mb_path <- microbenchmark(
  R_path = fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq = lam_seq1, n_lambda = length(lam_seq1), eps = 1e-8),
  C_path = fitLASSOstandardized_seq_c(Xtilde, Ytilde, lambda_seq = lam_seq1, eps = 1e-8),
  times = 20
)
print(mb_path)


# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################

# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression

# Make sure riboflavin$x is treated as matrix later in the code for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

# Standardize the data
out <- standardizeXY(riboflavin$x, riboflavin$y)

# This is just to create lambda_seq, can be done faster, but this is simpler
outl <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)

# The code below should assess your speed improvement on riboflavin data
microbenchmark(
  fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq),
  fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq),
  times = 10
)
