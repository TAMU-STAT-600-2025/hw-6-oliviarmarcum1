# [ToDo] Standardize X and Y: center both X and Y; scale centered X
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y){
  # [ToDo] Center Y
  Ymean <- mean(Y)
  Ytilde <- as.numeric(Y - Ymean)
  
  # [ToDo] Center and scale X
  n <- nrow(X)
  Xmeans <- colMeans(X)
  Xc <- sweep(X, 2, Xmeans, FUN = "-")
  # weights = sqrt( (X_j^T X_j)/n ) after centering, before scaling
  weights <- sqrt(colSums(Xc^2) / n)
  weights[weights == 0] <- 1  # safe-guard for constant columns
  Xtilde <- sweep(Xc, 2, weights, FUN = "/")
  
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

# [ToDo] Soft-thresholding of a scalar a at level lambda 
# [OK to have vector version as long as works correctly on scalar; will only test on scalars]
soft <- function(a, lambda){
  sign(a) * pmax(abs(a) - lambda, 0)
}

# [ToDo] Calculate objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda){
  n <- nrow(Xtilde)
  r <- as.numeric(Ytilde - Xtilde %*% beta)
  (sum(r^2) / (2*n)) + lambda * sum(abs(beta))
}

# [ToDo] Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1 (vector)
# lamdba - tuning parameter
# beta_start - p vector, an optional starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.001){
  #[ToDo]  Check that n is the same between Xtilde and Ytilde
  if (nrow(Xtilde) != length(Ytilde)) stop("Xtilde and Ytilde have incompatible n.")
  
  #[ToDo]  Check that lambda is non-negative
  if (lambda < 0) stop("lambda must be non-negative.")
  
  #[ToDo]  Check for starting point beta_start. 
  # If none supplied, initialize with a vector of zeros.
  # If supplied, check for compatibility with Xtilde in terms of p
  p <- ncol(Xtilde)
  beta <- if (is.null(beta_start)) rep(0, p) else {
    if (length(beta_start) != p) stop("beta_start length must equal ncol(Xtilde).")
    as.numeric(beta_start)
  }
  
  #[ToDo]  Coordinate-descent implementation. 
  # Stop when the difference between objective functions is less than eps for the first time.
  # For example, if you have 3 iterations with objectives 3, 1, 0.99999,
  # your should return fmin = 0.99999, and not have another iteration
  n <- nrow(Xtilde)
  r <- as.numeric(Ytilde - Xtilde %*% beta)  # residuals
  prev_obj <- Inf
  maxit <- 10000
  for (iter in 1:maxit) {
    # one full sweep over coordinates
    for (j in 1:p) {
      xj <- Xtilde[, j]
      r <- r + xj * beta[j]                        # add back old contribution
      rho <- sum(xj * r) / n                       # (1/n) * x_j^T r
      beta[j] <- soft(rho, lambda)                 # unit scaling => denom = 1
      r <- r - xj * beta[j]                        # subtract new contribution
    }
    obj <- (sum(r^2) / (2*n)) + lambda * sum(abs(beta))
    if (abs(prev_obj - obj) < eps) {
      return(list(beta = beta, fmin = obj))
    }
    prev_obj <- obj
  }
  list(beta = beta, fmin = (sum(r^2) / (2*n)) + lambda * sum(abs(beta)))
}
  # Return 
  # beta - the solution (a vector)
  # fmin - optimal function value (value of objective at beta, scalar)


# [ToDo] Fit LASSO on standardized data for a sequence of lambda values. Sequential version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence,
#             is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized_seq <- function(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Check that n is the same between Xtilde and Ytilde
  if (nrow(Xtilde) != length(Ytilde)) stop("Xtilde and Ytilde have incompatible n.")
  
  # [ToDo] Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >= 0,
  # and make sure the values are sorted from largest to smallest.
  # If none of the supplied values satisfy the requirement,
  # print the warning message and proceed as if the values were not supplied.
  if (!is.null(lambda_seq)) {
    lambda_seq <- sort(unique(lambda_seq[lambda_seq >= 0]), decreasing = TRUE)
    if (length(lambda_seq) == 0) {
      warning("No non-negative values in supplied lambda_seq; computing automatically.")
      lambda_seq <- NULL
    }
  }
  
  
  # If lambda_seq is not supplied, calculate lambda_max 
  # (the minimal value of lambda that gives zero solution),
  # and create a sequence of length n_lambda as
  if (is.null(lambda_seq)) {
    n <- nrow(Xtilde)
    lam_max <- max(abs(as.numeric(crossprod(Xtilde, Ytilde)) / n))
    lam_min <- lam_max * 1e-3                 # ratio, not absolute 0.01
    lambda_seq <- exp(seq(log(lam_max), log(lam_min), length.out = n_lambda))
  }
  
  # [ToDo] Apply fitLASSOstandardized going from largest to smallest lambda 
  # (make sure supplied eps is carried over). 
  # Use warm starts strategy discussed in class for setting the starting values.
  p <- ncol(Xtilde)
  L <- length(lambda_seq)
  beta_mat <- matrix(0, nrow = p, ncol = L)
  fmin_vec <- numeric(L)
  beta_start <- rep(0, p)
  for (i in 1:L) {
    fit <- fitLASSOstandardized(Xtilde, Ytilde, lambda_seq[i], beta_start = beta_start, eps = eps)
    beta_mat[, i] <- fit$beta
    fmin_vec[i] <- fit$fmin
    beta_start <- fit$beta  # warm start
  }
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec)
}

# [ToDo] Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Center and standardize X,Y based on standardizeXY function
  stan <- standardizeXY(X, Y)
  Xtilde <- stan$Xtilde; Ytilde <- stan$Ytilde
  
  # [ToDo] Fit Lasso on a sequence of values using fitLASSOstandardized_seq
  # (make sure the parameters carry over)
  path <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
  lambda_seq <- path$lambda_seq
  
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector
  # for each lambda
  p <- ncol(X)
  L <- length(lambda_seq)
  beta_mat <- matrix(0, nrow = p, ncol = L)
  beta0_vec <- numeric(L)
  for (i in seq_len(L)) {
    beta_orig  <- path$beta_mat[, i] / stan$weights
    beta0_orig <- as.numeric(stan$Ymean - sum(stan$Xmeans * beta_orig))
    beta_mat[, i]  <- beta_orig
    beta0_vec[i]   <- beta0_orig
  }
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec)
}


# [ToDo] Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# fold_ids - (optional) vector of length n specifying the folds assignment (from 1 to max(folds_ids)), if supplied the value of k is ignored 
# eps - precision level for convergence assessment, default 0.001
cvLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, k = 5, fold_ids = NULL, eps = 0.001){
  # [ToDo] Fit Lasso on original data using fitLASSO
  whole <- fitLASSO(X, Y, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
  if (is.null(whole$beta_mat) || is.null(whole$beta0_vec))
    stop("fitLASSO did not return beta_mat/beta0_vec; check fitLASSO implementation.")
  
  lambda_seq <- whole$lambda_seq
  L <- length(lambda_seq)
  n <- nrow(X)
  p <- ncol(X)
  
  # [ToDo] If fold_ids is NULL, split the data randomly into k folds.
  # If fold_ids is not NULL, split the data according to supplied fold_ids.
  if (is.null(fold_ids)) {
    perm <- sample.int(n)
    sizes <- rep(floor(n / k), k)
    if ((n %% k) > 0) sizes[seq_len(n %% k)] <- sizes[seq_len(n %% k)] + 1
    starts <- cumsum(c(1, head(sizes, -1)))
    ends   <- cumsum(sizes)
    fold_ids <- integer(n)
    for (j in 1:k) fold_ids[perm[starts[j]:ends[j]]] <- j
  } else {
    if (length(fold_ids) != n) stop("fold_ids must have length n.")
    k <- max(fold_ids)
  }
  # [ToDo] Calculate LASSO on each fold using fitLASSO,
  # and perform any additional calculations needed for CV(lambda) and SE_CV(lambda)
  mse <- matrix(NA_real_, L, k)
  for (fold in 1:k) {
    idx_te <- which(fold_ids == fold)
    idx_tr <- setdiff(seq_len(n), idx_te)
    
    fit_tr <- fitLASSO(X[idx_tr, , drop = FALSE], Y[idx_tr],
                       lambda_seq = lambda_seq,  # reuse SAME grid
                       n_lambda   = length(lambda_seq),
                       eps        = eps)
    
    if (is.null(fit_tr$beta_mat) || is.null(fit_tr$beta0_vec))
      stop("fitLASSO (fold) did not return beta_mat/beta0_vec; check fitLASSO implementation.")
    
    # vectorized predictions for speed
    Xte <- X[idx_te, , drop = FALSE]
    Yte <- Y[idx_te]
    Yhat_mat <- matrix(fit_tr$beta0_vec, nrow = length(Yte), ncol = L, byrow = TRUE) +
      Xte %*% fit_tr$beta_mat
    mse[, fold] <- colMeans((Yhat_mat - Yte)^2)
  }
  cvm <- rowMeans(mse)
  cvse <- apply(mse, 1, function(z) sd(z) / sqrt(k))
  
  # [ToDo] Find lambda_min
  i_min <- which.min(cvm)
  lambda_min <- lambda_seq[i_min]
  
  # [ToDo] Find lambda_1SE
  thresh <- cvm[i_min] + cvse[i_min]
  i_1se <- max(which(cvm <= thresh))
  lambda_1se <- lambda_seq[i_1se]
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # fold_ids - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  list(
    lambda_seq = lambda_seq,
    beta_mat   = whole$beta_mat,   # <-- defined here from 'whole'
    beta0_vec  = whole$beta0_vec,  # <-- defined here from 'whole'
    fold_ids   = fold_ids,
    lambda_min = lambda_min,
    lambda_1se = lambda_1se,
    cvm        = cvm,
    cvse       = cvse
  )
}

