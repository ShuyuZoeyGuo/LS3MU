#' @title Laplacian Selector for Uncertain Matrix
#'
#' @description
#' For high dimensional sparse regression with additive errors in potential variables, this 
#' Laplacian Spike & Slab Selector provides a possible solution of variable selection 
#' under EM framework, where both unknown true variables and spike-slab indicators are treated as 
#' latent variables. The path of regression coefficients is obtained following a decreasing 
#' list of spike parameters. 
#'
#' @import stats
#' @import MASS
#'
#' @param Z \eqn{n * p} covariate matrix with additive error
#' in each column (possibly \eqn{p>n}). It is assumed the
#' measurement error model takes the form \eqn{Z = X + \Xi},
#' where \eqn{X} is the unknown true design matrix and
#' \eqn{\Xi} is the matrix of i.i.d Gaussian measurement error of mean \eqn{0}
#' and same variance within each column.
#'
#' @param y Numeric response vector from \eqn{n} observations.
#'
#' @param tauList Vector of length \eqn{p}. Corresponding variances of additive
#' measurement error for p columns.
#'
#' @param spike_params Vector of of length \eqn{L}. Decreasing scale parameters
#' of spike prior for \code{beta}.
#' \code{spike_params} should be less than \code{slab_param}. If not specified,
#' \code{spike_params} will be assigned default values. See 'Details' for
#' more information.
#'
#' @param slab_param Numeric. Scale parameter \eqn{\lambda_1} of slab prior for
#' \code{beta}. If not specified, \code{slab_param = spike_params[1] * 10}
#' if \code{spike_params} is given; Otherwise, \code{slab_param = 1} by default.
#'
#' @param beta_init Vector. Initial value of regression coefficients \code{beta}.
#' \code{beta_init = 0} by default.
#'
#' @param sigma_update Logical. Whether the variance of model error is updated or not.
#' Default is \code{TRUE}.
#'
#' @param sigma_init Numeric. The initial value of standard deviation of model error.
#' If not specified, \code{sigma} is initialized according to \code{sd(y)}.
#'
#' @param theta_init Numeric. The initial value of prior proportion of nonzero
#' `beta`. \code{theta_init} must be in \code{(0,1]}.
#' Default is \code{0.5}.
#'
#' @param return_g Logical. Default is \code{FALSE}. If specified \code{TRUE},
#' return a list containing the expected value of log likelihood function w.r.t
#' the current conditional distribution of latent variables.
#'
#' @param a,b Numeric parameters of Beta prior distribution of \code{theta}, where
#' \code{theta ~ Beta(a,b)}. \code{a = 1} and \code{b = p} where \code{p = ncol(Z)}
#' by default.
#'
#' @param omega,kappa Numeric parameters of Inverse Gamma prior distribution of
#' \code{sigma^2}, where \code{sigma^2 ~ IG(omega/2, omega*kappa/2)}.
#' \code{omega = 1} and \code{kappa = 1} by default.
#'
#' @param tolerance Numeric. Criterion for early stopping at each \code{spike_param}. If
#' \eqn{\left\| \beta_{old} - \beta_{new} \right\|_2 < \code{tolerance}},
#' then break the iteration at the current \code{spike_param}.
#'
#' @param max_iter Integer. The maximum iteration number at each \code{spike_param}.
#'
#'
#'
#' @returns
#' `lsum` returns a list containing the following values:
#'
#' \item{beta_path}{\eqn{L * p} matrix. Each row is the `beta` fitted at corresponding
#' \code{spike_param}.}
#'
#' \item{beta_indices}{Vector. Indices of selected nonzero regression parameters.}
#'
#' \item{beta_values}{Vector. Values of selected nonzero regression parameters.}
#'
#' \item{beta_output}{Vector of length \eqn{p}. Full output `beta` including zeros.}
#'
#' \item{sigma_path}{Vector of Length \eqn{L}. Estimated \code{sigma} at each
#' \code{spike_param}.}
#'
#' \item{theta_path}{Vector of Length \eqn{L}. Estimated \code{theta} at each
#' \code{spike_param}.}
#'
#' \item{g_List}{List of Length \eqn{L}. Each element of `g_List` is a list containing
#' values of maximized function \eqn{g} over iterations at corresponding
#' \code{spike_param}.}
#'
#' \item{iter_nums}{Vector of Length \eqn{L}. Number of Iterations at each
#' \code{spike_param}.}
#'
#'
#'
#' @details
#'
#' Since 'lsum' is built based on EM algorithm, it is possible to converge to local rather 
#' than global maximum values. Thus the result can be sensitive to the initial choices of
#' `beta`, `sigma` and `theta`.
#'
#' In addition, the value of spike and slab parameters can be crucial
#' to variable selection result in practice. We set some default values for
#' \code{spike_params} and \code{slab_param}.
#' If \code{spike_param} is not specified and \code{slab_param} is given,
#' \code{spike_params = exp(seq(log(slab_L),by=-0.3,length.out=20))[-1]} by default.
#' If `slab_param` is also not specified, the `slab_param = 1` by default and the 
#' `spike_param` is set the same as the previous way. Users can assign any valid values 
#' to \code{spike_params} and \code{slab_param}.
#'
#' @examples
#' require(MASS)
#' n = 200
#' p = 500
#'
#' set.seed(1234)
#' beta_true = c(-3, 2, -1.5, -2, 3, rep(0,p-5))
#' X = matrix(rnorm(n*p, 0, 1), nrow = n)
#' epsilon = rnorm(n, 0, 1)
#' y = X %*% beta_true+epsilon
#'
#' tau = sample(seq(0.5,0.9,by = 0.1), size = p, replace = TRUE)
#' Xi = mvrnorm(n, rep(0,p), diag(tau))
#' Z = X + Xi
#'
#' lsum_fit = lsum(Z, y, tauList = tau)
#'
#' @export
lsum = function(Z, y, tauList, spike_params, slab_param, beta_init, 
                sigma_update = TRUE, sigma_init, theta_init = 0.5, 
                return_g = FALSE, a = 1, b = NULL, omega = 1, kappa = 1, 
                tolerance = 0.01, max_iter = 500){
  
  n = nrow(Z)
  p = ncol(Z)
  
  if ((theta_init > 1) | (theta_init <= 0)){
    stop("'theta_init' should in (0,1]!")
  }
  
  if(missing(beta_init)){
    beta_init = rep(0,p)
  }
  
  
  if(missing(spike_params)){
    if(missing(slab_param)){
      slab_param = 1
    }
    spike_params = exp(seq(log(slab_param),by=-0.3,length.out=20))[-1]
  }else{
    if(identical(sort(spike_params,decreasing = TRUE),spike_params)==FALSE){
      warning("spike_params should be an monotone decreasing sequence.", immediate. = TRUE)
      spike_params = sort(spike_params, decreasing = TRUE)
    }
    if(spike_params[length(spike_params)]<=0){
      stop("Min spike_param should be positive!")
    }else{
      if(missing(slab_param)){
        slab_param = spike_params[1]*10
      }
      if(slab_param <= spike_params[1]){
        stop("Max spike_param should be less than slab_param!")
      }
    }
  }
  
  
  if(missing(sigma_init)){
    sigma_over = sd(y)
    df = 3
    sigma_init = sqrt(qchisq(0.99,df,lower.tail = FALSE)*sigma_over^2/(df+2))
  }
  
  if(is.null(b)){b=p}
  
  
  # call algorithms from cpp file
  lsum_fit = Laplacian_path(Z, y, tauList, spike_params, slab_param, beta_init,
                             sigma_update, sigma_init, theta_init, a, b,
                             omega, kappa, tolerance, max_iter, return_g)
  lsum_fit$beta_indices = which(lsum_fit$beta_output!=0)
  lsum_fit$beta_values = lsum_fit$beta_output[lsum_fit$beta_indices]
  return(lsum_fit)
}