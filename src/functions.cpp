//[[Rcpp::depends(RcppEigen)]]
#include <math.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <string.h>
#include "functions.h"

using namespace std;
using namespace Eigen;
using namespace Rcpp;
using Eigen::MatrixXd;

//'@title posteriorX
//'
//'@description
//'In ordinary regression model, assuming each column of the design matrix \code{Z} contains additive
//'Gaussian measurement errors with known variance collected in \code{tauList}, function 'posteriorX'
//'returns posterior expectation as well as some related values for the true design matrix \eqn{X}.
//'
//'More specifically, if we denote the measurement error matrix as \eqn{\Xi},
//'then the measurement error model takes the form: \eqn{Z = X + \Xi}.
//'Since we have \eqn{y_i|\bold{x_i},\bold{\beta},\sigma^2 \sim N(\bold{x_i^\top }\bold{\beta},\sigma^2)}
//'and \eqn{\bold{x_i}|\bold{z_i}, \Lambda \sim N(\bold{z_i}, \Lambda)}
//'where \eqn{\Lambda} is the diagonal matrix consists of \code{tauList},
//'\code{posteriorX} computes the expectation of \eqn{X}, variance of \eqn{\bold{x_i}} and
//'expectation of \eqn{X^\top X} conditioning on \eqn{\bold{\beta}}, \eqn{\bold{y}}, \eqn{Z},
//'\eqn{\Lambda} and \eqn{\sigma}.
//'
//'This is part of E-step in the 'lsum' function.
//'
//'
//'
//'@param Z \eqn{n * p} covariate matrix with additive errors
//'in each column (possibly \eqn{p>n}). It is assumed the
//'measurement error model takes the form \eqn{Z = X + \Xi},
//'where \eqn{X} is the unknown true design matrix and
//'\eqn{\Xi} is the matrix of i.i.d Gaussian measurement errors of mean \eqn{0}
//'and same variance within each column.
//'
//'@param y Numeric response vector from \eqn{n} observations.
//'
//'@param beta Regression parameter vector of length \eqn{p}.
//'
//'@param sigma Numeric. The standard deviance of regression model error terms.
//'
//'@param tauList Vector of length \eqn{p}. Corresponding variances of additive
//'measurement error for p columns.
//'
//'
//'@returns
//''posteriorX' returns a list containing the following values:
//'
//'\item{hatX}{\eqn{n * p} matrix. The expectation of \eqn{X} conditional on
//'\eqn{\bold{\beta}}, \eqn{\bold{y}}, \eqn{Z}, \eqn{\Lambda} and \eqn{\sigma}.}
//'
//'\item{hatSigma}{\eqn{p * p} matrix. The variance of \eqn{\bold{x_i}} conditional on
//'\eqn{\bold{\beta}}, \eqn{\bold{y}}, \eqn{Z}, \eqn{\Lambda} and \eqn{\sigma},
//'which is equivalent for each \eqn{\bold{x_i}}.}
//'
//'\item{hatXX}{\eqn{p * p} matrix. The conditional expectation of \eqn{X^\top X}.
//'\code{hatXX = t(hatX) * hatX + hatSigma}.}
//'
//'@examples
//'require(MASS)
//'n = 200
//'p = 500
//'
//'set.seed(1234)
//'beta_true = c(-3, 2, -1.5, -2, 3, rep(0,p-5))
//'X = matrix(rnorm(n*p, 0, 1), nrow = n)
//'epsilon = rnorm(n, 0, 1)
//'y = X %*% beta_true+epsilon
//'
//'tau = sample(seq(0.5,0.9,by = 0.1), size = p, replace = TRUE)
//'Xi = mvrnorm(n, rep(0,p), diag(tau))
//'Z = X + Xi
//'
//'post_X = posteriorX(Z, y, beta = beta_true, sigma = 1, tauList = tau)
//'
//'
//'@useDynLib LS3MU
//'@import Rcpp
//'@import RcppEigen
//'@export
// [[Rcpp::export]]
List posteriorX(const Eigen::MatrixXd Z, const Eigen::VectorXd y,
                const Eigen::VectorXd beta, const float sigma,
                const Eigen::VectorXd tauList) {
  // Z: n*p matrix.
  // tauList: List of length p. Variance of additive measurement error.
  int p = Z.cols();
  int n = Z.rows();
  
  Eigen::MatrixXd HatSigma(p,p), HatX(n,p), HatXX(p,p), tauMat(p,p);
  
  tauMat = tauList.asDiagonal();
  HatSigma = beta * beta.transpose() / pow(sigma,2) + tauMat.inverse();
  HatSigma = HatSigma.inverse();
  
  // E(X|.): n*p matrix where each row is E(x_i|.)
  HatX = (y * beta.transpose() / pow(sigma,2) + Z * tauMat.inverse()) * HatSigma;
  // E(X^T X|.):summation of E(x_i x_i^T|.)
  HatXX = n * HatSigma + HatX.transpose() * HatX;
  
  return List::create(
    _["hatX"] = HatX,
    _["hatXX"] = HatXX,
    _["hatSigma"] = HatSigma
  );
}


//'@useDynLib LS3MU
//'@import Rcpp
//'@import RcppEigen
//'@export
// [[Rcpp::export]]
float comp_g(const bool Laplacian, const Eigen::VectorXd y, const Eigen::MatrixXd Z,
             const Eigen::VectorXd tauList, const Eigen::MatrixXd hatSigma,
             const Eigen::MatrixXd hatX, const Eigen::MatrixXd hatXX,
             const NumericVector p_k, const float spike_param,
             const float slab_param, const float omega,
             const float kappa, const float a, const float b, const float sigma,
             const float theta, const Eigen::VectorXd beta){
  int p = Z.cols();
  int n = Z.rows();
  //int d = tauList.size();
  
  float g, Const1, Const2, Const3, Const_LG, obj_func_wo_pen, pen;
  
  MatrixXd tauMat = tauList.asDiagonal();
  
  //Const1: sum E[log f(x_i|z_i,Lambda)]
  Rcpp::NumericVector taulist (wrap(tauList));
  Const1 = -n*p/2*log(2*M_PI) -n/2* sum(log(taulist)) -0.5*( tauMat.inverse()*(n*hatSigma+(Z-hatX).transpose()*(Z-hatX)) ).trace();
  Const2 = -n*p/2*log(2*M_PI) - n/2* log(hatSigma.determinant()) - n*p/2;
  Const3 = sum(p_k * log(p_k)) + sum((1-p_k)*log(1-p_k));
  obj_func_wo_pen = -1/2/pow(sigma,2)*((y.transpose() * y)(0)- 2*(y.transpose()*hatX*beta)(0)
                                         +(beta.transpose()* hatXX* beta)(0)) -(n+omega+2)*log(sigma)-omega*kappa/2/pow(sigma,2)+(p-sum(p_k)+b-1)* log(1-theta)+(a-1+sum(p_k))*log(theta);
  
  Rcpp::NumericVector beta_NV(wrap(beta));
  NumericVector pen_params;
  pen_params = 1/spike_param*(1-p_k) + 1/slab_param *p_k;
  if (Laplacian==true){
    Const_LG = -log(2)*p - (p-sum(p_k))*log(spike_param) - sum(p_k)*log(slab_param);
    pen = -sum(pen_params*abs(beta_NV));
  }else{
    Const_LG = -p/2*log(2*M_PI) - 0.5*(p-sum(p_k))*log(spike_param) - 0.5*sum(p_k)*log(slab_param);
    pen = - 0.5 * sum(pen_params*pow(beta_NV,2));
  }
  g =-n/2*log(2*M_PI)+ Const1-Const2 - Const3 + Const_LG +obj_func_wo_pen + pen;
  return g;
  
}