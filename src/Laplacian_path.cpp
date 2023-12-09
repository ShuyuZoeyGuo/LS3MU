//[[Rcpp::depends(RcppEigen)]]
#include <math.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <string.h>
#include "functions.h"

using namespace std;
using namespace Eigen;
using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Rcpp::as;

//'@useDynLib LS3MU
//'@import Rcpp
//'@import RcppEigen
// [[Rcpp::export]]
List Laplacian_path(const Eigen::MatrixXd Z, const Eigen::VectorXd y, const Eigen::VectorXd tauList,
                    const NumericVector spike_params, const float slab_param, const Eigen::VectorXd beta_init,
                    const bool sigma_update, const float sigma_init, const float theta_init,
                    const float a, const float b, const float omega, const float kappa,
                    const float tolerance, const int max_iter, const bool return_g){
  
  VectorXd beta = beta_init;
  
  int n = Z.rows();
  int p = Z.cols();
  int L = spike_params.length();
  
  float sigma, sigma_sqr, theta, diff, spike_param, u_j,  g_k;
  MatrixXd hatX(n,p), hatXX(p,p), hatSigma(p,p), beta_path(L,p);
  NumericVector p_k, pen_params, g_l, iter_nums(L), sigma_path(L), theta_path(L);
  VectorXd beta_old;
  List g_List = List(L);
  
  for (int l = 0; l < L; l++){
    sigma = sigma_init;
    theta = theta_init;
    if (return_g == true){
      g_l = 0;
    }
    spike_param = spike_params[l];
    
    for (int k = 0; k < max_iter; k++){
      // E step for X
      // function posteriorX in functions.cpp
      Rcpp::List postX = posteriorX(Z, y, beta, sigma, tauList);
      hatX = postX["hatX"];
      hatXX = postX["hatXX"];
      hatSigma = postX["hatSigma"];
      
      // E step for gamma
      Rcpp::NumericVector beta_NV(wrap(beta));
      p_k = theta / (theta + (1-theta) * slab_param/spike_param * exp(-(1/spike_param-1/slab_param)*abs(beta_NV)));
      pen_params = (1-p_k)/spike_param + p_k/slab_param;
      
      ////  M-step  ////
      
      // coordinate-wisely update beta
      beta_old = beta;
      for (int j = 0; j < p; j++){
        u_j = (hatX.col(j).transpose() * y)(0) - (hatXX.row(j) * beta)(0) + hatXX(j,j)*beta(j);
        if (fabs(u_j) <= pow(sigma,2) * pen_params[j]){
          beta(j) = 0;
        }else if (u_j > 0){
          beta(j) = 1/hatXX(j,j) * (u_j - pow(sigma,2)*pen_params[j]);
        }else{
          beta(j) = -1/hatXX(j,j) * (fabs(u_j)-pow(sigma,2) * pen_params[j]);
        }
      }
      
      // update theta
      theta = (sum(p_k)+a-1) / (p+a+b-2);
      
      // update sigma
      if (sigma_update==true){
        sigma_sqr = 1/(n+omega+2)*(kappa*omega + (y.transpose() * y)(0) - 2* (y.transpose() * hatX * beta)(0)
                                     + (beta.transpose() * hatXX * beta)(0));
        sigma = pow(sigma_sqr,0.5);
      }
      
      // compute g
      // function comp_g in functions.cpp
      if (return_g == true){
        g_k = comp_g(true, y, Z, tauList, hatSigma, hatX, hatXX, p_k,
                     spike_param, slab_param, omega, kappa, a, b, sigma, theta, beta);
        g_l.push_back(g_k);
      }
      
      //break or not
      diff = pow(((beta-beta_old).transpose() * (beta-beta_old))(0), 0.5);
      if ((diff <= tolerance) || (k==max_iter-1)){
        beta_path.row(l) = beta;
        sigma_path(l) = sigma;
        theta_path(l) = theta;
        if (return_g == true){
          g_l.erase(g_l.begin());
          g_List(l) = g_l;
        }
        iter_nums(l) = k+1;
        break;
      }
    }
  }
  
  // obtain the final beta
  VectorXd beta_output = beta;
  
  Rcpp::NumericMatrix beta_path_(wrap(beta_path));
  Rcpp::NumericVector beta_output_(wrap(beta_output));
  if (return_g == true){
    return List::create(
      _["beta_path"] = beta_path_,
      _["beta_output"] = beta_output_,
      _["sigma_path"] = sigma_path,
      _["theta_path"] = theta_path,
      _["g_List"] = g_List,
      _["iter_nums"] = iter_nums,
      _["spike_params"] = spike_params,
      _["slab_param"] = slab_param);
  }else{
    return List::create(
      _["beta_path"] = beta_path_,
      _["beta_output"] = beta_output_,
      _["sigma_path"] = sigma_path,
      _["theta_path"] = theta_path,
      _["iter_nums"] = iter_nums,
      _["spike_params"] = spike_params,
      _["slab_param"] = slab_param);
  }
}