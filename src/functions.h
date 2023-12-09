#ifndef functions_H
#define functions_H

//[[Rcpp::depends(RcppEigen)]]
#include <math.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <string.h>
using namespace std;
using namespace Eigen;
using namespace Rcpp;


Rcpp::List posteriorX(const Eigen::MatrixXd Z, const Eigen::VectorXd y,
                      const Eigen::VectorXd beta, const float sigma,
                      const Eigen::VectorXd TauList);

float comp_g(const bool Laplacian, const Eigen::VectorXd y, const Eigen::MatrixXd Z,
             const Eigen::VectorXd tauList, const Eigen::MatrixXd hatSigma,
             const Eigen::MatrixXd hatX, const Eigen::MatrixXd hatXX,
             const NumericVector p_k, const float spike_param,
             const float slab_param, const float omega,
             const float kappa, const float a, const float b, const float sigma,
             const float theta, const Eigen::VectorXd beta);

#endif