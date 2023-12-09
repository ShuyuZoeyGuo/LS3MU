// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Laplacian_path
List Laplacian_path(const Eigen::MatrixXd Z, const Eigen::VectorXd y, const Eigen::VectorXd tauList, const NumericVector spike_params, const float slab_param, const Eigen::VectorXd beta_init, const bool sigma_update, const float sigma_init, const float theta_init, const float a, const float b, const float omega, const float kappa, const float tolerance, const int max_iter, const bool return_g);
RcppExport SEXP _LS3MU_Laplacian_path(SEXP ZSEXP, SEXP ySEXP, SEXP tauListSEXP, SEXP spike_paramsSEXP, SEXP slab_paramSEXP, SEXP beta_initSEXP, SEXP sigma_updateSEXP, SEXP sigma_initSEXP, SEXP theta_initSEXP, SEXP aSEXP, SEXP bSEXP, SEXP omegaSEXP, SEXP kappaSEXP, SEXP toleranceSEXP, SEXP max_iterSEXP, SEXP return_gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type tauList(tauListSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type spike_params(spike_paramsSEXP);
    Rcpp::traits::input_parameter< const float >::type slab_param(slab_paramSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< const bool >::type sigma_update(sigma_updateSEXP);
    Rcpp::traits::input_parameter< const float >::type sigma_init(sigma_initSEXP);
    Rcpp::traits::input_parameter< const float >::type theta_init(theta_initSEXP);
    Rcpp::traits::input_parameter< const float >::type a(aSEXP);
    Rcpp::traits::input_parameter< const float >::type b(bSEXP);
    Rcpp::traits::input_parameter< const float >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const float >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const float >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< const int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const bool >::type return_g(return_gSEXP);
    rcpp_result_gen = Rcpp::wrap(Laplacian_path(Z, y, tauList, spike_params, slab_param, beta_init, sigma_update, sigma_init, theta_init, a, b, omega, kappa, tolerance, max_iter, return_g));
    return rcpp_result_gen;
END_RCPP
}
// posteriorX
List posteriorX(const Eigen::MatrixXd Z, const Eigen::VectorXd y, const Eigen::VectorXd beta, const float sigma, const Eigen::VectorXd tauList);
RcppExport SEXP _LS3MU_posteriorX(SEXP ZSEXP, SEXP ySEXP, SEXP betaSEXP, SEXP sigmaSEXP, SEXP tauListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const float >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type tauList(tauListSEXP);
    rcpp_result_gen = Rcpp::wrap(posteriorX(Z, y, beta, sigma, tauList));
    return rcpp_result_gen;
END_RCPP
}
// comp_g
float comp_g(const bool Laplacian, const Eigen::VectorXd y, const Eigen::MatrixXd Z, const Eigen::VectorXd tauList, const Eigen::MatrixXd hatSigma, const Eigen::MatrixXd hatX, const Eigen::MatrixXd hatXX, const NumericVector p_k, const float spike_param, const float slab_param, const float omega, const float kappa, const float a, const float b, const float sigma, const float theta, const Eigen::VectorXd beta);
RcppExport SEXP _LS3MU_comp_g(SEXP LaplacianSEXP, SEXP ySEXP, SEXP ZSEXP, SEXP tauListSEXP, SEXP hatSigmaSEXP, SEXP hatXSEXP, SEXP hatXXSEXP, SEXP p_kSEXP, SEXP spike_paramSEXP, SEXP slab_paramSEXP, SEXP omegaSEXP, SEXP kappaSEXP, SEXP aSEXP, SEXP bSEXP, SEXP sigmaSEXP, SEXP thetaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const bool >::type Laplacian(LaplacianSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type tauList(tauListSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type hatSigma(hatSigmaSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type hatX(hatXSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type hatXX(hatXXSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type p_k(p_kSEXP);
    Rcpp::traits::input_parameter< const float >::type spike_param(spike_paramSEXP);
    Rcpp::traits::input_parameter< const float >::type slab_param(slab_paramSEXP);
    Rcpp::traits::input_parameter< const float >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const float >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const float >::type a(aSEXP);
    Rcpp::traits::input_parameter< const float >::type b(bSEXP);
    Rcpp::traits::input_parameter< const float >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const float >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(comp_g(Laplacian, y, Z, tauList, hatSigma, hatX, hatXX, p_k, spike_param, slab_param, omega, kappa, a, b, sigma, theta, beta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LS3MU_Laplacian_path", (DL_FUNC) &_LS3MU_Laplacian_path, 16},
    {"_LS3MU_posteriorX", (DL_FUNC) &_LS3MU_posteriorX, 5},
    {"_LS3MU_comp_g", (DL_FUNC) &_LS3MU_comp_g, 17},
    {NULL, NULL, 0}
};

RcppExport void R_init_LS3MU(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}