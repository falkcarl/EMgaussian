// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// imp1matcov
arma::mat imp1matcov(Rcpp::NumericMatrix D, const arma::colvec& muest, const arma::mat& sigest);
RcppExport SEXP _EMgaussian_imp1matcov(SEXP DSEXP, SEXP muestSEXP, SEXP sigestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type muest(muestSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigest(sigestSEXP);
    rcpp_result_gen = Rcpp::wrap(imp1matcov(D, muest, sigest));
    return rcpp_result_gen;
END_RCPP
}
// imp2matcov
void imp2matcov(Rcpp::NumericMatrix D, const arma::mat& sigest, arma::mat& t2);
RcppExport SEXP _EMgaussian_imp2matcov(SEXP DSEXP, SEXP sigestSEXP, SEXP t2SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigest(sigestSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type t2(t2SEXP);
    imp2matcov(D, sigest, t2);
    return R_NilValue;
END_RCPP
}
// EMcyclecov
Rcpp::List EMcyclecov(const Rcpp::NumericMatrix& D, const arma::colvec& muest, const arma::mat& sigest);
RcppExport SEXP _EMgaussian_EMcyclecov(SEXP DSEXP, SEXP muestSEXP, SEXP sigestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type muest(muestSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigest(sigestSEXP);
    rcpp_result_gen = Rcpp::wrap(EMcyclecov(D, muest, sigest));
    return rcpp_result_gen;
END_RCPP
}
// nllcov
double nllcov(const arma::mat dat, const arma::colvec mu, const arma::mat sig);
RcppExport SEXP _EMgaussian_nllcov(SEXP datSEXP, SEXP muSEXP, SEXP sigSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type dat(datSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type sig(sigSEXP);
    rcpp_result_gen = Rcpp::wrap(nllcov(dat, mu, sig));
    return rcpp_result_gen;
END_RCPP
}
// imp1matprec
arma::mat imp1matprec(Rcpp::NumericMatrix D, const arma::colvec& muest, const arma::mat& kest);
RcppExport SEXP _EMgaussian_imp1matprec(SEXP DSEXP, SEXP muestSEXP, SEXP kestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type muest(muestSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type kest(kestSEXP);
    rcpp_result_gen = Rcpp::wrap(imp1matprec(D, muest, kest));
    return rcpp_result_gen;
END_RCPP
}
// imp2matprec
void imp2matprec(Rcpp::NumericMatrix D, const arma::mat& kest, arma::mat& t2);
RcppExport SEXP _EMgaussian_imp2matprec(SEXP DSEXP, SEXP kestSEXP, SEXP t2SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type kest(kestSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type t2(t2SEXP);
    imp2matprec(D, kest, t2);
    return R_NilValue;
END_RCPP
}
// EMcycleprec
Rcpp::List EMcycleprec(const Rcpp::NumericMatrix& D, const arma::colvec& muest, const arma::mat& kest);
RcppExport SEXP _EMgaussian_EMcycleprec(SEXP DSEXP, SEXP muestSEXP, SEXP kestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type muest(muestSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type kest(kestSEXP);
    rcpp_result_gen = Rcpp::wrap(EMcycleprec(D, muest, kest));
    return rcpp_result_gen;
END_RCPP
}
// nllprec
double nllprec(const arma::mat dat, const arma::colvec mu, const arma::mat K);
RcppExport SEXP _EMgaussian_nllprec(SEXP datSEXP, SEXP muSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type dat(datSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(nllprec(dat, mu, K));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EMgaussian_imp1matcov", (DL_FUNC) &_EMgaussian_imp1matcov, 3},
    {"_EMgaussian_imp2matcov", (DL_FUNC) &_EMgaussian_imp2matcov, 3},
    {"_EMgaussian_EMcyclecov", (DL_FUNC) &_EMgaussian_EMcyclecov, 3},
    {"_EMgaussian_nllcov", (DL_FUNC) &_EMgaussian_nllcov, 3},
    {"_EMgaussian_imp1matprec", (DL_FUNC) &_EMgaussian_imp1matprec, 3},
    {"_EMgaussian_imp2matprec", (DL_FUNC) &_EMgaussian_imp2matprec, 3},
    {"_EMgaussian_EMcycleprec", (DL_FUNC) &_EMgaussian_EMcycleprec, 3},
    {"_EMgaussian_nllprec", (DL_FUNC) &_EMgaussian_nllprec, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_EMgaussian(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
