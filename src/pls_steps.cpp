// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


static inline void outerUpdateCovarianceMatrices(
  const arma::mat& S,
  const arma::mat& lambda,
  arma::mat& C,
  arma::mat& SC
) {
  const arma::uword p = S.n_rows, k = lambda.n_cols;
  const arma::mat S_lower = S * lambda;

  C = lambda.t() * S_lower;
  SC.submat(0, 0, p - 1, p - 1)         = S;
  SC.submat(0, p, p - 1, p + k - 1)     = S_lower;
  SC.submat(p, 0, p + k - 1, p - 1)     = S_lower.t();
  SC.submat(p, p, p + k - 1, p + k - 1) = C;
}


// [[Rcpp::export]]
Rcpp::List estimatePLS_Step0_5_Cpp(
  arma::mat lambda,
  arma::mat gamma,
  arma::mat S,
  arma::mat C,
  arma::mat SC,
  const Rcpp::List& R_IndsIdxLVs,
  const arma::uvec& lvColIdx,
  const arma::uvec& modeB,
  const arma::umat& preds, // model@matrices$preds.cfa for CFA models
  const arma::umat& succs, // model@matrices$succs.cfa for CFA models
  const double tolerance,
  const int maxiter
) {

  const arma::uword p = S.n_rows, k = lambda.n_cols, nlin = lvColIdx.n_elem;

  std::vector<arma::uvec> indsIdxLVs(nlin);
  for (arma::uword j = 0; j < nlin; j++) {
    indsIdxLVs[j] = Rcpp::as<arma::uvec>(R_IndsIdxLVs[j]);
  }
 
  // Step 0
  for (arma::uword j = 0; j < nlin; j++) {
    const arma::uvec J = indsIdxLVs[j];
    if (J.is_empty()) continue;

    const arma::uword l = lvColIdx[j];

    const arma::vec wj(J.n_elem, arma::fill::ones);
    const arma::mat Sjj = SC.submat(J, J);
   
    const double std = std::sqrt(arma::as_scalar(wj.t() * Sjj * wj));
    if (std > 0.0) lambda.submat(J, arma::uvec{l}) = wj / std;
  }

  // Collected outer weights
  arma::vec w0 = vectorise(lambda);
  arma::vec w1 = w0;

  outerUpdateCovarianceMatrices(S, lambda, C, SC);

  bool converged = false;
  int iter;
  for (iter = 0; iter < maxiter; iter++) {
    
    // Step 1
    for (arma::uword j = 0; j < nlin; j++) {
      const arma::uword l = lvColIdx[j];
      const arma::uvec L{l};

      const arma::uvec pj = arma::find(preds.col(l) > 0);
      const arma::uvec sj = arma::find(succs.col(l) > 0);

      for (arma::uword i = 0; i < sj.n_elem; i++) {
        gamma(sj[i], l) = C(l, sj[i]);
      }

      if (pj.n_elem) {
        gamma.submat(pj, L) = arma::solve(
          C.submat(pj, pj), C.submat(pj, L)
        );
      }

      const double normalize = std::sqrt(arma::as_scalar(
        gamma.col(l).t() * C * gamma.col(l) 
      ));

      if (normalize > 0)
        gamma.col(l) = gamma.col(l) / normalize;
    }

    // Step 2
    if (gamma.n_rows > 1) {
      const arma::mat S_LambdaGamma = (S * lambda) * gamma;
      C = gamma.t() * C * gamma;
      SC.submat(0, 0, p - 1, p - 1)         = S;
      SC.submat(0, p, p - 1, p + k - 1)     = S_LambdaGamma;
      SC.submat(p, 0, p + k - 1, p - 1)     = S_LambdaGamma.t();
      SC.submat(p, p, p + k - 1, p + k - 1) = C;
    }

    // Step 3
    for (arma::uword j = 0; j < nlin; j++) {
      if (indsIdxLVs[j].n_elem <= 0) continue;

      const arma::uvec J = indsIdxLVs[j];
      const arma::uword l = lvColIdx[j];
      
      const arma::mat Sjj = S.submat(J, J);


      // Weights mode A
      arma::vec wj = SC.submat(J, arma::uvec{p + l});

      if (modeB[j]) {
        // Mode B uses weights from mode A
        wj = arma::solve(Sjj, wj);
      }

      wj /= std::sqrt(arma::as_scalar(wj.t() * Sjj * wj));
      lambda.submat(J, arma::uvec{l}) = wj; 
    }

    // Step 4
    outerUpdateCovarianceMatrices(S, lambda, C, SC);

    // Step 5
    w1 = vectorise(lambda);
    const arma::vec diff = arma::abs(w1 - w0);

    if (diff.is_empty() || !diff.is_finite()) {
      converged = false;
      break;
    }

    if (diff.max() <= tolerance) {
      converged = true;
      break;
    }

    w0 = w1;
  }

  return Rcpp::List::create(
    Rcpp::_["lambda"]      = lambda,
    Rcpp::_["gamma"]       = gamma,
    Rcpp::_["C"]           = C,
    Rcpp::_["SC"]          = SC,
    Rcpp::_["convergence"] = converged,
    Rcpp::_["iterations"]  = iter < maxiter ? iter + 1L : maxiter
  );
}


// [[Rcpp::export]]
Rcpp::List estimatePLS_Step6_Cpp(
  const arma::mat& X,
  const arma::mat& W,
  const Rcpp::List& prodElemsIdx, // g x m, 0-indexed score cols
  const arma::uvec& prodColIdx,   // 0-indexed target col
  const bool standardize
) {
  arma::mat F = X * W;
        
  if (standardize) for (arma::uword j = 0; j < F.n_cols; ++j) {
    const double mu = arma::mean(F.col(j));
    const double sigma = arma::stddev(F.col(j));

    F.col(j) -= mu;
    if (sigma > 0.0) F.col(j) /= sigma;
  }

  for (arma::uword t = 0; t < prodColIdx.n_elem; ++t) {
    arma::uvec elems = Rcpp::as<arma::uvec>(prodElemsIdx[(int)t]);
    F.col(prodColIdx[t]) = F.col(elems[0]) % F.col(elems[1]);

    for (arma::uword c = 2; c < elems.n_elem; c++) // n-way interaction terms
      F.col(prodColIdx[t]) = F.col(prodColIdx[t]) % F.col(elems[c]);
  }

  return Rcpp::List::create(Rcpp::_["F"] = F, Rcpp::_["C"] = arma::cov(F));
}
