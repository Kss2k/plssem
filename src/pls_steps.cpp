// PLS-SEM outer-weight iteration (steps 0-5) in C++.
//
// Steps 0-5 depend only on the indicator correlation matrix `S` and the model
// graph -- never on the raw data -- so their cost is independent of the sample
// size. In R they are thousands of tiny matrix operations per fit, which is
// interpreter overhead rather than arithmetic; MC-PLS runs this loop
// `mc.reps` times per Robbins-Monro iteration, so it dominates. This is a
// direct transcription of `estimatePLS_Step0()` .. `estimatePLS_Step5()`.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// C = L'SL ; SC = [[S, SL], [ (SL)', C ]]     (steps 0 and 4)
static inline void recomputeCSC(const arma::mat& S, const arma::mat& lambda,
                                arma::mat& C, arma::mat& SC) {
  const arma::uword p = S.n_rows, k = lambda.n_cols;
  const arma::mat SL = S * lambda;
  C = lambda.t() * SL;
  SC.submat(0, 0, p - 1, p - 1)         = S;
  SC.submat(0, p, p - 1, p + k - 1)     = SL;
  SC.submat(p, 0, p + k - 1, p - 1)     = SL.t();
  SC.submat(p, p, p + k - 1, p + k - 1) = C;
}

// `getNonZeroElems()`: column-major order, dropping NA and exact zeros.
static inline arma::vec nonZeroElems(const arma::mat& x) {
  std::vector<double> out;
  out.reserve(x.n_elem);
  for (arma::uword j = 0; j < x.n_cols; ++j)
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      const double v = x(i, j);
      if (!std::isnan(v) && v != 0.0) out.push_back(v);
    }
  return arma::vec(out.data(), out.size());
}

// [[Rcpp::export]]
List plsStep05Cpp(const arma::mat& S,
                  arma::mat lambda,
                  arma::mat gamma,
                  arma::mat C,
                  arma::mat SC,
                  const List& indsLvs,       // 0-based row indices into S
                  const arma::uvec& lvCol,   // 0-based column of each linear LV
                  const arma::uvec& modeB,   // 1 = mode B, 0 = mode A
                  const arma::umat& preds,
                  const arma::umat& succs,
                  arma::vec outerWeights,
                  const double tol,
                  const int maxIter) {

  const arma::uword p = S.n_rows, k = lambda.n_cols, nlin = lvCol.n_elem;

  std::vector<arma::uvec> inds(nlin);
  for (arma::uword a = 0; a < nlin; ++a)
    inds[a] = as<arma::uvec>(indsLvs[a]);

  // ---- Step 0: unit weights, normalised with the INCOMING (stale) SC ----
  for (arma::uword a = 0; a < nlin; ++a) {
    if (inds[a].n_elem == 0) continue;
    const double s = arma::accu(SC.submat(inds[a], inds[a]));
    lambda.submat(inds[a], arma::uvec{lvCol[a]}).fill(1.0 / std::sqrt(s));
  }
  recomputeCSC(S, lambda, C, SC);

  int iter = 0;
  bool converged = false;

  for (int it = 1; it <= maxIter; ++it) {
    iter = it;

    // ---- Step 1: inner (path) weights ----
    for (arma::uword a = 0; a < nlin; ++a) {
      const arma::uword c = lvCol[a];

      for (arma::uword b = 0; b < nlin; ++b)
        if (succs(b, c)) gamma(lvCol[b], c) = C(c, lvCol[b]);

      std::vector<arma::uword> pv;
      for (arma::uword b = 0; b < nlin; ++b)
        if (preds(b, c)) pv.push_back(lvCol[b]);

      if (!pv.empty()) {
        const arma::uvec pr(pv.data(), pv.size());
        const arma::uvec scpr = pr + p;
        gamma.submat(pr, arma::uvec{c}) =
          arma::solve(SC.submat(scpr, scpr), SC.submat(scpr, arma::uvec{p + c}));
      }

      const double sf =
        std::sqrt(arma::as_scalar(gamma.col(c).t() * C * gamma.col(c)));
      if (sf != 0.0) gamma.col(c) /= sf;
    }

    // ---- Step 2: rotate C/SC by the inner weights ----
    if (gamma.n_rows > 1) {
      const arma::mat SLg = (S * lambda) * gamma;
      C = gamma.t() * C * gamma;
      SC.submat(0, 0, p - 1, p - 1)         = S;
      SC.submat(0, p, p - 1, p + k - 1)     = SLg;
      SC.submat(p, 0, p + k - 1, p - 1)     = SLg.t();
      SC.submat(p, p, p + k - 1, p + k - 1) = C;
    }

    // ---- Step 3: outer weights ----
    for (arma::uword a = 0; a < nlin; ++a) {
      if (inds[a].n_elem == 0) continue;
      const arma::uword c = lvCol[a];
      const arma::mat Sjj = SC.submat(inds[a], inds[a]);
      arma::vec wj = SC.submat(inds[a], arma::uvec{p + c});
      if (modeB[a]) wj = arma::solve(Sjj, wj);
      wj /= std::sqrt(arma::as_scalar(wj.t() * Sjj * wj));
      lambda.submat(inds[a], arma::uvec{c}) = wj;
    }

    // ---- Step 4 ----
    recomputeCSC(S, lambda, C, SC);

    // ---- Step 5: convergence on the outer weights ----
    const arma::vec nw = nonZeroElems(lambda);
    bool cv = (outerWeights.n_elem == nw.n_elem);
    if (cv)
      for (arma::uword i = 0; i < nw.n_elem; ++i)
        if (!(std::fabs(outerWeights[i] - nw[i]) < tol)) { cv = false; break; }
    outerWeights = nw;

    if (cv) { converged = true; break; }
  }

  return List::create(_["lambda"] = lambda, _["gamma"] = gamma, _["C"] = C,
                      _["SC"] = SC, _["outerWeights"] = outerWeights,
                      _["convergence"] = converged, _["iterations"] = iter);
}

// ---------------------------------------------------------------------------
// Step 6: factor scores, interaction/quadratic terms, and their covariance.
//
// F = X W, then each product term is the elementwise product of two score
// columns, then C = cov(F) (denominator n - 1, matching Rfast::cova).
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List plsStep6Cpp(const arma::mat& X,          // n x p indicator data
                       const arma::mat& W,          // p x k outer weights
                       const arma::umat& prodElems, // 2 x m, 0-based score cols
                       const arma::uvec& prodCol,   // 0-based target col
                       const bool doStandardise) {

  arma::mat F = X * W;

  if (doStandardise) {
    for (arma::uword j = 0; j < F.n_cols; ++j) {
      const double mu = arma::mean(F.col(j));
      const double sd = arma::stddev(F.col(j));      // n - 1 denominator
      F.col(j) -= mu;
      if (sd > 0.0) F.col(j) /= sd;
    }
  }

  for (arma::uword t = 0; t < prodCol.n_elem; ++t)
    F.col(prodCol[t]) = F.col(prodElems(0, t)) % F.col(prodElems(1, t));

  return Rcpp::List::create(Rcpp::_["F"] = F, Rcpp::_["C"] = arma::cov(F));
}

// ---------------------------------------------------------------------------
// Ordinalisation: bucket a numeric vector by its own empirical quantiles.
//
// Equivalent to `findInterval(x, collapse::fquantile(x, probs))` with type-7
// quantiles. MC-PLS runs this over every simulated column on every
// Robbins-Monro iteration, and the R version pays for a full quantile
// computation plus a separate vectorised binary search.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::IntegerVector ordinalizeCpp(const Rcpp::NumericVector& x,
                                  Rcpp::NumericVector probs) {

  const R_xlen_t n = x.size();
  Rcpp::IntegerVector out(n);
  const double* px = REAL(x);

  std::vector<double> p;
  p.reserve(probs.size());
  for (R_xlen_t i = 0; i < probs.size(); ++i)
    if (probs[i] < 1.0) p.push_back(probs[i]);
  std::sort(p.begin(), p.end());

  const std::size_t k = p.size();
  if (k == 0 || n == 0) return out;

  std::vector<double> buf;
  buf.reserve(n);
  bool anyNA = false;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (ISNAN(px[i])) { anyNA = true; continue; }
    buf.push_back(px[i]);
  }
  if (buf.empty()) {
    for (R_xlen_t i = 0; i < n; ++i) out[i] = NA_INTEGER;
    return out;
  }

  // Only the order statistics bracketing each type-7 quantile are needed, so
  // select them with nth_element (linear) rather than sorting the whole vector.
  const double m = static_cast<double>(buf.size());
  const std::size_t last = buf.size() - 1;
  std::vector<double> h(k);
  std::vector<std::size_t> need;
  need.reserve(2 * k);
  for (std::size_t j = 0; j < k; ++j) {
    h[j] = (m - 1.0) * p[j];
    std::size_t li = static_cast<std::size_t>(std::floor(h[j]));
    if (li > last) li = last;
    need.push_back(li);
    if (li < last) need.push_back(li + 1);
  }
  std::sort(need.begin(), need.end());
  need.erase(std::unique(need.begin(), need.end()), need.end());

  std::size_t lo = 0;
  for (std::size_t t = 0; t < need.size(); ++t) {
    std::nth_element(buf.begin() + lo, buf.begin() + need[t], buf.end());
    lo = need[t] + 1;
  }

  std::vector<double> brk(k);
  for (std::size_t j = 0; j < k; ++j) {
    const double fl = std::floor(h[j]);
    std::size_t li = static_cast<std::size_t>(fl);
    if (li >= last) { brk[j] = buf[last]; continue; }
    brk[j] = buf[li] + (h[j] - fl) * (buf[li + 1] - buf[li]);
  }

  // findInterval(): count of breaks <= x. `k` is small and the breaks are
  // sorted, so sum the comparisons branchlessly instead of looping/searching.
  int* po = INTEGER(out);
  const double* pb = brk.data();
  if (anyNA) {
    for (R_xlen_t i = 0; i < n; ++i) {
      const double v = px[i];
      if (ISNAN(v)) { po[i] = NA_INTEGER; continue; }
      int c = 0;
      for (std::size_t j = 0; j < k; ++j) c += (pb[j] <= v);
      po[i] = c;
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      const double v = px[i];
      int c = 0;
      for (std::size_t j = 0; j < k; ++j) c += (pb[j] <= v);
      po[i] = c;
    }
  }
  return out;
}
