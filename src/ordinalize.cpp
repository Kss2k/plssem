#include "Rcpp.h"

// [[Rcpp::export]]
Rcpp::IntegerVector ordinalizeVectorCpp(
    const Rcpp::NumericVector& x,
    Rcpp::NumericVector probs,
    double ztol = 0.001
) {
  if (!R_finite(ztol) || ztol < 0 || ztol >= 0.5)
    Rcpp::stop("ztol must be finite and in [0, 0.5)");

  const R_xlen_t n = x.size();
  Rcpp::IntegerVector out(n);
  const double* px = REAL(x);

  std::vector<double> p;
  p.reserve(probs.size());
  for (R_xlen_t i = 0; i < probs.size(); ++i) {
    if (!R_finite(probs[i])) Rcpp::stop("probs must be finite");
    p.push_back(std::min(std::max(probs[i], ztol), 1.0 - ztol));
  }

  std::sort(p.begin(), p.end());

  const std::size_t k = p.size();
  Rcpp::NumericVector thresholds(k, NA_REAL);
  out.attr("tau") = thresholds;
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

  double* brk = REAL(thresholds);
  for (std::size_t j = 0; j < k; ++j) {
    const double fl = std::floor(h[j]);
    std::size_t li = static_cast<std::size_t>(fl);
    if (li >= last) { brk[j] = buf[last]; continue; }
    brk[j] = buf[li] + (h[j] - fl) * (buf[li + 1] - buf[li]);
  }

  // findInterval(): count of breaks <= x. k is small and the breaks are
  // sorted, so sum the comparisons branchlessly instead of looping/searching.
  int* po = INTEGER(out);
  const double* pb = brk;

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
