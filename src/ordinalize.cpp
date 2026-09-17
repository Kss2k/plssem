#include "Rcpp.h"

// [[Rcpp::export]]
Rcpp::IntegerVector ordinalizeVectorCpp(const Rcpp::NumericVector& x, Rcpp::NumericVector probs) {

  const R_xlen_t n = x.size();
  Rcpp::IntegerVector out(n);
  const double* px = REAL(x);

  std::vector<double> p;
  p.reserve(probs.size());
  for (R_xlen_t i = 0; i < probs.size(); ++i) {
    if (probs[i] < 1.0) p.push_back(probs[i]);
  }

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

  // findInterval(): count of breaks <= x. k is small and the breaks are
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
