#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::IntegerVector match_approx_cpp(Rcpp::NumericVector x,
				     Rcpp::NumericVector table,
				     int nomatch,
				     Rcpp::NumericVector tolerance) {
  int n = x.size();
  int m = table.size();
  Rcpp::IntegerVector out(n, nomatch);
  for (int i = 0; i < n; ++i) {
    int j = 0;
    while (j < m) {
      if (std::abs(x[i] - table[j]) <= tolerance[i]) {
        out[i] = j + 1;
        break;
      }
      ++j;
    }
  }
  return out;
}
