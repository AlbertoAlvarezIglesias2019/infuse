#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

//' Predict Probabilities (Strictly Less Than)
//'
//' @description Map new numeric values to cumulative probabilities using
//' the step function defined by previously calculated ECDF or Kaplan-Meier results.
//' This function implements a left-continuous mapping: \eqn{F(z-) = P(X < z)}.
//'
//' @param z A numeric vector of query values.
//' @param x_sorted A numeric vector of observed values, sorted in non-decreasing order.
//' @param p_sorted A numeric vector of cumulative probabilities corresponding to \code{x_sorted}.
//'
//' @return A numeric vector representing \eqn{P(X < z)}.
//'
//' @examples
//'
//' # Basic Usage
//' x <- c(4, 3, 2, 5, 34, 52, 64, 87, 23)
//' tf <- c_ecdf_plus(x)
//' x_s <- x[tf$o]
//' p_s <- tf$p[tf$o]
//' c_ecdf_predict(3.01,x_s, p_s)
//' c_ecdf_predict(3,x_s, p_s)
//' c_ecdf_predict(2.99,x_s, p_s)
//'
//' c_ecdf_predict_less(3.01,x_s, p_s)
//' c_ecdf_predict_less(3,x_s, p_s)
//' c_ecdf_predict_less(2.99,x_s, p_s)
//'
//' @export
// [[Rcpp::export]]
 NumericVector c_ecdf_predict_less(NumericVector z, NumericVector x_sorted, NumericVector p_sorted) {
   int nz = z.size();
   int n = x_sorted.size();
   NumericVector out(nz);

   // 1. Find the boundary of non-NA values (O(n) once)
   int n_valid = 0;
   for (int j = 0; j < n; ++j) {
     if (NumericVector::is_na(x_sorted[j])) break;
     n_valid++;
   }

   // If there are no valid numbers, return NAs
   if (n_valid == 0) {
     return rep(NA_REAL, nz);
   }

   for(int i = 0; i < nz; ++i) {
     if (NumericVector::is_na(z[i])) {
       out[i] = NA_REAL;
       continue;
     }

     // lower_bound finds the first element >= z
     auto it = std::lower_bound(x_sorted.begin(), x_sorted.begin() + n_valid, z[i]);

     if (it == x_sorted.begin()) {
       // All observed x are >= z, so no values are strictly less than z
       out[i] = 0.0;
     } else {
       // Move back one step to get the element where x_j < z
       int idx = std::distance(x_sorted.begin(), it) - 1;
       out[i] = p_sorted[idx];
     }
   }

   return out;
 }

