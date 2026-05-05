#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

//' Predict Probabilities from ECDF or Kaplan-Meier Results
//'
//' @description Map new numeric values to cumulative probabilities using
//' the step function defined by previously calculated ECDF or Kaplan-Meier results.
//' This function implements a right-continuous mapping: \eqn{F(z) = P(X \le z)}.
//'
//' @param z A numeric vector of query values (points at which to evaluate the function).
//' @param x_sorted A numeric vector of observed values, sorted in non-decreasing order.
//' Usually obtained via \code{x[eee$o]} from \code{c_ecdf_plus}.
//' @param p_sorted A numeric vector of cumulative probabilities corresponding to \code{x_sorted}.
//' Usually obtained via \code{p[eee$o]} from \code{c_ecdf_plus}.
//'
//' @return A numeric vector of the same length as \code{z} containing the
//' mapped cumulative probabilities. Returns \code{0.0} for values smaller than
//' the minimum of \code{x_sorted} and the last calculated probability for
//' values larger than the maximum of \code{x_sorted}.
//'
//' @details
//' The function uses \code{std::upper_bound} (binary search) to achieve
//' \eqn{O(\log n)} lookup time per query point, making it highly efficient for
//' large datasets. It correctly handles \code{NA} values in the query vector \code{z}.
//'
//' @examples
//' # 1. Fit the model
//' x <- c(1, 2, 3, 4, 5)
//' eee <- c_ecdf_plus(x)
//'
//' # 2. Prepare sorted inputs for prediction
//' x_s <- x[eee$o]
//' p_s <- eee$p[eee$o]
//'
//' # 3. Predict for new values
//' z_new <- c(0.5, 2, 2.5, 10)
//' c_ecdf_predict(z_new, x_s, p_s)
//' # Expected: 0.0, 0.4, 0.4, 1.0
//'
//' times <- c(1, 3, 5)
//' probs <- c(0.2, 0.5, 0.9)
//' points <- c(0, 2, 4, 6)
//' # Map points to probabilities with an offset
//' c_ecdf_predict(points,times, probs)
//' # Returns: c(0, 0.2, 0.5, 0.9)
//'
//'
//' # Basic Usage
//' x <- c(4, 3, 2, 5, 34, 52, 64, 87, 23)
//' tf <- c_ecdf_plus(x)
//' x_s <- x[tf$o]
//' p_s <- tf$p[tf$o]
//' c_ecdf_predict(3,x_s, p_s)
//' c_ecdf_predict(2.99,x_s, p_s)
//'
//' # Comparison with base R ecdf
//' z <- c(2, 4, 23, 100)
//' base_r <- ecdf(x)(z)
//' custom_cpp <- c_ecdf_predict(z,x_s, p_s)
//' all.equal(base_r, custom_cpp)
//'
//'
//' @export
// [[Rcpp::export]]
NumericVector c_ecdf_predict(NumericVector z, NumericVector x_sorted, NumericVector p_sorted) {
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


  // Ensure we handle the step function logic: P(X <= z)
  for(int i = 0; i < nz; ++i) {
    if (NumericVector::is_na(z[i])) {
      out[i] = NA_REAL;
      continue;
    }

    // upper_bound finds the first element > z
    auto it = std::upper_bound(x_sorted.begin(), x_sorted.begin() + n_valid, z[i]);

    if (it == x_sorted.begin()) {
      out[i] = 0.0; // z is smaller than any observed x
    } else {
      // Move back one step to get the element where x_j <= z
      int idx = std::distance(x_sorted.begin(), it) - 1;
      out[i] = p_sorted[idx];
    }
  }

  return out;
}
