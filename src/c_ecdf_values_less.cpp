#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

//' Map Values to Empirical Distribution Steps (Strictly Less Than)
//'
//' @description Evaluates a step function (like an ECDF) at specific points
//' using the logic \eqn{P(X < z)}. This results in a left-continuous
//' function, where the value at a jump point is the value of the lower step.
//'
//' @param uni A sorted numeric vector of unique time points or breakpoints.
//' @param w A numeric vector of weights or probabilities associated with
//' each interval.
//' @param z A numeric vector of points at which the function should be evaluated.
//' @param ofset A numeric value subtracted from each element of \code{z}.
//' Defaults to 0.
//'
//' @return A numeric vector of the same length as \code{z} containing the
//' mapped values from \code{w}.
//'
//' @details
//' This function uses \code{std::lower_bound} to find the first index where
//' \code{uni >= target}. Consequently, it maps to the weight associated with
//' values strictly less than the target.
//'
//' For a given target \eqn{t = z_i - \text{ofset}}, the result is:
//' \itemize{
//'   \item 0 if \eqn{t \le \text{uni}[0]}
//'   \item \eqn{w[j-1]} where \eqn{j} is the index such that \eqn{\text{uni}[j-1] < t \le \text{uni}[j]}
//' }
//'
//' @examples
//' \dontrun{
//' # Standard ECDF (X <= z) vs Strict ECDF (X < z)
//' times <- c(1, 3, 5)
//' weights <- c(0.2, 0.5, 0.9)
//'
//' # At the jump point (z = 3):
//' # c_ecdf_values would return 0.5
//' # c_ecdf_values_less returns 0.2
//' c_ecdf_values_less(times, weights, 3, ofset = 0)
//' }
//'
//' @export
// [[Rcpp::export]]
NumericVector c_ecdf_values_less(NumericVector uni,NumericVector w, NumericVector z,double ofset=0){

  int nz = z.size();
  int n_uni = uni.size();
  NumericVector out(nz);

  // We loop once through z
  for (int i = 0; i < nz; ++i) {

    // 1. Check for NA/NaN values first
    if (NumericVector::is_na(z[i])) {
      out[i] = NA_REAL;
      continue;
    }

    // 1. Apply offset 'on the fly' to avoid creating the zz vector
    double target = z[i] - ofset;

    // 2. Perform binary search directly here
    // This is what findInterval does internally
    auto pos = std::lower_bound(uni.begin(), uni.end(), target);
    int idx = std::distance(uni.begin(), pos);

    // 3. Map to w (equivalent to your if/else logic)
    if (idx == 0) {
      out[i] = 0;
    } else if (idx == n_uni) {
      // If target is >= max(uni), return 1
      out[i] = 1.0;
    } else {
      out[i] = w[idx - 1];
    }
  }

  return out;
}
