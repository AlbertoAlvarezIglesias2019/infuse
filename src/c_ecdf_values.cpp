#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

//' Map Values to Empirical Distribution Steps
//'
//' @description This function evaluates a step function (like an ECDF or Kaplan-Meier
//' estimate) at specific points. For each value in \code{z}, it finds the
//' corresponding value in \code{w} based on the intervals defined by \code{uni}.
//'
//' @param uni A sorted numeric vector of unique time points or breakpoints.
//' @param w A numeric vector of weights or probabilities associated with each
//' interval in \code{uni}. Must have the same length as \code{uni}.
//' @param z A numeric vector of points at which the function should be evaluated.
//' @param ofset A numeric value subtracted from each element of \code{z}
//' before evaluation. Defaults to 0.
//'
//' @return A numeric vector of the same length as \code{z} containing the
//' mapped values from \code{w}.
//'
//' @details
//' The function performs an "on-the-fly" lookup using a binary search
//' (\code{std::upper_bound}). For a given target \eqn{t = z_i - \text{ofset}},
//' the result is:
//' \itemize{
//'   \item 0 if \eqn{t < \text{uni}[0]}
//'   \item \eqn{w[j-1]} where \eqn{j} is the index such that \eqn{\text{uni}[j-1] \le t < \text{uni}[j]}
//' }
//' This is optimized to avoid temporary vector allocations and is much faster
//' than manual looping in R.
//'
//' @examples
//' \dontrun{
//' times <- c(1, 3, 5)
//' probs <- c(0.2, 0.5, 0.9)
//' points <- c(0, 2, 4, 6)
//' # Map points to probabilities with an offset
//' c_ecdf_values(times, probs, points, ofset = 0)
//' # Returns: c(0, 0.2, 0.5, 0.9)
//' }
//'
//' @export
// [[Rcpp::export]]
NumericVector c_ecdf_values(NumericVector uni,NumericVector w, NumericVector z,double ofset=0){

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
    auto pos = std::upper_bound(uni.begin(), uni.end(), target);
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
