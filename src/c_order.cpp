#include <Rcpp.h>
#include <numeric>   // for std::iota
#include <algorithm> // for std::stable_sort

using namespace Rcpp;

struct IndexedValue {
  double value;
  int index;
};


//' Optimized Order Function
//'
//' @description This function mimics the behavior of R's \code{order()} function
//' using C++'s \code{std::stable_sort}. It returns the 1-based indices that
//' would result in a sorted version of the input vector.
//'
//' @param x A numeric vector to be ordered.
//'
//' @return An integer vector of 1-based indices (matching R's indexing convention).
//'
//' @details
//' The function uses \code{std::stable_sort} to ensure that the relative
//' order of elements with equal values is preserved, making it perfectly
//' compatible with the default behavior of R's \code{order(..., method = "shell")}.
//'
//' @examples
//' # Create a vector with ties
//' x <- c(3, 2, 4, 3, 2, 1, 4)
//'
//' # Get indices from C++
//' cpp_idx <- winDwinR::c_order(x)
//'
//' # Get indices from base R
//' r_idx <- order(x)
//'
//' # Compare results
//' identical(cpp_idx, r_idx)
//'
//' # Use indices to sort the original vector
//' x[cpp_idx]
//'
//'
//' # CHECK SPEED
//' large_data <- rnorm(1000)
//' # Time the C++ Wrapper
//' start_time <- Sys.time()
//'
//'    res_cpp <- winDwinR::c_order(large_data)
//'
//' end_time <- Sys.time()
//' t1 <- as.numeric(end_time - start_time)
//'
//' start_time <- Sys.time()
//'
//'    res_r <- order(large_data)
//'
//' end_time <- Sys.time()
//' t2 <- as.numeric(end_time - start_time)
//' cat("\n Speed C++ ",t1,"\n Speed R ",t2,"\n Ratio",t2/t1,"\n")
//'
//' @export
//'
// [[Rcpp::export]]
IntegerVector c_order(NumericVector x) {
  int n = x.size();
  std::vector<IndexedValue> iv(n);

  // 1. Pack values and indices together
  for (int i = 0; i < n; ++i) {
    iv[i].value = x[i];
    iv[i].index = i + 1; // Pre-convert to R's 1-based indexing
  }

  // 2. Sort the structs
  // This is faster because value and index are next to each other in memory
  std::stable_sort(iv.begin(), iv.end(), [](const IndexedValue& a, const IndexedValue& b) {
    return a.value < b.value;
  });

  // 3. Unpack indices
  IntegerVector out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = iv[i].index;
  }

  return out;
}
