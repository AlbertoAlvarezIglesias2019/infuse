#include <Rcpp.h>
using namespace Rcpp;

//' Find Interval Indices (C++ Implementation)
//'
//' @description For each element in x, this function finds the index of the
//' interval in breaks that contains the value. It uses a binary search
//' algorithm (\code{std::upper_bound}) for high performance.
//'
//' @param x A numeric vector of values to be localized.
//' @param breaks A numeric vector of sorted breakpoints defining the intervals.
//'
//' @return An integer vector of the same length as \code{x}. Each value
//' represents the number of breakpoints that are less than or equal to the
//' corresponding value in \code{x}.
//'
//' @details
//' The function mimics the behavior of \code{findInterval(x, breaks)} in
//' base R. Specifically:
//' \itemize{
//'   \item If \code{x[i] < breaks[1]}, it returns 0.
//'   \item If \code{breaks[j] <= x[i] < breaks[j+1]}, it returns j.
//'   \item If \code{x[i] >= breaks[last]}, it returns the length of breaks.
//' }
//'
//' Because it uses \code{std::upper_bound}, the complexity is \eqn{O(N \cdot \log M)},
//' where \eqn{N} is the length of \code{x} and \eqn{M} is the length of \code{breaks}.
//'
//' @examples
//' \dontrun{
//' x <- c(0, 5, 12, 25)
//' breaks <- c(1, 10, 20)
//' # Returns: c(0, 1, 1, 3)
//' winDwinR::c_findInterval(x, breaks)
//'
//' findInterval(x, breaks)
//' }
//'
//' # CHECK SPEED
//' x <- rnorm(1000)
//' breaks <- rnorm(100)
//' # Time the C++ Wrapper
//' start_time <- Sys.time()
//'
//'    res_cpp <- winDwinR::c_findInterval(x,breaks)
//'
//' end_time <- Sys.time()
//' t1 <- as.numeric(end_time - start_time)
//'
//' start_time <- Sys.time()
//'
//'    res_r <- findInterval(x,sort(breaks) )
//'
//' end_time <- Sys.time()
//' t2 <- as.numeric(end_time - start_time)
//' cat("\n Speed C++ ",t1,"\n Speed R ",t2,"\n Ratio",t2/t1,"\n")
//'
//' @export
//'
// [[Rcpp::export]]
IntegerVector c_findInterval(NumericVector x, NumericVector breaks) {

  // Use clone() to create a copy that belongs ONLY to this function
  NumericVector sorted_breaks = clone(breaks);

  // 2. Sort the breaks (O(M log M))
  // This ensures std::upper_bound (binary search) works correctly
  std::sort(sorted_breaks.begin(), sorted_breaks.end());


 IntegerVector out(x.size());

 NumericVector::iterator it, pos;
 IntegerVector::iterator out_it;

 for(it = x.begin(), out_it = out.begin(); it != x.end();
 ++it, ++out_it) {
   // std::upper_bound performs the binary search
   pos = std::upper_bound(sorted_breaks.begin(), sorted_breaks.end(), *it);
   // std::distance calculates the index based on the iterator position
   *out_it = std::distance(sorted_breaks.begin(), pos);
 }

 return out;
}
