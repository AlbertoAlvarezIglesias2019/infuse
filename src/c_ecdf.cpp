#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

//' Empirical Cumulative Distribution Function (ECDF)
//'
//' @description Calculates the jump points of the empirical cumulative 
//' distribution function. It identifies unique sorted values and their 
//' corresponding cumulative probabilities.
//'
//' @param x A numeric vector of observations.
//'
//' @return A named list containing two components:
//' \itemize{
//'   \item \code{uni}: A numeric vector of the unique sorted values of \code{x}.
//'   \item \code{w}: A numeric vector of the cumulative probabilities \eqn{P(X \le x)}.
//' }
//'
//' @details
//' The function sorts the input vector and identifies unique values. The 
//' probability for each unique value \eqn{x_j} is calculated as \eqn{F(x_j) = i/n}, 
//' where \eqn{i} is the count of all elements \eqn{\le x_j}. This is equivalent 
//' to the "knots" and "weights" of a standard ECDF.
//'
//' @examples
//' \dontrun{
//' x <- c(1, 2, 2, 3, 5)
//' res <- c_ecdf(x)
//' 
//' # Compare with base R ecdf()
//' fit_base <- ecdf(x)
//' 
//' # Check unique values (knots)
//' all.equal(res$uni, knots(fit_base))
//' 
//' # Check probabilities at those knots
//' all.equal(res$w, fit_base(knots(fit_base)))
//' 
//' # Example with more data
//' set.seed(123)
//' y <- rnorm(100)
//' res_y <- c_ecdf(y)
//' plot(res_y$uni, res_y$w, type = "s", main = "C++ ECDF")
//' }
//' 
//' data_points <- c(5, 1, 3, 1, 5, 5)
//' result <- c_ecdf(data_points)
//'
//' # Access unique values
//' print(result$uni)
//'
//' # Access cumulative weights
//' print(result$w)
//'
//' # compare to ecdf
//' x <- c(1,2,4,2,2,35,6,3,2,3,1,3)
//' c_ecdf(x)
//' unique(x)
//' ecdf(x)(unique(x))
//'
//' # CHECK SPEED
//' large_data <- rnorm(100000)
//' # Time the C++ Wrapper
//' start_time <- Sys.time()
//'
//'    res_cpp <- c_ecdf(large_data)
//'
//' end_time <- Sys.time()
//' t1 <- as.numeric(end_time - start_time)
//'
//' start_time <- Sys.time()
//'
//'    res_r <- ecdf(large_data)(unique(large_data))
//'
//' end_time <- Sys.time()
//' t2 <- as.numeric(end_time - start_time)
//' cat("\n Speed C++ ",t1,"\n Speed R ",t2,"\n Ratio",t2/t1,"\n")
//' @export
// [[Rcpp::export]]
List c_ecdf(NumericVector x){

  int n_orig = x.size();
  
  // 1. Filter out NAs first
  // LogicalVector is_na = is_na(x); // Rcpp way
  // However, for speed, we can do a single pass:
  std::vector<double> clean_x;
  clean_x.reserve(n_orig);

  for(int i = 0; i < x.size(); ++i) {
    if (!NumericVector::is_na(x[i])) {
      clean_x.push_back(x[i]);
    }
  }

  int n = clean_x.size();
  //if (n == 0) return NumericVector(0);
  if (n == 0) {
    return List::create(Named("uni") = NumericVector(0),
                        Named("w")   = NumericVector(0));
  }
  
  // 1. Sort the input
  std::sort(clean_x.begin(),clean_x.end());

  // 2. Identify unique values and counts
  // We'll store unique values and the cumulative index
  //NumericVector xu(n);
  //NumericVector wu(n);
  // 3. Identify unique values and counts
  // Using std::vector here is safer and cleaner when using push_back
  std::vector<double> xu;
  std::vector<double> wu;
  xu.reserve(n);
  wu.reserve(n);

  //int j = 0;
  for (int i = 0; i < n; ++i) {
    // If it's the last occurrence of a specific value
    if (i == n - 1 || clean_x[i] != clean_x[i+1]) {
      //xu[j] = clean_x[i];
      //// Cumulative probability is (index + 1) / total N
      //wu[j] = static_cast<double>(i + 1) / n;
      //j++;
      xu.push_back(clean_x[i]);
      // Cumulative probability is (current index + 1) / total N
      wu.push_back(static_cast<double>(i + 1) / n);
    }
  }

  // 3. Generate compact output
  //NumericVector out(2 * j);
  //for (int i = 0; i < j; ++i) {
  //  out[i] = xu[i];         // Values
  //  out[i + j] = wu[i];     // Probabilities
  //}

  //std::cout << j << std::endl;

  //return out;

  // 4. Return as a named R list
  return List::create(
    Named("uni") = wrap(xu),
    Named("w")   = wrap(wu)
  );
  
}
