#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// Forward Declarations
List c_ecdf_surv(NumericVector ttt, IntegerVector sss);
NumericVector c_ecdf_values(NumericVector uni, NumericVector w, NumericVector z, double ofset = 0);

//' Calculate the Weighted Mean Survival Time
//' 
//' @description A C++ implementation that calculates the mean by mapping 
//' sorted observations to their Kaplan-Meier cumulative incidence weights.
//' 
//' @param x A numeric vector of observed times.
//' @param xs An integer vector of status indicators (1 = event, 0 = censored).
//' 
//' @return A double representing the estimated mean.
//' 
//' @examples
//' \dontrun{
//' x <- c(0.2, 1, 2, 4, 2, 2, 35, 6, 3, 2, 3, 1, 3, 1.5, 45, 23, 12)
//' xs <- c(1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0)
//' c_mean(x, xs)
//' 
//' x <- rnorm(1000,40,2)
//' xs <- rbinom(1000,size=1,prob=0.6)
//' c_mean(x, xs)
//' 
//' x <- rnorm(1000,40,2)
//' xs <- rep(1,1000)
//' c_mean(x, xs)
//' mean(x)
//' }
//' 
//' @export
// [[Rcpp::export]]
 double c_mean(NumericVector x, IntegerVector xs) {
   // 1. Get the KM survival / Cumulative incidence list
   List km = c_ecdf_surv(x, xs);
   NumericVector km_uni = km["uni"];
   NumericVector km_w = km["w"];
   
   if (km_uni.size() == 0) return NA_REAL;
   
   // 2. Sort x (Equivalent to xsort <- sort(x))
   NumericVector xsort = clone(x);
   std::sort(xsort.begin(), xsort.end());
   
   // 3. Call your existing c_ecdf_values (returns NumericVector)
   // We pass ofset = 0 as per your R code logic
   NumericVector www = c_ecdf_values(km_uni, km_w, xsort, 0.0);
   
   // 4. Calculate weighted sum (Equivalent to lag logic: www - lag(www))
   double total_mean = 0.0;
   double prev_w = 0.0;
   
   for (int i = 0; i < xsort.size(); ++i) {
     // Handle NAs (na.rm = TRUE)
     if (NumericVector::is_na(xsort[i]) || NumericVector::is_na(www[i])) {
       continue;
     }
     
     double current_w = www[i];
     double weight = current_w - prev_w;
     
     total_mean += xsort[i] * weight;
     
     // Update prev_w for the next iteration (the "lag")
     prev_w = current_w;
   }
   
   return total_mean;
 }