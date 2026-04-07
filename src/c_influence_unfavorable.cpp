#include <Rcpp.h>
using namespace Rcpp;

// --- Forward Declarations ---
NumericVector c_ecdf_values(NumericVector uni, NumericVector w, NumericVector z, double ofset=0);
NumericVector c_ecdf_values_less(NumericVector uni, NumericVector w, NumericVector z, double ofset=0);

//' Calculate Loss Probability Point Estimate and Influence Functions
//'
//' This function computes the point estimate for the probability of an unfavorable
//' pair (Loss Probability, P(Y >= X + lambda)) and returns the influence functions
//' for both the treatment (X) and control (Y) groups.
//'
//' @param x NumericVector: Raw values of the treatment group.
//' @param y NumericVector: Raw values of the control group.
//' @param x_uni NumericVector: Sorted unique values of the treatment group (X).
//' @param y_uni NumericVector: Sorted unique values of the control group (Y).
//' @param x_w NumericVector: Cumulative ECDF weights for x_uni.
//' @param y_w NumericVector: Cumulative ECDF weights for y_uni.
//' @param lambda double: The clinically meaningful threshold for defining a loss.
//'
//' @return List: A list containing:
//'    \item{theta}{Point Estimate: P(Y >= X + lambda)}
//'    \item{IF_x}{Influence Function values for the treatment group}
//'    \item{IF_y}{Influence Function values for the control group}
//'
//' @details
//' The loss probability theta is calculated by integrating the ECDF of X evaluated
//' at (Y - lambda). The influence functions represent the centered contribution
//' of each observation to this probability.
//'
//' @export
// [[Rcpp::export]]
List c_influence_unfavorable(NumericVector x, NumericVector y,
                       NumericVector x_uni, NumericVector y_uni,
                       NumericVector x_w, NumericVector y_w,
                       double lambda) {

   int nx = x.size();
   int ny = y.size();

   // --- PART 1: Point Estimate Calculation ---
   // theta = P(Y >= X + lambda) = E_y[ F_x(y_j - lambda) ]
   NumericVector fx_yj_minus = c_ecdf_values(x_uni, x_w, y, lambda);
   NumericVector fy_xi_plus_less = c_ecdf_values_less(y_uni, y_w, x, -lambda);
   
   // Weights for y to integrate
   NumericVector weights_x = c_ecdf_values(x_uni, x_w, x);
   double theta = 0.0;

   for(int i = 0; i < nx; ++i) {
     double ww = (i == 0) ? weights_x[i] : (weights_x[i] - weights_x[i-1]);
     if (!NumericVector::is_na(fy_xi_plus_less[i])) {
       theta += ww * (1 - fy_xi_plus_less[i]);
     }
   }

   // --- PART 2: Influence Function for X ---
   // For Losses: IF_x = [1 - F_y_less(x_i + lambda)] - theta
   NumericVector h1x(nx);
   for(int i = 0; i < nx; ++i) {
     h1x[i] = (1.0 - fy_xi_plus_less[i]) - theta;
   }

   // --- PART 3: Influence Function for Y ---
   // For Losses: IF_y = F_x(y_j - lambda) - theta
   NumericVector h2y(ny);
   for(int j = 0; j < ny; ++j) {
     h2y[j] = fx_yj_minus[j] - theta;
   }

   // Return as a list to match the "Wins" function structure
   return List::create(
     Named("theta") = theta,
     Named("IF_x")  = h1x,
     Named("IF_y")  = h2y
   );
 }
