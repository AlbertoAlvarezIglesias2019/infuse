#include <Rcpp.h>
using namespace Rcpp;

// --- Forward Declarations ---
NumericVector c_ecdf_values(NumericVector uni,NumericVector w, NumericVector z,double ofset=0);
NumericVector c_ecdf_values_less(NumericVector uni,NumericVector w, NumericVector z,double ofset=0);

//' Calculate Win Probability Point Estimate and Influence Functions
//'
//' This function computes the point estimate for the probability of a favorable
//' pair (Win Probability, P(X >= Y + lambda)) and returns the influence functions
//' for both the treatment (X) and control (Y) groups. These influence functions
//' can be used for variance estimation and combining multiple outcomes.
//'
//' @param x NumericVector: Raw values of the treatment group.
//' @param y NumericVector: Raw values of the control group.
//' @param x_uni NumericVector: Sorted unique values of the treatment group (X).
//' @param y_uni NumericVector: Sorted unique values of the control group (Y).
//' @param x_w NumericVector: Cumulative ECDF weights for x_uni.
//' @param y_w NumericVector: Cumulative ECDF weights for y_uni.
//' @param lambda double: The clinically meaningful threshold for defining a win.
//'
//' @return List: A list containing:
//'    \item{theta}{Point Estimate: P(X >= Y + lambda)}
//'    \item{IF_x}{Influence Function values for the treatment group}
//'    \item{IF_y}{Influence Function values for the control group}
//'
//' @details
//' The win probability theta is calculated by integrating the ECDF of Y against
//' the empirical distribution of X. The influence functions are calculated as:
//' IF_x = F_Y(x - lambda) - theta and IF_y = [1 - F_X_less(y + lambda)] - theta.
//'
//' @export
// [[Rcpp::export]]
List c_influence_favorable(NumericVector x, NumericVector y,
                   NumericVector x_uni, NumericVector y_uni,
                   NumericVector x_w, NumericVector y_w,
                   double lambda) {

   int nx = x.size();
   int ny = y.size();

   // --- PART 1: Point Estimate Calculation ---
   // Evaluate ECDF of Y at each (x_i - lambda)
   //NumericVector fy_xi_minus = c_ecdf_values(y_uni, y_w, x, lambda);
   NumericVector fx_yj_plus_less = c_ecdf_values_less(x_uni, x_w, y, -lambda);
   
   // Calculate weights for each x_i to perform the integration
   //NumericVector weights_x = c_ecdf_values(x_uni, x_w, x);
   NumericVector weights_y = c_ecdf_values(y_uni, y_w, y);
   double theta = 0.0;

   //for(int i = 0; i < nx; ++i) {
   for(int j = 0; j < ny; ++j) {
     double ww = (j == 0) ? weights_y[j] : (weights_y[j] - weights_y[j-1]);
     if (!NumericVector::is_na(fx_yj_plus_less[j])) {
       theta += ww * (1 - fx_yj_plus_less[j]);
     }
   }

   
   // --- PART 2: Influence Function for X ---
   // IF_x(x_i) = F_Y(x_i - lambda) - theta
   NumericVector fy_xi_minus = c_ecdf_values(y_uni, y_w, x, lambda);
   NumericVector h1x(nx);
   for(int i = 0; i < nx; ++i) {
     h1x[i] = fy_xi_minus[i] - theta;
   }

   // --- PART 3: Influence Function for Y ---
   // IF_y(y_j) = [1 - F_X_less(y_j + lambda)] - theta
   // This represents the contribution of control observations to the variance
   //NumericVector fx_yj_plus_less = c_ecdf_values_less(x_uni, x_w, y, -lambda);
   NumericVector h2y(ny);
   for(int j = 0; j < ny; ++j) {
     h2y[j] = (1.0 - fx_yj_plus_less[j]) - theta;
   }

   // Return only the estimate and the IFs
   return List::create(
     Named("theta") = theta,
     Named("IF_x")  = h1x,
     Named("IF_y")  = h2y
   );
}
