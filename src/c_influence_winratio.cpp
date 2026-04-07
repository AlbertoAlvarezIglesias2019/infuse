#include <Rcpp.h>
using namespace Rcpp;

// --- Forward Declarations ---
NumericVector c_ecdf_values(NumericVector uni, NumericVector w, NumericVector z, double ofset=0);
NumericVector c_ecdf_values_less(NumericVector uni, NumericVector w, NumericVector z, double ofset=0);

//' Calculate Influence Functions for Win Ratio
//'
//' @description
//' This function computes the point estimate for the Win Ratio and returns
//' the individual-level influence functions for groups X and Y. These
//' components are calculated on the log-scale to maintain numerical symmetry
//' for standard error derivation.
//'
//' @param x NumericVector: Raw values for the first group (X).
//' @param y NumericVector: Raw values for the second group (Y).
//' @param x_uni NumericVector: Unique sorted values of group X.
//' @param y_uni NumericVector: Unique sorted values of group Y.
//' @param x_w NumericVector: Cumulative ECDF weights for unique X values.
//' @param y_w NumericVector: Cumulative ECDF weights for unique Y values.
//' @param lambda double: The margin of clinical significance.
//'
//' @return List: A named list containing:
//' \itemize{
//'    \item \code{theta}: The point estimate for the Log Win Ratio \eqn{\theta = \log(A) - \log(B)}.
//'    \item \code{IF_x}: Numeric vector of influence functions for group X (log-scale).
//'    \item \code{IF_y}: Numeric vector of influence functions for group Y (log-scale).
//' }
//'
//' @export
//'
// [[Rcpp::export]]
 List c_influence_winratio(NumericVector x, NumericVector y,
                           NumericVector x_uni, NumericVector y_uni,
                           NumericVector x_w, NumericVector y_w,
                           double lambda) {

   int nx = x.size();
   int ny = y.size();

   // --- PART 1: Point Estimate Components ---

   // A = P(X > Y + lambda)
   NumericVector fy_xi_minus = c_ecdf_values(y_uni, y_w, x, lambda);
   NumericVector weights_x = c_ecdf_values(x_uni, x_w, x);
   double A = 0.0;
   for(int i = 0; i < nx; ++i) {
     double ww = (i == 0) ? weights_x[i] : (weights_x[i] - weights_x[i-1]);
     if (!NumericVector::is_na(fy_xi_minus[i])) A += ww * fy_xi_minus[i];
   }

   // B = P(Y > X + lambda)
   NumericVector fx_yj_minus = c_ecdf_values(x_uni, x_w, y, lambda);
   NumericVector weights_y = c_ecdf_values(y_uni, y_w, y);
   double B = 0.0;
   for(int j = 0; j < ny; ++j) {
     double ww = (j == 0) ? weights_y[j] : (weights_y[j] - weights_y[j-1]);
     if (!NumericVector::is_na(fx_yj_minus[j])) B += ww * fx_yj_minus[j];
   }

   // Calculate theta as log(A) - log(B)
   double log_theta = std::log(A) - std::log(B);

   // --- PART 2: Influence Functions (Log-scale) ---

   // IF for X: (F_Y(x-L) / A) - ((1 - F_Y_less(x+L)) / B)
   NumericVector fy_xi_plus_less = c_ecdf_values_less(y_uni, y_w, x, -lambda);
   NumericVector IF_x(nx);
   for(int i = 0; i < nx; ++i) {
     IF_x[i] = (fy_xi_minus[i] / A) - ((1.0 - fy_xi_plus_less[i]) / B);
   }

   // IF for Y: ((1 - F_X_less(y+L)) / A) - (F_X(y-L) / B)
   NumericVector fx_yj_plus_less = c_ecdf_values_less(x_uni, x_w, y, -lambda);
   NumericVector IF_y(ny);
   for(int j = 0; j < ny; ++j) {
     IF_y[j] = ((1.0 - fx_yj_plus_less[j]) / A) - (fx_yj_minus[j] / B);
   }

   return List::create(
     Named("theta") = log_theta,
     Named("IF_x")  = IF_x,
     Named("IF_y")  = IF_y
   );
 }
