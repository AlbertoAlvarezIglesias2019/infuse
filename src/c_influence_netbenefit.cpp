#include <Rcpp.h>
using namespace Rcpp;

// --- Forward Declarations ---
NumericVector c_ecdf_values(NumericVector uni,NumericVector w, NumericVector z,double ofset=0);
NumericVector c_ecdf_values_less(NumericVector uni,NumericVector w, NumericVector z,double ofset=0);


//' Calculate Influence Functions for Net Benefit (Win Difference)
//'
//' @description
//' This function computes the point estimate for the Net Benefit (Win Difference)
//' and returns the individual-level influence functions for groups X and Y.
//' These components are used to derive standard errors and confidence intervals
//' for hierarchical or weighted outcomes.
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
//'    \item \code{theta}: The point estimate for Net Benefit \eqn{\theta = P(X > Y + \lambda) - P(Y > X + \lambda)}.
//'    \item \code{IF_x}: Numeric vector of influence functions for subjects in group X.
//'    \item \code{IF_y}: Numeric vector of influence functions for subjects in group Y.
//' }
//'
//' @details
//' The influence function for group X is defined as \eqn{IF_x(x) = h_1(x) - \theta},
//' and for group Y as \eqn{IF_y(y) = h_2(y) - \theta}.
//' The point estimate is calculated by integrating the weighted ECDFs using
//' the masses derived from \code{x_w} and \code{y_w}.
//'
//' @export
// [[Rcpp::export]]
List c_influence_netbenefit(NumericVector x, NumericVector y,
                            NumericVector x_uni, NumericVector y_uni,
                            NumericVector x_w, NumericVector y_w,
                   double lambda) {

  int nx = x.size();
  int ny = y.size();

  // --- PART 1: Point Estimate Calculation ---

  // F_Y(x_i - lambda)
  NumericVector fy_xi_minus = c_ecdf_values(y_uni, y_w, x, lambda);
  // Calculate weights for each x_i to perform the integration
  NumericVector weights_x = c_ecdf_values(x_uni, x_w, x);
  double A = 0.0;
  for(int i = 0; i < nx; ++i) {
    double ww = (i == 0) ? weights_x[i] : (weights_x[i] - weights_x[i-1]);
    if (!NumericVector::is_na(fy_xi_minus[i])) {
      A += ww * fy_xi_minus[i];
    }

  }

  // F_X(y_j - lambda)
  NumericVector fx_yj_minus = c_ecdf_values(x_uni, x_w, y, lambda);
  // Calculate weights for each y_j to perform the integration
  NumericVector weights_y = c_ecdf_values(y_uni, y_w, y);
  double B = 0.0;
  for(int j = 0; j < ny; ++j) {
    double ww = (j == 0) ? weights_y[j] : (weights_y[j] - weights_y[j-1]);
    if (!NumericVector::is_na(fx_yj_minus[j])) {
      B += ww * fx_yj_minus[j];
    }

  }

  double T_FnX_FmY = A - B;

  // --- PART 2: Influence Functions ---

  // For Group X: h_1_x = F_Y(x_i - lambda) - 1 + F_Y_less(x_i + lambda)
  NumericVector fy_xi_plus_less = c_ecdf_values_less(y_uni, y_w, x, -lambda);
  NumericVector IF_x(nx);

  for(int i = 0; i < nx; ++i) {
    // The IF is defined as h(x) - theta
    IF_x[i] = (fy_xi_minus[i] - 1.0 + fy_xi_plus_less[i]) - T_FnX_FmY;
  }

  // For Group Y: h_2_y = 1 - F_X_less(y_j + lambda) - F_X(y_j - lambda)
  NumericVector fx_yj_plus_less = c_ecdf_values_less(x_uni, x_w, y, -lambda);
  NumericVector IF_y(ny);

  for(int j = 0; j < ny; ++j) {
    // The IF is defined as h(y) - theta
    IF_y[j] = (1.0 - fx_yj_plus_less[j] - fx_yj_minus[j]) - T_FnX_FmY;
  }

  // Return as a List
  return List::create(
    Named("theta") = T_FnX_FmY,
    Named("IF_x")  = IF_x,
    Named("IF_y")  = IF_y
  );
}

