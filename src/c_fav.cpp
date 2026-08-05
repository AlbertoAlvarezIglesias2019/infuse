#include <Rcpp.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// Forward declarations of your existing functions so the compiler can link them
List c_ecdf_plus(NumericVector x, Nullable<IntegerVector> xs = R_NilValue);
NumericVector c_ecdf_predict_less(NumericVector z, NumericVector x_sorted, NumericVector p_sorted);

//' Calculate Favorable Outcome Estimand A with Censoring and Margin Shift
//'
//' @description
//' Computes the population-level favorable outcome probability \eqn{A = P(X \ge Y + \lambda)}
//' between two treatment arms ("X" and "Y"). It handles continuous or time-to-event
//' data with optional right-censoring via non-parametric survival/ECDF estimation.
//'
//' @param arm A \code{CharacterVector} specifying group allocation (must contain \code{"X"} or \code{"Y"}).
//' @param time A \code{NumericVector} of observed event or censoring times.
//' @param status An \code{IntegerVector} indicating event status (typically \code{1} for observed event, \code{0} for censored).
//' @param lambda_val A \code{double} specifying a margin or threshold shift (\eqn{\lambda}) applied to arm Y's outcomes.
//'
//' @details
//' The function executes the following steps:
//' \enumerate{
//'   \item \bold{Complete Case Filtering:} Automatically filters out observations where \code{arm}, \code{time}, or \code{status} contains \code{NA}.
//'   \item \bold{Group-wise Non-parametric Estimation:} Fits Kaplan-Meier / modified ECDF step functions using \code{\link{c_ecdf_plus}} separately for arm X and arm Y.
//'   \item \bold{Predictive Evaluation:} Evaluates the survival probability \eqn{S_X(t) = P(X \ge t)} at shifted time points \eqn{t_Y + \lambda} via \code{\link{c_ecdf_predict_less}}.
//'   \item \bold{Weighted Integration:} Computes the overall estimand \eqn{A = \sum S_X(Y_j + \lambda) \cdot w_{Y, j}} using the survival jumps/weights (\eqn{w_Y}).
//' }
//'
//' @return A \code{double} representing the favorable outcome probability \eqn{A}.
//' Returns \code{NA_REAL} if either treatment group contains zero complete cases.
//'
//' @export
// [[Rcpp::export]]
 double c_fav(CharacterVector arm, NumericVector time, IntegerVector status, double lambda_val) {
   int N = arm.size();

   // 1. Separate indices for complete cases only to maintain original row ordering
   std::vector<int> idx_X, idx_Y;
   for(int i = 0; i < N; i++) {
     // Skip row if arm, time, or status is NA
     if (CharacterVector::is_na(arm[i]) ||
         NumericVector::is_na(time[i]) ||
         IntegerVector::is_na(status[i])) {
       continue;
     }

     if (arm[i] == "X") {
       idx_X.push_back(i);
     } else if (arm[i] == "Y") {
       idx_Y.push_back(i);
     }
   }

   int nx = idx_X.size();
   int ny = idx_Y.size();

   // Safeguard: If we don't have valid data in both arms, return NA
   if (nx == 0 || ny == 0) {
     return NA_REAL;
   }

   // 2. Extract complete-case subsets natively in Rcpp
   NumericVector time_X(nx), time_Y(ny);
   IntegerVector status_X(nx), status_Y(ny);

   for(int i = 0; i < nx; i++) {
     time_X[i] = time[idx_X[i]];
     status_X[i] = status[idx_X[i]];
   }
   for(int i = 0; i < ny; i++) {
     time_Y[i] = time[idx_Y[i]];
     status_Y[i] = status[idx_Y[i]];
   }

   // 3. Call your ECDF functions with status
   List ppx = c_ecdf_plus(time_X, status_X);
   List ppy = c_ecdf_plus(time_Y, status_Y);

   NumericVector p_X = ppx["p"];
   IntegerVector o_X = ppx["o"];

   NumericVector p_Y = ppy["p"];
   NumericVector w_Y = ppy["w"];
   IntegerVector o_Y = ppy["o"];

   // 4. Create sorted subsets for the predict step (correcting R's 1-based index)
   NumericVector x_s(nx), px_s(nx);
   for(int i = 0; i < nx; i++) {
     x_s[i] = time_X[o_X[i] - 1];
     px_s[i] = p_X[o_X[i] - 1];
   }

   // 5. Predict S_X_ge
   NumericVector dy_time_plus_lambda(ny);
   for(int i = 0; i < ny; i++) {
     dy_time_plus_lambda[i] = time_Y[i] + lambda_val;
   }

   NumericVector pred_less_X = c_ecdf_predict_less(dy_time_plus_lambda, x_s, px_s);
   NumericVector S_X_ge(ny);
   for(int i = 0; i < ny; i++) {
     S_X_ge[i] = 1.0 - pred_less_X[i];
   }

   // 6. Calculate A (Weighted sum using Kapen-Meier/ECDF weights)
   double A = 0.0;
   for(int j = 0; j < ny; j++) {
     A += S_X_ge[j] * w_Y[j];
   }

   return A;
 }
