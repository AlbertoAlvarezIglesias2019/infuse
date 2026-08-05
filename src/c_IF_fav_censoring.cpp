#include <Rcpp.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// Forward declarations of your existing functions
List c_ecdf_plus(NumericVector x, Nullable<IntegerVector> xs = R_NilValue);
NumericVector c_ecdf_predict_less(NumericVector z, NumericVector x_sorted, NumericVector p_sorted);

//' Calculate Subject-Level Influence Functions for Favorable Outcome Estimand
//'
//' @description
//' Computes subject-level influence function (IF) values for the favorable outcome
//' estimand \eqn{A = P(X \ge Y + \lambda)} under right-censoring. The function handles
//' missing data (\code{NA} values) safely and maps computed influence values back
//' to the original row ordering of the input dataset.
//'
//' @param arm A \code{CharacterVector} specifying treatment allocation (\code{"X"} or \code{"Y"}).
//' @param time A \code{NumericVector} of observed event or censoring times.
//' @param status An \code{IntegerVector} of event status indicators (\code{1} for event, \code{0} for censored).
//' @param lambda_val A \code{double} specifying the margin/threshold shift (\eqn{\lambda}) applied to arm Y outcomes.
//'
//' @details
//' The influence function calculations account for Kaplan-Meier/ECDF estimation under right-censoring:
//' \enumerate{
//'   \item \bold{NA Handling & Indexing:} Identifies complete cases per arm (\code{idx_X}, \code{idx_Y}) while preserving original row positions. Returns an \code{NA}-filled column if either group lacks valid data.
//'   \item \bold{Cross-Arm Survival Prediction:} Evaluates step-function survival probabilities \eqn{S_X(Y + \lambda)} and \eqn{S_Y(X - \lambda)} using \code{\link{c_ecdf_predict_less}}.
//'   \item \bold{Arm X Influence Calculation (\eqn{\text{IF}_X}):} Evaluates the martingale-based influence components accounting for the probability at risk \eqn{H_X(t)} and cumulative hazard increments \eqn{d\Lambda_X(t)}.
//'   \item \bold{Arm Y Influence Calculation (\eqn{\text{IF}_Y}):} Computes corresponding martingale integral terms for control group observations.
//'   \item \bold{Row Reconstruction:} Re-aligns computed \eqn{\text{IF}_X} and \eqn{\text{IF}_Y} values into an \eqn{N}-length vector, placing \code{NA_REAL} in rows where inputs were incomplete.
//' }
//'
//' @return A \code{DataFrame} containing four columns matching the original length \eqn{N}:
//' \item{arm}{Original \code{arm} vector.}
//' \item{time}{Original \code{time} vector.}
//' \item{status}{Original \code{status} vector.}
//' \item{IF}{Computed subject-level influence function value (or \code{NA} for incomplete cases).}
//'
//' @name c_IF_fav_censoring
//' @export
// [[Rcpp::export]]
 DataFrame c_IF_fav_censoring(CharacterVector arm, NumericVector time, IntegerVector status, double lambda_val) {
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

   // Safeguard: If we don't have valid data in both arms, return early with NA values
   if (nx == 0 || ny == 0) {
     NumericVector IF_out(N, NA_REAL);
     return DataFrame::create(
       Named("arm") = arm,
       Named("time") = time,
       Named("status") = status,
       Named("IF") = IF_out,
       Named("stringsAsFactors") = false
     );
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

   // 3. Call your ECDF functions (now guaranteed to have no NAs)
   List ppx = c_ecdf_plus(time_X, status_X);
   List ppy = c_ecdf_plus(time_Y, status_Y);

   NumericVector p_X = ppx["p"];
   NumericVector w_X = ppx["w"];
   IntegerVector o_X = ppx["o"];

   NumericVector p_Y = ppy["p"];
   NumericVector w_Y = ppy["w"];
   IntegerVector o_Y = ppy["o"];

   // 4. Create sorted subsets for the predict step
   NumericVector x_s(nx), px_s(nx);
   for(int i = 0; i < nx; i++) {
     x_s[i] = time_X[o_X[i] - 1];
     px_s[i] = p_X[o_X[i] - 1];
   }

   NumericVector y_s(ny), py_s(ny);
   for(int i = 0; i < ny; i++) {
     y_s[i] = time_Y[o_Y[i] - 1];
     py_s[i] = p_Y[o_Y[i] - 1];
   }

   // 5. Predictions
   NumericVector dy_time_plus_lambda(ny);
   for(int i = 0; i < ny; i++) dy_time_plus_lambda[i] = time_Y[i] + lambda_val;

   NumericVector pred_less_X = c_ecdf_predict_less(dy_time_plus_lambda, x_s, px_s);
   NumericVector S_X_ge(ny);
   for(int i = 0; i < ny; i++) S_X_ge[i] = 1.0 - pred_less_X[i];

   NumericVector dx_time_minus_lambda(nx);
   for(int i = 0; i < nx; i++) dx_time_minus_lambda[i] = time_X[i] - lambda_val;

   NumericVector pred_less_Y = c_ecdf_predict_less(dx_time_minus_lambda, y_s, py_s);
   NumericVector S_Y(nx);
   for(int i = 0; i < nx; i++) S_Y[i] = 1.0 - pred_less_Y[i];


   // ==========================================
   // IF_A_X Calculation
   // ==========================================
   NumericVector phi_x(nx), H_X(nx), dLambda_X(nx), IF_X(nx);

   // Pass 1: Compute components
   for(int k = 0; k < nx; k++) {
     double x_i = time_X[k];
     double sum_phi = 0.0;

     for(int j = 0; j < ny; j++) {
       if(time_Y[j] >= std::max(0.0, x_i - lambda_val)) {
         sum_phi += w_Y[j] * S_X_ge[j];
       }
     }
     phi_x[k] = sum_phi;

     int count_ge = 0;
     for(int m = 0; m < nx; m++) {
       if(time_X[m] >= x_i) count_ge++;
     }
     H_X[k] = (double)count_ge / nx;
     dLambda_X[k] = status_X[k] / (nx * H_X[k]);
   }

   // Pass 2: Calculate IF
   for(int k = 0; k < nx; k++) {
     double x_start = time_X[k];
     double term1 = -status_X[k] * phi_x[k] / H_X[k];
     double term2 = 0.0;

     for(int m = 0; m < nx; m++) {
       if(time_X[m] <= x_start) {
         term2 += (phi_x[m] / H_X[m]) * dLambda_X[m];
       }
     }
     IF_X[k] = term1 + term2;
   }


   // ==========================================
   // IF_A_Y Calculation
   // ==========================================
   NumericVector psi_y(ny), H_Y(ny), dLambda_Y(ny), IF_Y(ny);

   // Pass 1: Compute components
   for(int k = 0; k < ny; k++) {
     double y_j = time_Y[k];
     double sum_psi = 0.0;

     for(int i = 0; i < nx; i++) {
       if(time_X[i] >= y_j + lambda_val) {
         sum_psi += w_X[i] * S_Y[i];
       }
     }
     psi_y[k] = sum_psi;

     int count_ge = 0;
     for(int m = 0; m < ny; m++) {
       if(time_Y[m] >= y_j) count_ge++;
     }
     H_Y[k] = (double)count_ge / ny;
     dLambda_Y[k] = status_Y[k] / (ny * H_Y[k]);
   }

   // Pass 2: Calculate IF
   for(int k = 0; k < ny; k++) {
     double y_start = time_Y[k];
     double term1 = status_Y[k] * psi_y[k] / H_Y[k];
     double term2 = 0.0;

     for(int m = 0; m < ny; m++) {
       if(time_Y[m] <= y_start) {
         term2 -= (psi_y[m] / H_Y[m]) * dLambda_Y[m];
       }
     }
     IF_Y[k] = term1 + term2;
   }

   // ==========================================
   // Recombine into original layout (Preserves NAs)
   // ==========================================
   NumericVector IF_out(N, NA_REAL); // Initialize all N positions to NA

   // Map calculated X values back to original rows
   for(int i = 0; i < nx; i++) {
     IF_out[idx_X[i]] = IF_X[i];
   }
   // Map calculated Y values back to original rows
   for(int i = 0; i < ny; i++) {
     IF_out[idx_Y[i]] = IF_Y[i];
   }

   // Return the final DataFrame matching input sizes
   return DataFrame::create(
     Named("arm") = arm,
     Named("time") = time,
     Named("status") = status,
     Named("IF") = IF_out,
     Named("stringsAsFactors") = false
   );
 }

