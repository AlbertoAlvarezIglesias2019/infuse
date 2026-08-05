#include <Rcpp.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// Forward declarations of your existing functions so the compiler can link them
List c_ecdf_plus(NumericVector x, Nullable<IntegerVector> xs = R_NilValue);
//NumericVector c_ecdf_predict(NumericVector z, NumericVector x_sorted, NumericVector p_sorted);
NumericVector c_ecdf_predict_less(NumericVector z, NumericVector x_sorted, NumericVector p_sorted);
NumericVector c_ecdf_predict(NumericVector z, NumericVector x_sorted, NumericVector p_sorted);


//' Calculate Uncensored Subject-Level Influence Functions for Favorable Outcome Estimand
//'
//' @description
//' Computes subject-level influence function (IF) values for the favorable outcome
//' estimand \eqn{A = P(X \ge Y + \lambda)} when outcomes are fully observed (uncensored).
//' Automatically handles \code{NA} values and restores computed influence values back to
//' the original row order of the input vectors.
//'
//' @param arm A \code{CharacterVector} specifying treatment allocation (\code{"X"} or \code{"Y"}).
//' @param time A \code{NumericVector} of fully observed outcomes/times (no status variable required).
//' @param lambda_val A \code{double} specifying the margin/threshold shift (\eqn{\lambda}) applied to arm Y outcomes.
//'
//' @details
//' In the absence of censoring, the empirical influence functions simplify to direct ECDF projections:
//' \enumerate{
//'   \item \bold{NA-Safe Filtering:} Identifies valid cases for arm X (\code{idx_X}) and arm Y (\code{idx_Y}). If either group has no complete cases, returns an \code{NA}-filled output table.
//'   \item \bold{Empirical CDF Estimation:} Fits unweighted ECDFs using \code{\link{c_ecdf_plus}} without status indicators.
//'   \item \bold{Cross-Arm Predictions:} Computes \eqn{S_X(Y_j + \lambda) = 1 - F_X(Y_j + \lambda)} for arm Y subjects and \eqn{F_Y(X_i - \lambda)} for arm X subjects using \code{\link{c_ecdf_predict_less}}.
//'   \item \bold{Point Estimate (\eqn{A}):} Calculates \eqn{A = \frac{1}{n_Y} \sum_{j=1}^{n_Y} S_X(Y_j + \lambda)}.
//'   \item \bold{Exact IF Evaluation:}
//'     \itemize{
//'       \item \bold{Arm X subjects:} \eqn{\text{IF}_{X, i} = F_Y(X_i - \lambda) - A}
//'       \item \bold{Arm Y subjects:} \eqn{\text{IF}_{Y, j} = S_X(Y_j + \lambda) - A}
//'     }
//'   \item \bold{Row Realignment:} Places calculated \eqn{\text{IF}} values into an \eqn{N}-length vector, setting missing/incomplete rows to \code{NA_REAL}.
//' }
//'
//' @return A \code{DataFrame} containing three columns matching the original length \eqn{N}:
//' \item{arm}{Original \code{arm} vector.}
//' \item{time}{Original \code{time} vector.}
//' \item{IF}{Computed subject-level influence function value (or \code{NA} for incomplete cases).}
//'
//' @name c_IF_fav_no_censoring
//' @export
// [[Rcpp::export]]
 DataFrame c_IF_fav_no_censoring(CharacterVector arm, NumericVector time, double lambda_val) {
   int N = arm.size();

   // 1. Separate indices for complete cases only to maintain original row ordering
   std::vector<int> idx_X, idx_Y;
   for(int i = 0; i < N; i++) {
     // Skip row if arm or time is NA
     if (CharacterVector::is_na(arm[i]) ||
         NumericVector::is_na(time[i])) {
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
       Named("IF") = IF_out,
       Named("stringsAsFactors") = false
     );
   }

   // 2. Extract complete-case subsets natively in Rcpp
   NumericVector time_X(nx), time_Y(ny);

   for(int i = 0; i < nx; i++) {
     time_X[i] = time[idx_X[i]];
   }
   for(int i = 0; i < ny; i++) {
     time_Y[i] = time[idx_Y[i]];
   }

   // 3. Call your ECDF functions without status
   // Because xs is Nullable and defaults to R_NilValue, c_ecdf_plus uses standard ECDF weights
   List ppx = c_ecdf_plus(time_X);
   List ppy = c_ecdf_plus(time_Y);

   NumericVector p_X = ppx["p"];
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

   // 5. Predict S_X_ge
   NumericVector dy_time_plus_lambda(ny);
   for(int i = 0; i < ny; i++) dy_time_plus_lambda[i] = time_Y[i] + lambda_val;

   NumericVector pred_less_X = c_ecdf_predict_less(dy_time_plus_lambda, x_s, px_s);
   NumericVector S_X_ge(ny);
   for(int i = 0; i < ny; i++) S_X_ge[i] = 1.0 - pred_less_X[i];

   // 6. Predict F_Y (Using c_ecdf_predict)
   NumericVector dx_time_minus_lambda(nx);
   for(int i = 0; i < nx; i++) dx_time_minus_lambda[i] = time_X[i] - lambda_val;

   NumericVector F_Y = c_ecdf_predict(dx_time_minus_lambda, y_s, py_s);

   // 7. Calculate A
   double A = 0.0;
   for(int j = 0; j < ny; j++) {
     A += S_X_ge[j] * w_Y[j];
   }

   // ==========================================
   // Recombine into original layout (Preserves NAs)
   // ==========================================
   NumericVector IF_out(N, NA_REAL); // Initialize all positions to NA

   // Map X values back to original rows (F_Y - A)
   for(int i = 0; i < nx; i++) {
     IF_out[idx_X[i]] = F_Y[i] - A;
   }
   // Map Y values back to original rows (S_X_ge - A)
   for(int i = 0; i < ny; i++) {
     IF_out[idx_Y[i]] = S_X_ge[i] - A;
   }

   // Return the final DataFrame matching original sizes (omitting status)
   return DataFrame::create(
     Named("arm") = arm,
     Named("time") = time,
     Named("IF") = IF_out,
     Named("stringsAsFactors") = false
   );
 }
