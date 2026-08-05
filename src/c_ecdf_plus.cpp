#include <Rcpp.h>
#include <numeric>
#include <algorithm>

using namespace Rcpp;

//' Flexible Empirical Cumulative Distribution Function (Unified CDF)
//'
//' @description A high-performance unified function to calculate the Empirical
//' Cumulative Distribution Function (ECDF) for any numerical variable. It
//' automatically generalizes to the Kaplan-Meier estimator when right-censoring
//' indicators are provided, making it suitable for both general statistics
//' and time-to-event analysis.
//'
//' @param x A numeric vector of observations (e.g., ages, heights, or survival times). NAs are ignored in calculations but preserved in output.
//' @param xs An optional integer vector of indicators. For survival data:
//' (1 = event, 0 = censored). For standard ECDF: leave NULL or provide all 1s.
//'
//' @return A named \code{List} containing:
//' \itemize{
//'   \item \code{p}: The cumulative probabilities \eqn{F(x)}. In the absence of
//'   censoring, this matches the standard ECDF.
//'   \item \code{w}: The probability mass or "weight" assigned to each observation.
//'   \item \code{o}: The 1-based ordering index used during the internal calculation.
//'   \item \code{pi}: Numeric vector mapping the proportion of the non-NA cohort
//'         remaining at risk (\eqn{\hat{\pi}(t)}) at that individual's time point.
//'   \item \code{lambda_jump}: Numeric vector containing the hazard jump (\eqn{d\hat{\Lambda}(t)}),
//'         equivalent to the status indicator divided by the current size of the risk set.
//' }
//'
//' @details
//' When \code{xs} is NULL, this function computes the standard ECDF by assigning
//' equal weights (\eqn{1/n}) to all observations. When \code{xs} is provided,
//' it applies the "Redistribute-to-the-Right" algorithm.
//'
//' To ensure consistency with survival analysis theory, ties in \code{x} are
//' resolved by placing events before censoring. The output is always mapped
//' back to the original input order of \code{x}.
//'
//' @examples
//'
//' x <- c(1, 2, 2, 3, 5)
//' c_ecdf_plus(x)
//'
//' # Example with more data
//' set.seed(123)
//' y <- rnorm(100)
//' res_y <- c_ecdf_plus(y)
//' or <- res_y$o
//' plot(y[or], res_y$p[or], type = "s", main = "C++ ECDF")
//'
//'
//' data_points <- c(5, 1, 3, 1, 5, 5)
//' result <- c_ecdf_plus(data_points)
//'
//' # Access cumulative weights
//' print(result$p[result$o])
//'
//' # compare to ecdf
//' x <- c(1,2,4,2,2,35,6,3,2,3,1,3)
//' c_ecdf_plus(x)
//' unique(x)
//' ecdf(x)(unique(x))
//'
//'
//' times <- c(10, 20, 20, 30)
//' status <- c(1, 0, 1, 1) # Note the tie at 20
//'
//' # Calculate KM stats
//' res <- c_ecdf_plus(times, status)
//'
//' # The weights for censored items will be 0
//' res$w
//'
//' # Example data: times and event status
//' times <- c(10, 20, 20, 35, 40, 50)
//' status <- c(1, 1, 0, 1, 0, 1)
//'
//' km_res <- c_ecdf_plus(times, status)
//'
//' # Plotting the cumulative incidence step function
//' plot(times, km_res$p[km_res$o], type = "s",
//'      xlab = "Time", ylab = "Cumulative Incidence",
//'      main = "KM Estimate")
//'
//'
//' # Test with ipilimumab data
//' data(ipilimumab)
//' dat <- ipilimumab
//' res_ipili <- c_ecdf_plus(dat$time, dat$event)
//' fit_ipili <- survfit(Surv(time, event) ~ 1, data = dat)
//'
//' # Check equivalence
//' plot(summary(fit_ipili)$time,1-summary(fit_ipili)$surv,type="s")
//'  or <- res_ipili$o
//' points(dat$time[or],res_ipili$p[or],col="red",type="s")
//'
//'
//'
//' @export
// [[Rcpp::export]]
List c_ecdf_plus(NumericVector x, Nullable<IntegerVector> xs = R_NilValue) {
  int n = x.size();
  IntegerVector delta(n);

  if (xs.isNull()) {
    delta = rep(1, n);
  } else {
    delta = as<IntegerVector>(xs);
  }

  // 1. Create index vector
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);

  // 2. Sorting logic:
  // - Put NAs at the very end
  // - Primary sort by x ascending
  // - Secondary sort by delta descending (events before censoring)
  std::sort(idx.begin(), idx.end(), [&](int i, int j) {
    bool na_i = NumericVector::is_na(x[i]);
    bool na_j = NumericVector::is_na(x[j]);
    if (na_i && !na_j) return false; // i is NA, j is not -> i comes after j
    if (!na_i && na_j) return true;  // i is not NA, j is -> i comes before j
    if (na_i && na_j) return i < j;  // both NA, keep original order

    if (x[i] != x[j]) return x[i] < x[j];
    return delta[i] > delta[j];
  });

  // 3. Count non-NA values
  int n_non_na = 0;
  for (int i = 0; i < n; ++i) {
    if (!NumericVector::is_na(x[idx[i]])) {
      n_non_na++;
    } else {
      break; // Since they are sorted, all subsequent are NA
    }
  }

  NumericVector p_sorted(n, NA_REAL);
  NumericVector w_sorted(n, NA_REAL);
  NumericVector pi_sorted(n, NA_REAL);          // Track sorted proportion at risk
  NumericVector lambda_jump_sorted(n, NA_REAL); // Track sorted hazard jumps
  IntegerVector o(n);

  // 4. Calculate only for non-NA values
  double current_S = 1.0;
  for (int i = 0; i < n_non_na; ++i) {
    int original_idx = idx[i];
    // Risk set only counts remaining non-NA values
    double risk_set = static_cast<double>(n_non_na - i);
    double hazard_ratio = static_cast<double>(delta[original_idx]) / risk_set;

    w_sorted[i] = hazard_ratio * current_S;
    current_S *= (1.0 - hazard_ratio);
    p_sorted[i] = 1.0 - current_S;
    o[i] = original_idx + 1;

    // Calculate new metrics
    // (Using n_non_na as the denominator since NAs are excluded from the risk set entirely)
    pi_sorted[i] = risk_set / static_cast<double>(n_non_na);
    lambda_jump_sorted[i] = hazard_ratio;
  }

  // Handle NA indices in the 'o' vector
  for (int i = n_non_na; i < n; ++i) {
    o[i] = idx[i] + 1;
  }

  // 5. Map results back to the original input order
  NumericVector p_final(n);
  NumericVector w_final(n);
  NumericVector pi_final(n);
  NumericVector lambda_jump_final(n);
  for (int i = 0; i < n; ++i) {
    if (i < n_non_na) {
      p_final[idx[i]] = p_sorted[i];
      w_final[idx[i]] = w_sorted[i];
      pi_final[idx[i]] = pi_sorted[i];
      lambda_jump_final[idx[i]] = lambda_jump_sorted[i];
    } else {
      p_final[idx[i]] = NA_REAL;
      w_final[idx[i]] = NA_REAL;
      pi_final[idx[i]] = NA_REAL;
      lambda_jump_final[idx[i]] = NA_REAL;
    }
  }

  return List::create(
    Named("p") = p_final,
    Named("w") = w_final,
    Named("o") = o,
    Named("pi") = pi_final,
    Named("lambda_jump") = lambda_jump_final
  );
}
