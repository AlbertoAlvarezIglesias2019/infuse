#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;


// Forward Declaration
IntegerVector c_order(NumericVector x);

//' Kaplan-Meier Product-Limit Estimator (C++ Implementation)
//'
//' @description Computes the Kaplan-Meier estimate of the cumulative incidence 
//' (1 - S(t)) for right-censored data and returns a structured list. This 
//' implementation is optimized for speed and correctly handles tied event times.
//'
//' @param ttt A numeric vector of observed times (time-to-event or time-to-censoring).
//' @param sss An integer vector of status indicators (1 = event/death, 0 = censored).
//'
//' @return A named list containing two components:
//' \itemize{
//'   \item \code{uni}: A numeric vector of the unique event times (sorted).
//'   \item \code{w}: A numeric vector of the cumulative incidence (1 - S(t)) 
//'   calculated at those times.
//' }
//'
//' @details
//' The function calculates survival probability using the formula:
//' \deqn{\hat{S}(t) = \prod_{i: t_i \le t} \left(1 - \frac{d_i}{n_i}\right)}
//' where \eqn{d_i} is the number of events at time \eqn{t_i} and \eqn{n_i}
//' is the number of subjects at risk just before \eqn{t_i}.
//'
//' Note: If an observation is censored at the same time an event occurs, the
//' censored observation is assumed to remain at risk for that event (standard
//' Kaplan-Meier convention).
//'
//' @examples
//' \dontrun{
//' # Example data: times and event status
//' times <- c(10, 20, 20, 35, 40, 50)
//' status <- c(1, 1, 0, 1, 0, 1)
//'
//' km_res <- c_ecdf_surv(times, status)
//'
//' # Plotting the cumulative incidence step function
//' plot(km_res$uni, km_res$w, type = "s",
//'      xlab = "Time", ylab = "Cumulative Incidence",
//'      main = "KM Estimate")
//'
//' # Comparison with survival package
//' x <- c(0.2, 1, 2, 4, 2, 2, 35, 6, 3, 2, 3, 1, 3, 1.5, 45, 23, 12)
//' xs <- c(1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0)
//' 
//' library(survival)
//' fit <- survfit(Surv(x, xs) ~ 1)
//' fit1 <- c_ecdf_surv(x, xs)
//'
//' data.frame(time1 = summary(fit)$time,
//'            time2 = fit1$uni,
//'            surv1 = 1 - summary(fit)$surv,
//'            surv2 = fit1$w)
//'
//' # Test with ipilimumab data
//' data(ipilimumab)
//' dat <- ipilimumab
//' res_ipili <- c_ecdf_surv(dat$time, dat$event)
//' fit_ipili <- survfit(Surv(time, event) ~ 1, data = dat)
//' 
//' # Check equivalence
//' all.equal(res_ipili$uni, summary(fit_ipili)$time)
//' all.equal(res_ipili$w, 1 - summary(fit_ipili)$surv)
//'
//' # ++++++ TEST SPEED
//' start.time <- Sys.time()
//' fit_bench <- survfit(Surv(x, xs) ~ 1)
//' end.time <- Sys.time()
//' time.taken1 <- as.numeric(end.time - start.time)
//'
//' start.time <- Sys.time()
//' tf_bench <- c_ecdf_surv(x, xs)
//' end.time <- Sys.time()
//' time.taken2 <- as.numeric(end.time - start.time)
//'
//' time.taken1
//' time.taken2
//' time.taken1 / time.taken2
//' }
//'
//' @export
// [[Rcpp::export]]
List c_ecdf_surv(NumericVector ttt, IntegerVector sss) {
   //int n = ttt.size();
   //if (n == 0) return NumericVector(0);

   int n_orig = ttt.size();
   //if (n_orig == 0) return NumericVector(0);
   if (n_orig == 0) {
     return List::create(Named("uni") = NumericVector(0),
                         Named("w")   = NumericVector(0));
   }
   
   // 1. Filter out observations with NA in time (ttt)
   // We keep only observations where ttt is not NA.
   std::vector<double> clean_t;
   std::vector<int> clean_s;
   clean_t.reserve(n_orig);
   clean_s.reserve(n_orig);

   for (int i = 0; i < n_orig; ++i) {
     if (!NumericVector::is_na(ttt[i]) && !IntegerVector::is_na(sss[i])) {
       clean_t.push_back(ttt[i]);
       clean_s.push_back(sss[i]);
     }
   }

   int n = clean_t.size();
   //if (n == 0) return NumericVector(0);
   if (n == 0) {
     return List::create(Named("uni") = NumericVector(0),
                         Named("w")   = NumericVector(0));
   }
   
   // Convert back to Rcpp types for use with your existing c_order function
   NumericVector t_vec = wrap(clean_t);
   IntegerVector s_vec = wrap(clean_s);


   // 2. Get ordering indices (Based on the cleaned data)
   IntegerVector ord = c_order(t_vec);

   std::vector<double> unique_times;
   std::vector<double> survival_probs;

   double current_surv = 1.0;
   int i = 0;

   while (i < n) {
     double current_time = t_vec[ord[i] - 1];
     int events = 0;
     int n_at_risk = n - i;

     // Handle tied time points
     while (i < n && t_vec[ord[i] - 1] == current_time) {
       if (s_vec[ord[i] - 1] > 0) {
         events++;
       }
       i++;
     }

     // Kaplan-Meier Product-Limit Calculation
     if (events > 0) {
       current_surv *= (1.0 - static_cast<double>(events) / n_at_risk);
       unique_times.push_back(current_time);
       survival_probs.push_back(1.0 - current_surv);
     }
   }

   // 3. Prepare output (Times in first half, 1-Surv in second half)
   //int k = unique_times.size();
   //NumericVector out(2 * k);
   //for (int j = 0; j < k; ++j) {
   //   out[j] = unique_times[j];
   //   out[j + k] = 1.0 - survival_probs[j];
   // }

   //return out;
   
   // 3. Return as a named R list
   return List::create(
     Named("uni") = wrap(unique_times),
     Named("w")   = wrap(survival_probs)
   );
 }
