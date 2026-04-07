#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

//' Extend the tail of a distribution using Generalized Pareto Distribution (GPD)
//' 
//' @description This function takes existing observation vectors and extends 
//' them by generating extra values based on GPD parameters. It combines the 
//' original data with newly generated tail data.
//' 
//' @param uni A NumericVector of existing observation values.
//' @param w A NumericVector of existing probabilities/weights.
//' @param mlesi Double; the scale parameter for the GPD.
//' @param mlesh Double; the shape parameter for the GPD.
//' @param extra_n Double; the number of extra points to generate for the tail.
//' 
//' @return A named list containing:
//' \itemize{
//'   \item \code{uni}: A combined numeric vector of original and extended values.
//'   \item \code{w}: A combined numeric vector of associated probabilities/weights.
//' }
//' 
//' 
//' @examples
//' \dontrun{
//' data(ipilimumab)
//' dat <- ipilimumab
//' 
//' ### Read the data
//' X <- dat %>% filter(arm=="ipilimumab") %>% pull(time)
//' Xs <- dat %>% filter(arm=="ipilimumab") %>% pull(event)
//' 
//' ### Remove the last few events to illustrate the functionality
//' wher <- X>26
//' Xs[wher] <- 0
//' 
//' ### Estimate and plot the Kaplan-Meier estimation of the ECDF
//' fit <- c_ecdf_surv(X,Xs)
//' plot(w~uni,data=fit,type="s",ylim=c(0,1))
//' abline(h=1,lty=3)
//' 
//' ### Estimate the scale and shape parameters for the GPD
//' p <- param_gpd(X,Xs)
//' 
//' ### Extend the tail
//' uni <- fit$uni
//' w <- fit$w
//' mlesi <- p$mlesi
//' mlesh <- p$mlesh
//' ta <- c_extendtail(uni,w,mlesi,mlesh,p$extra_n)
//' ta
//' 
//' ### plot the extension
//' plot(w~uni,data=ta,type="s",col="red",ylim=c(0,1),xlim = c(0,max(ta$uni)))
//' points(w~uni,data=fit,type="s")
//' abline(h=1,lty=3)
//'
//' ###############
//' ### Example 2
//' ###############
//' ### Read the data
//' X <- dat %>% filter(arm=="placebo") %>% pull(time)
//' Xs <- dat %>% filter(arm=="placebo") %>% pull(event)
//' 
//' ### Remove the last few events to illustrate the functionality
//' wher <- X>10
//' Xs[wher] <- 0
//' 
//' ### Estimate aplan-Meier estimation of the ECDF
//' fit <- c_ecdf_surv(X,Xs)
//' 
//' ### Estimate the scale and shape parameters for the GPD
//' p <- param_gpd(X,Xs)
//' 
//' ### Extend the tail
//' ta <- c_extendtail(fit$uni,fit$w,p$mlesi,p$mlesh,extra_n=10)
//' ta
//' 
//' ### plot the extension
//' plot(w~uni,data=ta,type="s",col="red",ylim=c(0,1),xlim = c(0,max(ta$uni)))
//' points(w~uni,data=fit,type="s")
//' abline(h=1,lty=3)
//' 
//' ### Extend the tail but only with 5 values
//' ta <- c_extendtail(fit$uni,fit$w,p$mlesi,p$mlesh,extra_n=5)
//' ta
//' 
//' ### plot the extension
//' plot(w~uni,data=ta,type="s",col="red",ylim=c(0,1),xlim = c(0,max(ta$uni)))
//' points(w~uni,data=fit,type="s")
//' abline(h=1,lty=3)
//' }
//' @export
// [[Rcpp::export]]
 List c_extendtail(NumericVector uni,
                            NumericVector w, 
                            double mlesi,
                            double mlesh,
                            double extra_n) {
   
   int n = uni.size();
   int i_extra_n = (int)extra_n;
   
   // --- Internal Helpers ---
   auto internal_pmax = [](double x1, double x2) {
     return (x1 >= x2) ? x1 : x2;
   };
   
   auto internal_qgpd = [&](double p, double loc, double scale, double shape) {
     double out;
     p = 1.0 - p; 
     // lambda is 0 based on original code defaults
     if (shape == 0) {
       out = loc - scale * std::log(p);
     } else {
       out = loc + scale * (std::pow(p, -shape) - 1.0) / shape;
     }
     return out;
   };
   
   auto internal_pgpd = [&](double q, double loc, double scale, double shape) {
     double out;
     double q_std = internal_pmax(q - loc, 0.0) / scale;
     if (shape == 0) {
       out = 1.0 - std::exp(-q_std);
     } else {
       out = internal_pmax(1.0 + shape * q_std, 0.0);
       out = 1.0 - std::pow(out, -1.0 / shape);
     }
     return out;
   };
   
   // --- Logic ---
   
   // If input is empty, return empty list
   if (n == 0 && i_extra_n == 0) {
     return List::create(Named("uni") = NumericVector(0),
                         Named("w")   = NumericVector(0));
   }
   
   
   // Maximum observed time (the threshold v)
   //double v = *std::max_element(uni.begin(), uni.end());
   double v = (n > 0) ? *std::max_element(uni.begin(), uni.end()) : 0.0;
   
   // Estimated survival prob at threshold v
   //double max_w = *std::max_element(w.begin(), w.end());
   double max_w = (n > 0) ? *std::max_element(w.begin(), w.end()) : 0.0;
   double survInV = 1.0 - max_w;
   
   //NumericVector extra_uni(i_extra_n);
   //NumericVector extra_w(i_extra_n);
   std::vector<double> extra_uni_vec;
   std::vector<double> extra_w_vec;
   extra_uni_vec.reserve(i_extra_n);
   extra_w_vec.reserve(i_extra_n);
   
   //for (int i = 1; i <= i_extra_n; ++i) {
  //   // Generate quantile values for the tail extension
  //   double p_val = (double)i / extra_n - 1.0 / (extra_n * 2.0);
  //   extra_uni[i - 1] = internal_qgpd(p_val, v, mlesi, mlesh);
  // }
   
  // for (int i = 0; i < i_extra_n; ++i) {
  //   // Generate probability values for the tail extension
  //   extra_w[i] = 1.0 - survInV * (1.0 - internal_pgpd(extra_uni[i], v, mlesi, mlesh));
  // }
   
   for (int i = 1; i <= i_extra_n; ++i) {
     // Generate quantile values for the tail extension
     double p_val = (double)i / extra_n - 1.0 / (extra_n * 2.0);
     double q_val = internal_qgpd(p_val, v, mlesi, mlesh);
     
     extra_uni_vec.push_back(q_val);
     
     // Generate probability values for the tail extension
     double w_val = 1.0 - survInV * (1.0 - internal_pgpd(q_val, v, mlesi, mlesh));
     extra_w_vec.push_back(w_val);
   }
   
   
   // Ensure the tail closes at 1
   if (i_extra_n > 0) {
     extra_w_vec[i_extra_n - 1] = 1.0;
   }
   
   
//   // Prepare output: combined length (n + extra_n) * 2
//   int total_len = n + i_extra_n;
//   NumericVector out(total_len * 2);
//   
//   // Fill original and extra values for observations
//   for (int i = 0; i < n; ++i) {
//     out[i] = uni[i];
//   }
//   for (int i = 0; i < i_extra_n; ++i) {
//     out[n + i] = extra_uni[i];
//   }
//   
//   // Fill original and extra values for weights (offset by total_len)
//   for (int i = 0; i < n; ++i) {
//     out[total_len + i] = w[i];
//   }
//   for (int i = 0; i < i_extra_n; ++i) {
//     out[total_len + n + i] = extra_w[i];
//   }
   
//   return out;

  // Combine original and extra data
  int total_len = n + i_extra_n;
   NumericVector final_uni(total_len);
   NumericVector final_w(total_len);
   
   // Fill final vectors
   for (int i = 0; i < n; ++i) {
     final_uni[i] = uni[i];
     final_w[i] = w[i];
   }
   for (int i = 0; i < i_extra_n; ++i) {
     final_uni[n + i] = extra_uni_vec[i];
     final_w[n + i] = extra_w_vec[i];
   }
   
   // Return as a named R list
   return List::create(
     Named("uni") = final_uni,
     Named("w")   = final_w
   );
 }