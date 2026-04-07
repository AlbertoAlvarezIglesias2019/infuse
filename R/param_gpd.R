#' Maximum Likelihood Estimation for GPD Parameters
#'
#' @description
#' Estimates the scale and shape parameters of the Generalised Pareto Distribution (GPD)
#' using maximum likelihood. This function specifically handles right-censored
#' survival data to model the "tail" of the distribution beyond a specific threshold.
#'
#' @param x A numeric vector of time-to-event values.
#' @param xs A numeric vector of censoring indicators (1 = event, 0 = censored).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{mlesi}: Maximum Likelihood Estimate of the \bold{scale} parameter (\eqn{\sigma}).
#'   \item \code{mlesh}: Maximum Likelihood Estimate of the \bold{shape} parameter (\eqn{\xi}).
#'   \item \code{u}: The threshold used for exceedance calculation.
#'   \item \code{extra_n}: The suggested number of points to append to the tail for extension.
#' }
#'
#' @details
#' \bold{Threshold Selection:} The threshold \code{u} is set as the 80th percentile
#' of observed event times (\code{xs == 1}). Only values exceeding this threshold
#' are used to fit the GPD.
#'
#' \bold{Optimization:} Uses the Nelder-Mead algorithm via \code{optim} to
#' minimize the negative log-likelihood produced by \code{llik}.
#'
#' \bold{Tail Extension (\code{extra_n}):} If the survival curve does not reach
#' zero (\code{max(w) < 1}), \code{extra_n} calculates the number of points
#' required to complete the tail based on the proportion of observations
#' at or beyond the maximum observed event time.
#'
#' @examples
#' library(survival)
#' library(dplyr)
#' data(ipilimumab)
#'
#' # Extract treated arm data
#' x <- ipilimumab %>% filter(arm == "ipilimumab") %>% pull(time)
#' xs <- ipilimumab %>% filter(arm == "ipilimumab") %>% pull(event)
#'
#' # Standard estimation
#' param_gpd(x, xs)
#'
#' @export
param_gpd <- function(x, xs) {

  # --- 1. Threshold Selection ---
  # u is the 0.8 quantile of observed events
  event_times <- x[xs > 0]
  if (length(event_times) == 0) stop("No events observed; cannot estimate GPD parameters.")

  u <- sort(event_times)[floor(length(event_times) * 0.8)]

  # --- 2. Data Preparation ---
  # Subset data to exceedances (the tail)
  is_tail <- x > u
  x_tail  <- x[is_tail] - u  # Shift to start at 0
  xs_tail <- xs[is_tail]

  # --- 3. Maximum Likelihood Estimation ---
  # Initial guesses for scale (0.1) and shape (0.1)
  mle_fit <- optim(
    par     = c(0.1, 0.1),
    fn      = llik,
    times   = x_tail,
    status  = xs_tail,
    method  = "Nelder-Mead",
    control = list(maxit = 10000)
  )

  mlesi <- mle_fit$par[1]
  mlesh <- mle_fit$par[2]

  # --- 4. Tail Extension Calculation (extra_n) ---
  fit <- c_ecdf_surv(x, xs)
  surv_at_max <- 1 - max(fit$w)

  if (surv_at_max == 0) {
    extra_n <- 0
  } else {
    # Find max observed event time
    max_event_time <- max(x[xs == 1], na.rm = TRUE)
    # Calculate extra_n based on proportion of points at/above max event
    prop_at_end <- mean(x >= max_event_time, na.rm = TRUE)
    extra_n     <- round(prop_at_end * sum(!is.na(x)), 0)
  }

  return(list(
    mlesi   = mlesi,
    mlesh   = mlesh,
    u       = u,
    extra_n = extra_n
  ))
}
