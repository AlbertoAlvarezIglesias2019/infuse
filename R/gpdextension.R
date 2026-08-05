#' Tail Extension Using Generalized Pareto Distribution (GPD)
#'
#' @description
#' This function performs tail extrapolation for survival/time-to-event data.
#' It estimates GPD parameters using Maximum Likelihood Estimation (MLE) on
#' exceedances past an automatically calculated threshold (the 80th percentile
#' of observed events). It then mathematically projects new event times for
#' censored observations that occur after the last observed event time, while
#' perfectly preserving the original input ordering of the data.
#'
#' @param x A numeric vector representing time-to-event or survival times.
#' @param xs A numeric or integer vector representing censoring status
#'   (1 = observed event, 0 = censored).
#'
#' @return A named list containing:
#' \item{mlesi}{Estimated GPD scale parameter at the baseline threshold `u`.}
#' \item{mlesh}{Estimated GPD shape parameter.}
#' \item{u}{The calculated baseline threshold (80th percentile of observed events).}
#' \item{extra_n}{The number of censored observations that were extrapolated.}
#' \item{x}{The output time vector, matching the original input order, with extrapolated values injected.}
#' \item{xs}{The output status vector, matching the original input order, with updated event statuses (1 for extrapolated points).}
#' \item{xtail}{A binary flag vector matching the original input order (1 if the row was extrapolated, 0 otherwise).}
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
gpdextension <- function(x, xs) {

  or <- order(x,-xs)

  ox <- x[or]
  oxs <- xs[or]

  # --- 1. Threshold Selection ---
  # u is the 0.8 quantile of observed events
  event_times <- ox[oxs > 0]
  if (length(event_times) == 0) stop("No events observed; cannot estimate GPD parameters.")

  u <- sort(event_times)[floor(length(event_times) * 0.8)]

  # --- 2. Data Preparation ---
  # Subset data to exceedances (the tail)
  is_tail <- ox > u
  is_tail[is.na(is_tail)]<-FALSE
  x_tail  <- ox[is_tail] - u  # Shift to start at 0
  xs_tail <- oxs[is_tail]

  # --- 3. Maximum Likelihood Estimation ---
  # Initial guesses for scale (0.1) and shape (0.1)
  mle_fit <- optim(
    par     = c(0.1, 0.1),
    fn      = llik,
    times   = x_tail[!is.na(x_tail)],
    status  = xs_tail[!is.na(xs_tail)],
    method  = "Nelder-Mead",
    control = list(maxit = 10000)
  )

  mlesi <- mle_fit$par[1]
  mlesh <- mle_fit$par[2]


  # --- 4. Tail Extension Calculation (extra_n) ---
  max_event_time <- max(ox[oxs == 1], na.rm = TRUE)
  wher <- ox >= max_event_time & oxs==0
  wher[is.na(wher)]<- FALSE

  extra_n <- sum(wher,na.rm=TRUE)
  oxtail <- as.numeric(wher)
  if (extra_n>0) {
    t_start <- ox[wher][1]
    u_star <- t_start
    mlesi_start <- mlesi + mlesh*(u_star-u)
    ppp <- 1:extra_n/(extra_n+1)
    if (!mlesh==0) {newox <- u_star + mlesi_start * ((1-ppp)^(-mlesh) - 1)/mlesh} else {newox <- u_star -mlesi_start*log(1-ppp)}
    ox[wher] <- newox
    oxs[wher] <- 1
    }

  op <- order(or)

  return(list(
    mlesi   = mlesi,
    mlesh   = mlesh,
    u       = u,
    extra_n = extra_n,
    x = ox[op],
    xs = oxs[op],
    xtail = oxtail[op]
  ))
}
