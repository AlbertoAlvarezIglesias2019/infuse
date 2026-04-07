#' Log-Likelihood for Generalized Pareto Distribution with Censoring
#'
#' @description
#' Calculates the negative log-likelihood of a Generalized Pareto Distribution (GPD)
#' given survival times and censoring indicators. This is used internally for
#' maximum likelihood estimation of the tail extension parameters.
#'
#' @param pa A numeric vector of length 2: \code{pa[1]} is the scale parameter (sigma),
#'   and \code{pa[2]} is the shape parameter (xi).
#' @param times A numeric vector of observed times (exceedances).
#' @param status A numeric vector of censoring indicators (1 = observed, 0 = censored).
#'
#' @details
#' The likelihood handles censoring by using the PDF for observed events and the
#' Survival Function (1 - CDF) for censored observations.
#'
#' Constraints enforced:
#' \itemize{
#'   \item Scale (\eqn{\sigma}) > 0
#'   \item Shape (\eqn{\xi}) between -0.5 and 0.2
#'   \item Support: \eqn{1 + \xi \cdot \frac{t}{\sigma} > 0}
#' }
#'
#' @return The negative log-likelihood value (for minimization).
#' @noRd
llik <- function(pa, times, status) {

  sigma <- pa[1] # scale
  xi    <- pa[2] # shape

  # --- 1. Constraint Validation ---
  # Check parameter bounds and GPD support constraints
  # xi < -0.5: Usually too "thin-tailed" for stable estimation
  # xi > 0.2: Usually too "heavy-tailed" for clinical survival data
  if (sigma <= 0 || xi < -0.5 || xi > 0.2) {
    return(1e7) # Return large value to nudge optimizer away
  }

  # Support check: 1 + xi * t / sigma must be positive
  support <- 1 + (xi * times / sigma)
  if (any(support <= 0)) {
    return(1e7)
  }

  # --- 2. Likelihood Calculation ---
  # log(PDF) for status == 1
  # log(Survival) for status == 0
  log_pdf      <- -log(sigma) - (1/xi + 1) * log(support)
  log_survival <- -(1/xi) * log(support)

  nll <- -sum(status * log_pdf + (1 - status) * log_survival)

  # Safety check for non-finite results (Inf/NaN)
  if (!is.finite(nll)) return(1e7)

  return(nll)
}
