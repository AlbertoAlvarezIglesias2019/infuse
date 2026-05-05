#' infuse: Win Statistics Inference via Influence Functions
#'
#' @description
#' The \code{infuse} package provides a unified framework for estimating Net
#' Treatment Benefit (NTB) and Win Ratio (WR) across continuous, binary,
#' count, and time-to-event outcomes. By leveraging Empirical Survival
#' Functions and Influence Functions, \code{infuse} offers significant speed
#' advantages simplifies the analysis of win statistics.
#'
#' @section Getting Started:
#' The best way to learn \code{infuse} is through the Quick Start Guide:
#' \code{vignette("infuse_start_guide", package = "infuse")}
#'
#' @section Core Workflow:
#' \enumerate{
#'   \item \code{\link{infuse}}: Initialize data and pre-calculate distributions.
#'   \item \code{\link{realise}}: Extract statistical inference (CI, p-values, NPO).
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib infuse, .registration = TRUE
## usethis namespace: end
NULL


