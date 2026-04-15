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
#'   \item \code{\link{sow}}: Initialize data and pre-calculate distributions.
#'   \item \code{\link{harvest}}: Calculate point estimates and influence functions.
#'   \item \code{\link{reap}}: Extract statistical inference (CI, p-values, NPO).
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
#' @useDynLib infuse, .registration = TRUE
## usethis namespace: end
NULL


