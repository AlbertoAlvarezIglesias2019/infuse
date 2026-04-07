#' Anti-Epileptic Drug Trial Data (Progabide)
#'
#' @description
#' Data from a randomized clinical trial investigating the effect of the
#' anti-epileptic drug Progabide as an adjuvant therapy for patients
#' suffering from partial seizures.
#'
#' @details
#' In this trial, 59 patients were randomized to receive either Progabide
#' or a placebo in addition to their standard chemotherapy.
#'
#' \itemize{
#'   \item \bold{Baseline:} A baseline seizure count was recorded for the
#'   8-week period prior to randomization.
#'   \item \bold{Treatment Phase:} The number of seizures was recorded
#'   during four consecutive two-week periods following randomization.
#' }
#'
#' This dataset is in a "long" format, where each patient has four
#' corresponding rows representing each post-randomization period.
#'
#' @format A data frame with 236 rows and 6 columns:
#' \describe{
#'   \item{subject}{Unique patient identifier.}
#'   \item{treatment}{Treatment allocation: \code{"Progabide"} or \code{"Placebo"}.}
#'   \item{base}{The baseline seizure count recorded during the 8-week pre-trial period.}
#'   \item{age}{Age of the participant in years.}
#'   \item{period}{The observation period number (\code{1}, \code{2}, \code{3}, or \code{4}).}
#'   \item{seizure.rate}{The number of seizures suffered during the specific two-week period.}
#' }
#'
#' @source Thall, P.F. and Vail, S.C. (1990). Some covariance models for
#'   longitudinal count data with overdispersion. \emph{Biometrics}, 46, 657-671.
#' @keywords datasets
"epilepsy"
