#' Respiratory Illness Clinical Trial Data
#'
#' @description
#' Data from a clinical trial involving 111 patients with respiratory illness
#' from two different centers. Patients were randomized to receive either a
#' placebo or an active treatment, with respiratory status assessed at
#' baseline and over four subsequent visits.
#'
#' @details
#' This is a longitudinal dataset where each patient is represented by 4 rows
#' (one for each post-baseline visit). The primary outcome is the respiratory
#' status, a binary categorical variable.
#'
#' \itemize{
#'   \item \bold{Treatment Groups:} Patients received either Active treatment
#'   (\code{"A"}) or Placebo (\code{"P"}).
#'   \item \bold{Assessment:} Respiratory status was determined by clinicians
#'   and categorized as \code{1} (Good) or \code{0} (Poor).
#' }
#'
#' @format A data frame with 444 rows and 8 columns:
#' \describe{
#'   \item{center}{The clinical center involved in the study (coded as \code{1} or \code{2}).}
#'   \item{id}{Unique patient identifier.}
#'   \item{treat}{Treatment allocation: \code{"A"} for Active treatment and \code{"P"} for Placebo.}
#'   \item{sex}{Factor indicating the sex of the participant (\code{"M"} or \code{"F"}).}
#'   \item{age}{Age of the participant in years.}
#'   \item{baseline}{Respiratory status at the start of the study: \code{1} = Good, \code{0} = Poor.}
#'   \item{visit}{The visit number (\code{1}, \code{2}, \code{3}, or \code{4}).}
#'   \item{outcome}{The primary outcome: Respiratory status at the specific visit: \code{1} = Good, \code{0} = Poor.}
#' }
#'
#' @source Clinical trial data often used in longitudinal and categorical data analysis.
#' @keywords datasets
"respiratory"
