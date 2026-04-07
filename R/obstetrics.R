#' Periodontal Treatment and Birth Outcomes (Obstetrics)
#'
#' @description
#' Data from a randomized controlled trial investigating whether treating
#' maternal periodontal disease reduces the risk of pre-term birth and
#' low birth weight.
#'
#' @details
#' Periodontal disease is an inflammatory condition characterized by the
#' destruction of tissue and/or bone around the teeth. A major component of
#' the disease is oral colonization by gram-negative bacteria; systemic release
#' of cytokines and/or lipopolysaccharides from these bacteria may impact
#' fetal condition.
#'
#' In this trial:
#' \itemize{
#'   \item \bold{Treatment Group (T):} Received periodontal treatment, oral
#'   hygiene instruction, and tooth polishing at follow-ups.
#'   \item \bold{Control Group (C):} Underwent only brief oral examinations.
#' }
#'
#' Several observational studies have suggested an association between maternal
#' periodontal disease and pre-term birth, which this RCT aimed to validate.
#'
#' @format A data frame with 814 rows and 6 columns:
#' \describe{
#'   \item{Pid}{Unique participant identifier.}
#'   \item{Clinic}{The clinical site where the patient was treated (KY, MN, MS, NY).}
#'   \item{Group}{Treatment allocation: \code{"T"} for Treated and \code{"C"} for Controls.}
#'   \item{Age}{Age of the participant in years.}
#'   \item{PretermYN}{Binary indicator for pre-term birth (Yes/No).}
#'   \item{Birthweight}{The birth weight of the infant measured in grams.}
#' }
#'
#' @source Randomized Controlled Trial data.
#' @keywords datasets
"obstetrics"
