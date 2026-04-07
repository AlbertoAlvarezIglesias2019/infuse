#' Progression-Free Survival with Ipilimumab in Prostate Cancer
#'
#' @description
#' Data from a randomized clinical trial assessing the efficacy of ipilimumab
#' following radiotherapy in patients with metastatic castration-resistant
#' prostate cancer that had progressed after docetaxel chemotherapy.
#'
#' @details
#' Ipilimumab is a fully human monoclonal antibody that binds cytotoxic
#' T-lymphocyte antigen 4 (CTLA-4) to enhance antitumour immunity.
#'
#' In this study, participants were randomized (1:1) to receive bone-directed
#' radiotherapy (8 Gy in one fraction) followed by either:
#' \itemize{
#'   \item \bold{Ipilimumab:} 10 mg/kg every 3 weeks for up to four doses.
#'   \item \bold{Placebo:} Administered on the same schedule as the active arm.
#' }
#'
#' The primary outcome for this specific dataset is Progression-Free Survival (PFS).
#'
#' @format A data frame with 799 rows and 3 columns:
#' \describe{
#'   \item{arm}{Treatment allocation: \code{"ipilimumab"} or \code{"placebo"}.}
#'   \item{event}{Censoring indicator: \code{1} = progression/event, \code{0} = censored.}
#'   \item{time}{Time to progression or censoring, measured in months.}
#' }
#'
#' @source Clinical trial data involving metastatic castration-resistant prostate cancer.
#' @keywords datasets
"ipilimumab"
