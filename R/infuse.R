#' Prepare Data for Win Statistics Calculations
#'
#' @description
#' A highly optimized preprocessing utility that "infuses" the analysis by parsing outcomes,
#' partitioning data into treatment/control groups, and pre-calculating favorable
#' and unfavorable point estimates alongside their corresponding influence functions.
#' It utilizes \code{data.table} for fast, in-place memory modifications and produces
#' a structured S3 object optimized for downstream statistical functions.
#'
#' @param formu A formula object (e.g., \code{Score(5, "high") ~ Group}).
#'   Outcomes on the LHS can be wrapped in function-like calls to specify
#'   clinical thresholds (\code{lambda}) and direction (\code{direction}).
#'   Supports multiple outcomes (joined by \code{+}) on the LHS and exactly
#'   one binary grouping variable on the RHS.
#' @param data A data frame or \code{data.table} containing the variables specified in \code{formu}.
#' @param treated A character string specifying the level of the grouping variable to be
#'   treated as the "Treated" (Group X) arm. Defaults to \code{"Treated"}.
#'
#' @return An S3 object of class \code{"infuse"}. A list containing three stacked \code{data.table}s
#' and a \code{groups} attribute:
#' \itemize{
#'   \item \code{data}: Patient-level data combined across all outcomes. Contains:
#'     \itemize{
#'       \item \code{name}: The name of the outcome variable (factor).
#'       \item \code{type}: Detected data type (e.g., "binary", "numeric", "survival").
#'       \item \code{lambda}: Threshold defining the margin of clinical importance.
#'       \item \code{direction}: Direction of superiority ("high" or "low").
#'       \item \code{p_i_d}: Unique participant ID within their arm.
#'       \item \code{arm}: Group assignment (factor: "X" for treated, "Y" for control).
#'       \item \code{time}: The observed value or time-to-event.
#'       \item \code{status}: Status indicator (1 = observed, 0 = censored/missing).
#'       \item \code{time_ext}: Generalized Pareto Distribution (GPD) extended time.
#'       \item \code{status_ext}: GPD extension status indicator.
#'       \item \code{f_pe, u_pe}: The estimated favorable (A) and unfavorable (B) point estimates.
#'       \item \code{f_ifval, u_ifval}: Individual patient-level influence function values for favorable and unfavorable estimates.
#'     }
#'   \item \code{metadata}: Variable-level configurations. Contains:
#'     \itemize{
#'       \item \code{name}: The name of the outcome variable.
#'       \item \code{comment}: Internal metadata regarding default settings applied during parsing.
#'     }
#'   \item \code{tail}: Generalized Pareto Distribution (GPD) tail extension parameters. Contains:
#'     \itemize{
#'       \item \code{arm}: Treatment arm.
#'       \item \code{mlesi, mlesh, extra_n}: Maximum likelihood parameters for the GPD fit.
#'       \item \code{name}: The corresponding outcome variable.
#'     }
#' }
#'
#' @details
#' \bold{Performance Optimization:} This function utilizes \code{data.table} to perform
#' in-place calculations (\code{:=}) and fast row-binding (\code{rbindlist}). This eliminates
#' unnecessary memory copies, making the point estimate and influence function calculations
#' extremely fast.
#'
#' \bold{Outcome Specification:} The function supports four main types:
#' \code{numeric}, \code{binary}, \code{ordered}, and \code{survival}.
#'
#' \bold{In-Formula Arguments:} Thresholds and directions are defined directly
#' in the formula:
#' \itemize{
#'   \item \code{Binary}: Must specify the "superior" level, e.g., \code{Cured("Yes")}.
#'   \item \code{Numeric}: Requires \code{lambda} and optionally \code{direction}
#'         (defaults to \code{"high"}). E.g., \code{Weight(500, "high")}.
#'   \item \code{Survival}: Uses \code{Surv(time, status)(lambda, direction)}.
#' }
#'
#' @examples
#' # ---------------------------------------------------------
#' # Example 1: Continuous Data with Thresholds in Formula
#' # ---------------------------------------------------------
#' data(obstetrics)
#' # Treated patients "win" if birthweight is 500g higher
#' infused_obs <- infuse(Birthweight(500, high) ~ Group,
#'                     data = obstetrics,
#'                     treated = "T")
#'
#' # ---------------------------------------------------------
#' # Example 2: Multiple Outcomes (Binary & Continuous)
#' # ---------------------------------------------------------
#' # For PretermYN, "Yes" is superior. For Birthweight, 500g margin is used.
#' infused_multi <- infuse(PretermYN("Yes") + Birthweight(500) ~ Group,
#'                       data = obstetrics,
#'                       treated = "T")
#'
#' @import data.table
#' @importFrom survival Surv
#' @export
infuse <- function(formu, data, treated = "Treated") {

  DD <- prepareData(formu, data, treated )

  # 2. Vectorized Arm Logic: Compute arm1 and arm2 for the entire dataset at once
  DD$data[, arm1 := fifelse(direction == "high", as.character(arm), fifelse(arm == "X", "Y", "X"))]
  DD$data[, arm2 := fifelse(arm1 == "X", "Y", "X")]

  # 3. Grouped Calculations: Apply the custom functions by endpoint 'name'
  DD$data[, c("f_pe", "f_ifval", "u_pe", "u_ifval") := {

    # Extract scalar properties for the current group (since they are constant per 'name')
    lam <- lambda[1]
    typ <- type[1]

    ttt <- time
    ttt[status_ext==1] <- time_ext[status_ext==1]
    sss <- status
    sss[status_ext==1] <- 1

    # -- FAVORABLE --
    A <- c_fav(arm1, ttt, sss, lambda_val = lam)

    if (typ == "survival") {
      #if_A <- c_IF_fav_censoring(arm1, time, status, lambda_val = lam)$IF
      if_A <- c_IF_fav_censoring(arm1, ttt, sss, lambda_val = lam)$IF
    } else {
      if_A <- c_IF_fav_no_censoring(arm1, time, lambda_val = lam)$IF
    }

    # -- UNFAVORABLE --
    B <- c_fav(arm2, ttt, sss, lambda_val = lam)

    if (typ == "survival") {
      #if_B <- c_IF_fav_censoring(arm2, time, status, lambda_val = lam)$IF
      if_B <- c_IF_fav_censoring(arm2, ttt, sss, lambda_val = lam)$IF
    } else {
      if_B <- c_IF_fav_no_censoring(arm2, time, lambda_val = lam)$IF
    }

    # Return as a list. data.table will automatically map these to the 4 new columns.
    # It will safely recycle scalars (A and B) to the length of vectors (if_A, if_B).
    list(A, if_A, B, if_B)

  }, by = name]

  DD$data[,`:=`(arm1=NULL,arm2=NULL)]

  arm_var <- all.vars(formu[[3]])[1]
  groups_cmt <- paste0("Grouping variable: ",arm_var,";\nGroups: X = ",treated,"; Y = ",setdiff(data[[arm_var]], treated))


  structure(
    list(
      data=DD$data,
      metadata = DD$metadata,
      tail = DD$tail
    ),
    class = "infuse",
    groups = groups_cmt
  )
}
