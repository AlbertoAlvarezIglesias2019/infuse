#' Prepare Data for Win/Loss Probability Calculations
#'
#' @description
#' A preprocessing utility that "sows" the analysis by parsing outcomes
#' (continuous or survival), partitioning data into treatment/control groups,
#' and pre-calculating Empirical Cumulative Distribution Functions (ECDF).
#' It produces a structured S3 object optimized for the \code{harvest} function.
#'
#' @param formu A formula object (e.g., \code{Surv(time, status) + score ~ treatment}).
#'   Supports multiple outcomes on the LHS and exactly one binary grouping variable on the RHS.
#' @param data A data frame containing the variables specified in \code{formu}.
#' @param treated A string specifying the level of the grouping variable to be
#'   treated as the "Treated" (Group X) arm. Defaults to \code{"Treated"}.
#' @param good A character vector (\code{"high"} or \code{"low"}) indicating
#'   clinical favorability. Values are recycled if shorter than the number of outcomes.
#'
#' @return An S3 object of class \code{"sow"}. A nested list containing:
#' \itemize{
#'   \item \code{xid, yid}: Participant IDs for Treated (X) and Control (Y) groups.
#'   \item \code{Variables}: A named list for each outcome containing pre-calculated
#'   ECDFs, weights, and GPD tail parameters if applicable.
#' }
#'
#' @details
#' \bold{Data Partitioning:} The function assigns internal IDs (\code{p_i_d})
#' to maintain the joint distribution (correlation) between multiple outcomes
#' for the same individual during later resampling or influence function steps.
#'
#' \bold{Survival Tail Extension:} If a survival curve does not reach zero
#' (censoring at the end), a Generalized Pareto Distribution (GPD) is used
#' to extend the tail via \code{param_gpd} and \code{c_extendtail}.
#'
#' @examples
#' # ---------------------------------------------------------
#' # Example 1: Survival Data (Ipilimumab)
#' # ---------------------------------------------------------
#' data(ipilimumab)
#' # Higher time is "good" for Progression-Free Survival
#' sown_surv <- sow(Surv(time, event) ~ arm,
#'                  data = ipilimumab,
#'                  treated = "ipilimumab",
#'                  good = "high")
#'
#' print(sown_surv)
#' summary(sown_surv)
#' plot(sown_surv)
#'
#' # ---------------------------------------------------------
#' # Example 2: Continuous Data (Obstetrics)
#' # ---------------------------------------------------------
#' data(obstetrics)
#' # High birthweight is generally considered "good"
#' sown_obs <- sow(Birthweight ~ Group,
#'                 data = obstetrics,
#'                 treated = "T",
#'                 good = "high")
#'
#' summary(sown_obs)
#' plot(sown_obs)
#'
#' # ---------------------------------------------------------
#' # Example 3: Longitudinal Binary Data (Respiratory)
#' # ---------------------------------------------------------
#' data(respiratory)
#' # We analyze status at the final visit (Visit 4)
#' resp_v4 <- respiratory[respiratory$visit == 4, ]
#' sown_resp <- sow(outcome ~ treat,
#'                  data = resp_v4,
#'                  treated = "A",
#'                  good = "high")
#'
#' summary(sown_resp)
#' # Note: For binary outcomes, the plot shows a single-step ECDF
#' plot(sown_resp)
#'
#' # ---------------------------------------------------------
#' # Example 4: Multiple Outcomes (Epilepsy)
#' # ---------------------------------------------------------
#' data(epilepsy)
#' # Looking at seizure rate at Period 4. Lower rate is "good".
#' epi_v4 <- epilepsy[epilepsy$period == 4, ]
#' sown_epi <- sow(seizure.rate ~ treatment,
#'                 data = epi_v4,
#'                 treated = "Progabide",
#'                 good = "low")
#'
#' print(sown_epi)
#' summary(sown_epi)
#'
#' @importFrom stringr str_split str_trim str_detect str_remove_all
#' @importFrom survival Surv
#' @export
sow <- function(formu, data, treated = "Treated", good = "high") {

  data <- as.data.frame(data)

  # --- 1. Identify Grouping (RHS) ---
  arm_var <- all.vars(formu[[3]])[1]
  if (length(unique(data[[arm_var]])) != 2) {
    stop("Grouping variable must have exactly two levels.")
  }

  # Assign IDs within each arm to maintain correlation across multiple outcomes
  #data$p_i_d <- sequence(table(data[[arm_var]]))

  data$p_i_d <-0
  is_treated <- data[[arm_var]] == treated
  data[is_treated,]$p_i_d <- 1:sum(is_treated)
  data[!is_treated,]$p_i_d <- 1:sum(!is_treated)

  x_id <- data$p_i_d[is_treated]
  y_id <- data$p_i_d[!is_treated]

  # --- 2. Identify Outcomes (LHS) ---
  lhs_string <- deparse(formu[[2]])
  outcomes_raw <- stringr::str_trim(unlist(stringr::str_split(lhs_string, "\\+")))
  n_outcomes <- length(outcomes_raw)

  # Recycle 'good' argument if user provided fewer values than outcomes
  good_vec <- rep(good, length.out = n_outcomes)

  # --- 3. Process Each Outcome ---
  processed_vars <- lapply(seq_along(outcomes_raw), function(k) {

    current_out <- outcomes_raw[k]
    is_surv <- stringr::str_detect(current_out, "Surv")

    # Extract raw vectors and status
    if (is_surv) {
      out_type <- "Survival"
      tmp_f <- as.formula(paste(current_out, "~ 1"))
      surv_obj <- model.frame(tmp_f, data = data)[[1]]
      time_var  <- surv_obj[, 1]
      status_var <- surv_obj[, 2]
    } else {
      if (length(unique(na.omit(data[[current_out]]))) == 2) {out_type <- "binary"} else {out_type <- class(data[[current_out]])}
      time_var  <- data[[current_out]]
      status_var <- rep(1, length(time_var)) # 1 = observed for continuous
    }

    # Split into X (Treated) and Y (Control)
    x_val <- time_var[is_treated]
    y_val <- time_var[!is_treated]
    x_stat <- status_var[is_treated]
    y_stat <- status_var[!is_treated]

    x_stat[is.na(x_val)]<-NA
    y_stat[is.na(y_val)]<-NA

    # Handle Character/Factor inputs by converting to numeric
    if (is.character(x_val) || is.factor(x_val)) {
      x_val <- as.numeric(as.factor(x_val)) - 1
      y_val <- as.numeric(as.factor(y_val)) - 1
    }

    # --- 4. Calculate ECDFs ---
    if (out_type == "Survival") {
      ecdf_x <- c_ecdf_surv(x_val, x_stat)
      ecdf_y <- c_ecdf_surv(y_val, y_stat)
    } else {
      ecdf_x <- c_ecdf(x_val)
      ecdf_y <- c_ecdf(y_val)
    }

    # --- 5. GPD Tail Extension (Survival Only) ---
    tail_params <- list(x = list(mlesi=0, mlesh=0, extra_n=0),
                        y = list(mlesi=0, mlesh=0, extra_n=0))

    if (out_type == "Survival") {
      # Process Treated
      p_x <- param_gpd(x_val, x_stat)
      if (p_x$extra_n > 0) {
        ecdf_x <- c_extendtail(ecdf_x$uni, ecdf_x$w, p_x$mlesi, p_x$mlesh, p_x$extra_n)
        idx <- which(x_val > max(x_val[x_stat == 1], na.rm = TRUE))
        if (length(idx) > 0) x_val[idx] <- tail(ecdf_x$uni, length(idx))
      }

      # Process Control
      p_y <- param_gpd(y_val, y_stat)
      if (p_y$extra_n > 0) {
        ecdf_y <- c_extendtail(ecdf_y$uni, ecdf_y$w, p_y$mlesi, p_y$mlesh, p_y$extra_n)
        idx <- which(y_val > max(y_val[y_stat == 1], na.rm = TRUE))
        if (length(idx) > 0) y_val[idx] <- tail(ecdf_y$uni, length(idx))
      }
      tail_params$x <- p_x
      tail_params$y <- p_y
    }

    # --- 6. Internal Weight Calculation ---
    calc_weights <- function(uni, w_vec, raw_vals) {
      ord <- order(raw_vals)
      cum_w <- c_ecdf_values(uni, w_vec, raw_vals[ord])
      # Calculate differences to get individual point weights
      w <- cum_w - dplyr::lag(cum_w, default = 0)
      return(w[order(ord)])
    }

    list(
      x = x_val, xs = x_stat, x_weights = calc_weights(ecdf_x$uni, ecdf_x$w, x_val),
      y = y_val, ys = y_stat, y_weights = calc_weights(ecdf_y$uni, ecdf_y$w, y_val),
      x_uni = ecdf_x$uni, x_w = ecdf_x$w,
      y_uni = ecdf_y$uni, y_w = ecdf_y$w,
      tail_params = tail_params,
      good = good_vec[k],
      type = out_type
    )
  })

  # --- 7. Finalize Names and Object ---
  # Cleans Surv(Time, Status) into just "Time" for list names
  clean_names <- stringr::str_remove_all(outcomes_raw, "Surv\\(|\\)|\\s") |>
    stringr::str_split(",") |>
    sapply(`[[`, 1)
  names(processed_vars) <- clean_names

  structure(
    list(
      xid = x_id,
      yid = y_id,
      Variables = processed_vars
    ),
    class = "sow"
  )
}
