#' Prepare Data for Win/Loss Probability Calculations
#'
#' @description
#' A preprocessing utility that "sows" the analysis by parsing outcomes,
#' partitioning data into treatment/control groups, and pre-calculating
#' Empirical Cumulative Distribution Functions (ECDF). It produces a
#' structured S3 object optimized for the \code{harvest} function.
#'
#' @param formu A formula object (e.g., \code{Score(5, "high") ~ Group}).
#'   Outcomes on the LHS can be wrapped in function-like calls to specify
#'   clinical thresholds (\code{lambda}) and direction (\code{good}).
#'   Supports multiple outcomes (joined by \code{+}) on the LHS and exactly
#'   one binary grouping variable on the RHS.
#' @param data A data frame containing the variables specified in \code{formu}.
#' @param treated A string specifying the level of the grouping variable to be
#'   treated as the "Treated" (Group X) arm. Defaults to \code{"Treated"}.
#'
#' @return An S3 object of class \code{"sow"}. A nested list containing:
#' \itemize{
#'   \item \code{xid, yid}: Participant IDs for Treated (X) and Control (Y) groups.
#'   \item \code{Variables}: A named list for each outcome containing:
#'     \itemize{
#'       \item \code{x, y}: Raw outcome values for Treated and Control.
#'       \item \code{xs, ys}: Status indicators (1 = observed, 0 = censored). Defaults to all 1 if outcome is not time-to-event.
#'       \item \code{x_weights, y_weights}: Probability weights assigned to each observation (usually 1/n unless it is a time-to-event outcome in which case the weights are based on Kaplan-Meier).
#'       \item \code{x_uni, y_uni}: Unique observed values used for the ECDF.
#'       \item \code{x_w, y_w}: The jump sizes (weights) associated with the unique values.
#'       \item \code{tail_params}: Parameters for the GPD tail extension (for survival data).
#'       \item \code{good}: The direction of superiority ("high" or "low").
#'       \item \code{lambda}: The threshold defining the margin of clinical importance.
#'       \item \code{comment}: Internal metadata regarding default settings applied.
#'       \item \code{name}: The name of the outcome variable.
#'       \item \code{type}: The detected data type (e.g., "binary", "numeric", "survival").
#'     }
#' }
#'
#' @details
#' \bold{Outcome Specification:} The function supports four main types:
#' \code{numeric}, \code{binary}, \code{ordered factor}, and \code{survival}.
#'
#' \bold{In-Formula Arguments:} Thresholds and directions are defined directly
#' in the formula:
#' \itemize{
#'   \item \code{Binary}: Must specify the "superior" level, e.g., \code{Cured("Yes")}.
#'         Sets \code{lambda = 0.5} internally.
#'   \item \code{Ordered}: Must specify direction, e.g., \code{Stage("low")}.
#'         Sets \code{lambda = 1} internally.
#'   \item \code{Numeric}: Requires \code{lambda} and optionally \code{good}
#'         (defaults to \code{"high"}). E.g., \code{Weight(500, "high")}.
#'   \item \code{Survival}: Uses \code{Surv(time, status)(lambda, direction)}.
#' }
#'
#' \bold{Data Partitioning:} The function assigns internal IDs (\code{p_i_d})
#' within each arm. This ensures that when multiple outcomes are analyzed,
#' the joint distribution (correlation) between outcomes for the same individual
#' is preserved.
#'
#' \bold{Tail Extension:} For survival data, the function automatically
#' calculates Generalized Pareto Distribution (GPD) parameters to handle
#' right-censored tails.
#'
#' @examples
#' # ---------------------------------------------------------
#' # Example 1: Continuous Data with Thresholds in Formula
#' # ---------------------------------------------------------
#' data(obstetrics)
#' # Treated patients "win" if birthweight is 500g higher
#' sown_obs <- sow(Birthweight(500, high) ~ Group,
#'                  data = obstetrics,
#'                  treated = "T")
#'
#' # ---------------------------------------------------------
#' # Example 2: Multiple Outcomes (Binary & Continuous)
#' # ---------------------------------------------------------
#' # For Preterm01, "0" is superior. For Birthweight, 500g margin is used.
#' sown_multi <- sow(Preterm01(0) + Birthweight(500) ~ Group,
#'                    data = obstetrics,
#'                    treated = "T")
#'
#' @importFrom survival Surv
#' @importFrom dplyr if_else lag
#' @export
sow <- function(formu, data, treated = "Treated") {

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

  # --- 2. Separates formula components ---
  lhs <- formu[[2]]
  rhs <- formu[[3]]
  get_calls <- function(x) {
    if (is.call(x) && x[[1]] == quote(`+`)) {
      return(c(get_calls(x[[2]]), get_calls(x[[3]])))
    }
    return(list(x))
  }
  calls <- get_calls(lhs)
  n_outcomes <- length(calls)

  # --- 2. Identify Outcomes (LHS) ---
  #lhs_string <- deparse(formu[[2]])
  #outcomes_raw <- stringr::str_trim(unlist(stringr::str_split(lhs_string, "\\+")))
  #n_outcomes <- length(outcomes_raw)

  # Recycle 'good' argument if user provided fewer values than outcomes
  #good_vec <- rep(good, length.out = n_outcomes)
  findgoodlambda <- function(ccc) {
    lccc <- length(ccc) - 1
    if (lccc==0) return(data.frame(lambda=NA,good="high"))
    if (lccc==1 & class(ccc[[2]])=="numeric") return(data.frame(lambda=ccc[[2]],good=NA))
    if (lccc==1 & class(ccc[[2]]) %in% c("name","character") ) return(data.frame(lambda=NA,good=as.character(ccc[[2]])))
    if (lccc==2 & class(ccc[[2]])=="numeric") return(data.frame(lambda=ccc[[2]],good=as.character(ccc[[3]])))
    if (lccc==2 & class(ccc[[2]]) %in% c("name","character")) return(data.frame(lambda=ccc[[3]],good=as.character(ccc[[2]])))
    return(NULL)
  }

  # --- 3. Process Each Outcome ---
  #processed_vars <- lapply(seq_along(outcomes_raw), function(k) {
  processed_vars <- lapply(1:n_outcomes, function(k) {

    if (class(calls[[k]])=="name") {
      calls[[k]] <- as.call(list(calls[[k]]))
    }
    current_call <- calls[[k]]

    is_surv <- "Surv" %in% as.character(current_call[[1]]) & length(as.character(current_call[[1]])) >1

    # Extract raw vectors and status
    if (is_surv) {
      out_type <- "survival"
      time_var <- data[[ as.character(current_call[[1]][[2]]) ]]
      status_var <- data[[ as.character(current_call[[1]][[3]]) ]]
      current_out <- as.character(current_call[[1]][[2]])
    } else {
      current_out <- as.character(current_call[[1]])
      if (length(unique(na.omit(data[[ as.character(current_call[[1]]) ]] ) )) == 2) {
        out_type <- "binary"} else {
          out_type <- class(data[[ as.character(current_call[[1]]) ]])[1] }
      time_var  <- data[[ as.character(current_call[[1]]) ]]
      status_var <- rep(1, length(time_var)) # 1 = observed for continuous
    }


    #######################
    ## Check variable type
    #######################
    ###### Throw an error if the variable is not one of the correct types
    if (!out_type %in% c("survival", "binary", "ordered","numeric","integer")) {
      stop(
        paste("Variable '",
              current_out,
              "' must be one of 'survival', 'binary', 'ordered', 'numeric', or 'integer'. You provided an outcome of type: '",
              out_type,"'",sep=""),
        call. = FALSE
      )
    }

    ##################################
    ### Check if there are no options
    ##################################
    ###### Throw an error if there are no options for binary outcomes
    if (out_type == "binary" & length(current_call) == 1) {
      stop(
        paste("Please specify which value of '",
              current_out,
              "' is superior; '",
              current_out,"(",unique(time_var)[1],")' or '",
              current_out,"(",unique(time_var)[2],")'",sep=""),
        call. = FALSE
      )
    }

    ###### Throw an error if there are no options for ordered outcomes
    if (out_type == "ordered" & length(current_call) == 1) {
      stop(
        paste("Please specify the direction of superior outcomes; '",
              current_out,"(high)' or '",
              current_out,"(low)'",sep=""),
        call. = FALSE
      )
    }

    ###### Throw an error if there are no options for any other outcome
    if (length(current_call) == 1) {
      stop(
        paste("Please specify the lambda, and the direction of superior outcomes; Examples: '",
              current_out,"(10,high)', '",
              current_out,"(20,low)', ",
              current_out,"(100,high)', etc..'",sep=""),
        call. = FALSE
      )
    }

    #######################################
    ### Check if there is only one option
    #######################################
    ###### Throw an error if there is one option for binary outcomes but is not one of the levels
    if (out_type == "binary" & length(current_call) == 2 & !(as.character(current_call[[2]]) %in% unique(time_var)) ) {
      stop(
        paste("The level you have selected is not one of the levels of '",current_out,
              "'\n Please specify which value of '",
              current_out,
              "' is superior; '",
              current_out,"(",unique(time_var)[1],")' or '",
              current_out,"(",unique(time_var)[2],")'",sep=""),
        call. = FALSE
      )
    }

    ###### Throw an error if there is one option for ordered outcomes but is not high or low
    if (out_type == "ordered" & length(current_call) == 2 & !(as.character(current_call[[2]]) %in% c("high","low")) ) {
      stop(
        paste("You need to specify the direction of superior outcomes; '",
              current_out,"(high)' or '",
              current_out,"(low)'",sep=""),
        call. = FALSE
      )
    }

    ###### Throw an error if there is one option for any other outcome and is numeric
    if ( !(out_type %in% c("ordered","binary")) & length(current_call) == 2 & !is.numeric(current_call[[2]]) ) {
      stop(
        paste("You need to specify the lambda for '",
              current_out,"; Example: '",
              current_out,"(20,",current_call[[2]],")' for lambda = 20",sep=""),
        call. = FALSE
      )
    }

    vardef <- findgoodlambda(current_call)
    vardef$comment <- NULL
    if (out_type == "binary" & length(current_call) == 2 & is.numeric(current_call[[2]])) {
      vardef$lambda <- 0.5
      vardef$good <- current_call[[2]]
    }



    if (out_type == "binary") {
      vardef$lambda <- 0.5
      vardef$comment <- paste(" (default); '",current_out," = ",vardef$good,"' is superior.",sep="")
      time_var <- as.character(time_var)
      time_var <- if_else(time_var==vardef$good,1,0)
      vardef$good <- "high"
    }

    if (out_type == "ordered") {
      vardef$lambda <- 1
      vardef$comment <- paste(" (default).",sep="")
    }

    if (!(out_type %in% c("ordered","binary")) & is.na(vardef$good)) {
      vardef$good <- "high"
      vardef$comment <- paste(". Direction = high chosen by default.",sep="")
    }


    # Split into X (Treated) and Y (Control)
    x_val <- time_var[is_treated]
    y_val <- time_var[!is_treated]
    x_stat <- status_var[is_treated]
    y_stat <- status_var[!is_treated]

    x_stat[is.na(x_val)]<-NA
    y_stat[is.na(y_val)]<-NA

    # Handle Character/Factor inputs by converting to numeric
    #if (is.character(x_val) || is.factor(x_val)) {
    #  x_val <- as.numeric(as.factor(x_val)) - 1
    #  y_val <- as.numeric(as.factor(y_val)) - 1
    #}
    #if (out_type == "binary") {
    #  if (is.character(x_val) || is.factor(x_val)) {
    #    x_val <- if_else(x_val==vardef$good,1,0)
    #    y_val <- if_else(y_val==vardef$good,1,0)
    #    vardef$good <- "high"}
    #}


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
      idx <- NULL
      if (p_x$extra_n > 0) {
        ecdf_x <- c_extendtail(ecdf_x$uni, ecdf_x$w, p_x$mlesi, p_x$mlesh, p_x$extra_n)
        idx <- which(x_val > max(x_val[x_stat == 1], na.rm = TRUE))
        if (length(idx) > 0) x_val[idx] <- tail(ecdf_x$uni, length(idx))
      }

      # Process Control
      p_y <- param_gpd(y_val, y_stat)
      idy <- NULL
      if (p_y$extra_n > 0) {
        ecdf_y <- c_extendtail(ecdf_y$uni, ecdf_y$w, p_y$mlesi, p_y$mlesh, p_y$extra_n)
        idy <- which(y_val > max(y_val[y_stat == 1], na.rm = TRUE))
        if (length(idy) > 0) y_val[idy] <- tail(ecdf_y$uni, length(idy))
      }
      tail_params$x <- p_x
      tail_params$y <- p_y
    }

    # --- 6. Internal Weight Calculation ---
    calc_weights <- function(uni, w_vec, raw_vals,is_event) {
      # 1. Standard calculation of cumulative weights
      ord <- order(raw_vals)
      cum_w <- c_ecdf_values(uni, w_vec, raw_vals[ord])

      # 2. Calculate differentials (the "jumps")
      # At this stage, ties usually have the full jump on the first instance
      # and 0 on subsequent instances due to the lag.
      w <- cum_w - dplyr::lag(cum_w, default = 0)

      # 3. Restore original order so we can match weights to raw_vals
      w_original_order <- w[order(ord)]

      # 4. Clean up is_event (Handle NAs and ensure it's logical)
      # This prevents the "NAs not allowed in subscripted assignments" error
      is_event_logical <- as.logical(is_event)
      is_event_logical[is.na(is_event_logical)] <- FALSE

      # 5. Initialize result vector
      final_w <- rep(0, length(raw_vals))

      # 6. Redistribute weights ONLY among valid events
      if(any(is_event_logical)) {
        # We use ave on the subset, then assign back to the logical mask
        final_w[is_event_logical] <- ave(w_original_order[is_event_logical],
                                         raw_vals[is_event_logical],
                                         FUN = function(x) sum(x) / length(x))
      }

      return(final_w)
    }


    x_weights <- calc_weights(ecdf_x$uni, ecdf_x$w, x_val,x_stat)
    y_weights <- calc_weights(ecdf_y$uni, ecdf_y$w, y_val,y_stat)

    ##########################################################
    ### Recalculates the weights in the extended tail
    ### at the moment they have 0 weights so I just re-weight
    ### to make sure that total weight is 1
    ##########################################################
    if (out_type == "Survival") {
      if (!is.null(idx)) {
        x_weights[idx] <- (1-sum(x_weights,na.rm=TRUE)) / sum(!is.na(x_val[idx]))
      }

      if (!is.null(idy)) {
        y_weights[idy] <- (1-sum(y_weights,na.rm=TRUE)) / sum(!is.na(y_val[idy]))
      }
    }

    list(
      x = x_val, xs = x_stat, x_weights = x_weights,
      y = y_val, ys = y_stat, y_weights = y_weights,
      x_uni = ecdf_x$uni, x_w = ecdf_x$w,
      y_uni = ecdf_y$uni, y_w = ecdf_y$w,
      tail_params = tail_params,
      #good = good_vec[k],
      good = vardef$good,
      lambda = vardef$lambda,
      comment = vardef$comment,
      name = current_out,
      type = out_type
    )
  })

  # --- 7. Finalize Names and Object ---
  # Cleans Surv(Time, Status) into just "Time" for list names
  #clean_names <- stringr::str_remove_all(outcomes_raw, "Surv\\(|\\)|\\s") |>
  #  stringr::str_split(",") |>
  #  sapply(`[[`, 1)
  clean_names <- sapply(processed_vars,function(x) x$name)

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
