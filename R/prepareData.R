#' Prepare Data for Endpoint Analysis
#'
#' Processes a given formula and dataset to standardize endpoints (survival, binary,
#' ordered, numeric, integer), define treatment arms, extract metadata, and apply
#' GPD extensions for survival data.
#'
#' @param formu A \code{\link[stats]{formula}} specifying the outcomes on the
#'   left-hand side and the grouping variable on the right-hand side. The LHS
#'   can contain multiple outcomes separated by \code{+} and wrapped in parameters,
#'   e.g., \code{bp(2) + Surv(time, status)(20, high) ~ treatment}.
#' @param data A data.frame or \code{data.table} containing the dataset.
#' @param treated A character string specifying the level in the grouping
#'   variable that represents the treated arm. Defaults to \code{"Treated"}.
#'
#' @return A named list containing three \code{data.table}s:
#' \describe{
#'   \item{data}{The processed dataset in long format. Columns include \code{name},
#'   \code{type}, \code{lambda}, \code{direction}, \code{p_i_d} (patient ID), \code{arm},
#'   \code{time}, \code{status}, and GPD extension variables \code{time_ext} and \code{status_ext}.}
#'   \item{metadata}{Extracted metadata mapping each outcome \code{name} to its corresponding \code{comment}.}
#'   \item{tail}{Tail parameters derived from the GPD extension for survival outcomes.}
#' }
#'
#' @import data.table
#' @export
prepareData <- function(formu, data, treated = "Treated") {

  # Coerce to data.table and copy to avoid modifying the user's global object by reference
  Data <- copy(as.data.table(data))


  # --- 1. Arm Identification & Partitioning ---
  arm_var <- all.vars(formu[[3]])[1]

  # uniqueN is data.table's highly optimized C-level unique counter
  if (uniqueN(Data[[arm_var]]) != 2) {
    stop("Grouping variable must have exactly two levels.")
  }

  # Create a fast sequential ID and the treatment arm categorization
  Data[, p_i_d := .I]
  Data[, arm := fifelse(get(arm_var) == treated, "X", "Y")]

  # --- 2. Separate formula components ---
  lhs <- formu[[2]]
  get_calls <- function(x) {
    if (is.call(x) && x[[1]] == quote(`+`)) {
      return(c(get_calls(x[[2]]), get_calls(x[[3]])))
    }
    return(list(x))
  }

  lhs_calls <- get_calls(lhs)
  n_outcomes <- length(lhs_calls)


  # Process each outcome
  processed_vars <- lapply(1:n_outcomes, function(k) {

    current_call <- lhs_calls[[k]]
    if (inherits(current_call, "name")) {
      current_call <- as.call(list(current_call))
    }

    # Extract components of the call
    call_inner <- current_call[[1]]
    call_char  <- as.character(call_inner)

    is_surv <- "Surv" %in% call_char && length(call_char) > 1

    # Identify variables without copying the whole dataset
    if (is_surv) {
      out_type <- "survival"
      time_var_name   <- call_char[2]
      status_var_name <- call_char[3]
      current_out     <- time_var_name

      # Extract only required columns
      D <- Data[, .(p_i_d, arm, time = get(time_var_name), status = get(status_var_name))]

    } else {
      current_out     <- call_char[1]
      time_var_name   <- current_out
      status_var_name <- NULL

      # Extract only required columns, initialize status safely
      D <- Data[, .(p_i_d, arm, time = get(time_var_name), status = NA_real_)]
      D[!is.na(time), status := 1]

      # Fast check for binary using uniqueN
      if (uniqueN(D[!is.na(time), time]) == 2) {
        out_type <- "binary"
      } else {
        out_type <- class(D$time)[1]
      }
    }

    D[, name := time_var_name]


    # 1. Validate variable type
    valid_types <- c("survival", "binary", "ordered", "numeric", "integer")
    if (!out_type %in% valid_types) {
      stop(paste0("Variable '", current_out, "' must be one of 'survival', 'binary', 'ordered', 'numeric', or 'integer'. You provided an outcome of type: '", out_type, "'"), call. = FALSE)
    }

    n_args <- length(current_call)

    # 2. Handle missing options (length == 1)
    if (n_args == 1) {
      if (out_type == "binary") {
        time_levels <- D[!is.na(time), unique(time)]
        stop(sprintf("Please specify which value of '%s' is superior; '%s(%s)' or '%s(%s)'",
                     current_out, current_out, time_levels[1], current_out, time_levels[2]), call. = FALSE)
      } else if (out_type == "ordered") {
        stop(sprintf("Please specify the direction of superior outcomes; '%s(high)' or '%s(low)'",
                     current_out, current_out), call. = FALSE)
      } else {
        stop(sprintf("Please specify the lambda and direction of superior outcomes. Example: '%s(10,high)', '%s(20,low)'",
                     current_out, current_out), call. = FALSE)
      }
    }

    # 3. Handle invalid single options (length == 2, Only one option)
    if (n_args == 2) {
      opt_val <- current_call[[2]]
      if (out_type == "binary" && !as.character(opt_val) %in% as.character(D$time)) {
        time_levels <- D[!is.na(time), unique(time)]
        stop(sprintf("The level you selected is not one of the levels of '%s'.\nPlease specify which value is superior; '%s(%s)' or '%s(%s)'",
                     current_out, current_out, time_levels[1], current_out, time_levels[2]), call. = FALSE)
      } else if (out_type == "ordered" && !as.character(opt_val) %in% c("high", "low")) {
        stop(sprintf("You need to specify the direction of superior outcomes; '%s(high)' or '%s(low)'",
                     current_out, current_out), call. = FALSE)
      } else if (!out_type %in% c("ordered", "binary") && !is.numeric(opt_val)) {
        stop(sprintf("You need to specify the lambda for '%s'; Example: '%s(20,%s)' for lambda = 20",
                     current_out, current_out, as.character(opt_val)), call. = FALSE)
      }
    }

    # 4. Extract parameters (Lambda and Direction)
    lam <- NA_real_
    gd <- "high"
    if (n_args == 2) {
      arg2 <- current_call[[2]]
      if (is.numeric(arg2)) {
        lam <- arg2; gd <- NA_character_
      } else {
        gd <- as.character(arg2)
      }
    } else if (n_args == 3) {
      arg2 <- current_call[[2]]
      arg3 <- current_call[[3]]
      if (is.numeric(arg2)) {
        lam <- arg2; gd <- as.character(arg3)
      } else {
        lam <- as.numeric(arg3); gd <- as.character(arg2)
      }
    }

    # Process types and build meta comments
    cmt <- NA_character_

    if (out_type == "binary") {
      if (n_args == 2 && is.numeric(current_call[[2]])) {
        gd <- as.character(current_call[[2]])
      }
      lam <- 0.5
      cmt <- sprintf("Lambda = 0.5 by default; '%s = %s' is superior", current_out, gd)
      D[, time := as.integer(as.character(time) == gd)]
      gd <- "high"

    } else if (out_type == "ordered") {
      lam <- 1
      cmt <- "Lambda = 1 by default"
      D[, timelab := as.character(time)]
      D[, time := as.numeric(time)]

    } else {
      if (is.na(gd)) {
        gd <- "high"
        cmt <- "Direction = high chosen by default"
      }
    }

    # Initialize tail params (fast row creation)
    tail_params <- data.table(arm = c("X", "Y"), mlesi = 0, mlesh = 0, extra_n = 0, name = current_out)

    D[,`:=`(status_ext = 0 , time_ext = 0)]
    if (is_surv) {
      # Optimize loop by using data.table integer indexing (.I) to safely assign back
      for (a in c("X", "Y")) {
        idx <- D[, .I[arm == a]]

        if (length(idx) > 0) {
          ppp <- gpdextension(D$time[idx], D$status[idx])

          # Update Tail Parameters table
          tail_params[arm == a, `:=`(mlesi = ppp$mlesi, mlesh = ppp$mlesh, extra_n = ppp$extra_n)]

          # Update Data directly via indices
          D[idx, `:=`(time_ext = ppp$x, status_ext = ppp$xtail)]
        }
      }
      D[status_ext == 0, time_ext := 0]
    }

  D[,`:=`(type = out_type , lambda = lam,direction = gd)]

  D_meta <- data.table(name = current_out, comment = cmt)

  list(D = D, D_meta = D_meta, D_tail = tail_params)
  })

  # Consolidate outputs via optimized rbindlist
  sow_data     <- rbindlist(lapply(processed_vars, "[[", "D"), fill = TRUE)
  sow_metadata <- rbindlist(lapply(processed_vars, "[[", "D_meta"), fill = TRUE)
  sow_tail     <- rbindlist(lapply(processed_vars, "[[", "D_tail"), fill = TRUE)

  # Enforce factor levels cleanly
  ordered_names <- unique(sow_metadata$name)
  sow_data[, `:=`(
    name = factor(name, levels = ordered_names),
    arm  = factor(arm, levels = c("X", "Y"))
  )]

  setcolorder(sow_data, c("name", "type","lambda", "direction","p_i_d","arm","time","status","time_ext","status_ext"))

  list(data = sow_data, metadata = sow_metadata, tail = sow_tail)


}
