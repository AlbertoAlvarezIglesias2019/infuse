#' Prepare Data for Win Statistics Calculations
#'
#' @description
#' A highly optimized preprocessing utility that "infuses" the analysis by parsing outcomes,
#' partitioning data into treatment/control groups, and pre-calculating
#' Empirical Cumulative Distribution Functions (ECDF). It utilizes \code{data.table}
#' for fast, in-place memory modifications and produces a structured S3 object
#' optimized for the \code{realise} function.
#'
#' @param formu A formula object (e.g., \code{Score(5, "high") ~ Group}).
#'   Outcomes on the LHS can be wrapped in function-like calls to specify
#'   clinical thresholds (\code{lambda}) and direction (\code{direction}).
#'   Supports multiple outcomes (joined by \code{+}) on the LHS and exactly
#'   one binary grouping variable on the RHS.
#' @param data A data frame or data.table containing the variables specified in \code{formu}.
#' @param treated A string specifying the level of the grouping variable to be
#'   treated as the "Treated" (Group X) arm. Defaults to \code{"Treated"}.
#'
#' @return An S3 object of class \code{"infuse"}. A list containing three stacked \code{data.table}s
#' and a \code{groups} attribute:
#' \itemize{
#'   \item \code{data}: Patient-level data combined across all outcomes. Contains:
#'     \itemize{
#'       \item \code{p_i_d}: Unique participant ID within their arm.
#'       \item \code{arm}: Group assignment (factor: "X" for treated, "Y" for control).
#'       \item \code{name}: The name of the outcome variable (factor).
#'       \item \code{time}: The observed value or time-to-event.
#'       \item \code{status}: Status indicator (1 = observed, 0 = censored/missing).
#'       \item \code{prob}: ECDF probability.
#'       \item \code{weight}: Probability weight (jump size) for the observation.
#'       \item \code{ord}: Original ordering index.
#'       \item \code{p_i_dtail}: Indicator for right-censored tail events (1 = tail, 0 = non-tail).
#'       \item \code{f, u}: Individual patient-level influence function values for favorable and unfavorable probabilities.
#'     }
#'   \item \code{metadata}: Variable-level configurations and overall probabilities. Contains:
#'     \itemize{
#'       \item \code{name}: The name of the outcome variable.
#'       \item \code{type}: Detected data type (e.g., "binary", "numeric", "survival").
#'       \item \code{lambda}: Threshold defining the margin of clinical importance.
#'       \item \code{direction}: Direction of superiority ("high" or "low").
#'       \item \code{comment}: Internal metadata regarding default settings applied.
#'       \item \code{f, u}: Overall favorable (A) and unfavorable (B) win probabilities for the outcome.
#'       \item \code{ssx, ssy}: Effective sample size (observed events) for the Treated (X) and Control (Y) arms.
#'     }
#'   \item \code{tail}: Generalized Pareto Distribution (GPD) tail extension parameters.
#' }
#'
#' @details
#' \bold{Performance Optimization:} This function utilizes \code{data.table} to perform
#' in-place calculations (\code{:=}) and fast row-binding (\code{rbindlist}). This eliminates
#' unnecessary memory copies, making the ECDF mapping and tail extension calculations extremely fast.
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
#'                    data = obstetrics,
#'                    treated = "T")
#'
#' # ---------------------------------------------------------
#' # Example 2: Multiple Outcomes (Binary & Continuous)
#' # ---------------------------------------------------------
#' # For PretermYN, "Yes" is superior. For Birthweight, 500g margin is used.
#' infused_multi <- infuse(PretermYN("Yes") + Birthweight(500) ~ Group,
#'                      data = obstetrics,
#'                      treated = "T")
#'
#' @import data.table
#' @importFrom survival Surv
#' @export
infuse <- function(formu, data, treated = "Treated") {

  library(data.table, quietly = TRUE)

  #data <- as.data.frame(data)
  Data <- as.data.table(data)


  # --- 1. Arm Identification & Partitioning ---
  arm_var <- all.vars(formu[[3]])[1]
  if (length(unique(Data[[arm_var]])) != 2) {
    stop("Grouping variable must have exactly two levels.")
  }

  # creates a counter 1,2,3... for each unique value of the column named in arm_var
  Data[, p_i_d := rowid(get(arm_var))]


  #################################
  ## Reconfigure Grouping variable
  #################################
  setnames(Data, arm_var, "arm")
  groups_cmt <- paste0("Grouping variable: ",arm_var,";\nGroups: X = ",treated,"; Y = ",Data[arm != treated, unique(arm)] )
  Data[, arm := fcase(arm == treated, "X",arm != treated, "Y")]

  # --- 2. Separates formula components ---
  lhs <- formu[[2]]
  get_calls <- function(x) {
    if (is.call(x) && x[[1]] == quote(`+`)) {
      return(c(get_calls(x[[2]]), get_calls(x[[3]])))
    }
    return(list(x))
  }
  lhs_calls <- get_calls(lhs)
  n_outcomes <- length(lhs_calls)


  # --- 3. Process Each Outcome ---
  processed_vars <- lapply(1:n_outcomes, function(k) {

    D <- copy(Data)

    if (class(lhs_calls[[k]])=="name") {
      lhs_calls[[k]] <- as.call(list(lhs_calls[[k]]))
    }
    current_call <- lhs_calls[[k]]
    current_out <- as.character(current_call[[1]])

    is_surv <- "Surv" %in% current_out & length(current_out) >1

    # Extract vectors and status names
    if (is_surv) {
      out_type <- "survival"
      current_out <- as.character(current_call[[1]][[2]])
      time_var_name <- current_out
      status_var_name <- as.character(current_call[[1]][[3]])
      setnames(D, c(time_var_name,status_var_name), c("time","status") )

      #D$time <- data[[ time_var_name ]]
      #D$status <- data[[ status_var_name ]]
    } else {
      if (length(unique(na.omit(D[[ current_out ]] ) )) == 2) {
        out_type <- "binary"} else {out_type <- class(Data[[ current_out ]])[1] }
      time_var_name <- current_out
      status_var_name <- NULL
      setnames(D, time_var_name, "time" )
      D[, status := NA_real_]
      D[!is.na(time), status := 1]
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
        stop(paste0("Please specify which value of '", current_out, "' is superior; '", current_out, "(", unique(time_var)[1], ")' or '", current_out, "(", unique(time_var)[2], ")'"), call. = FALSE)
      } else if (out_type == "ordered") {
        stop(paste0("Please specify the direction of superior outcomes; '", current_out, "(high)' or '", current_out, "(low)'"), call. = FALSE)
      } else {
        stop(paste0("Please specify the lambda, and the direction of superior outcomes; Examples: '", current_out, "(10,high)', '", current_out, "(20,low)', '", current_out, "(100,high)', etc..'"), call. = FALSE)
      }
    }

    # 3. Handle invalid options (length == 2, Only one option)
    if (n_args == 2) {
      opt_val <- current_call[[2]]
      if (out_type == "binary" && !as.character(opt_val) %in% unique(D$time)) {
        stop(
          paste0("The level you have selected is not one of the levels of '",current_out,"'\n Please specify which value of '",
                 current_out,"' is superior; '",current_out,"(",unique(time_var)[1],")' or '",current_out,
                 "(",unique(time_var)[2],")'"),call. = FALSE)
      } else if (out_type == "ordered" && !(as.character(current_call[[2]]) %in% c("high","low"))) {
        stop(
          paste0("You need to specify the direction of superior outcomes; '",current_out,"(high)' or '",
                 current_out,"(low)'"),call. = FALSE)
      } else if (!(out_type %in% c("ordered","binary")) && !is.numeric( opt_val )) {
        stop(
          paste0("You need to specify the lambda for '",current_out,"; Example: '",current_out,"(20,",
                 current_call[[2]],")' for lambda = 20"),call. = FALSE)
      }
    }

    # 4. Read the lambda and the good from the formula
    lam <- NA_real_
    gd <- "high"
    if (n_args==2) {
      arg2 <- current_call[[2]]
      if (is.numeric(arg2)) {lam <- arg2;gd <- NA_character_} else {gd <- as.character(arg2)}
    } else if (n_args == 3) {
      arg2 <- current_call[[2]]
      arg3 <- current_call[[3]]
      if (is.numeric(arg2)) {lam <- arg2;gd <- as.character(arg3)} else {lam <- arg3;gd <- as.character(arg2)}
    }

    # Initialize comment as NA (better than NULL for data frame columns)
    cmt <- NA_character_
    if(out_type == "binary") {
      # Use && for speed and n_args from the previous cache
      if (n_args == 2 && is.numeric(current_call[[2]])) {
        gd <- current_call[[2]]
      }
      lam <- 0.5
      cmt <- paste0("Lambda = 0.5 by default; '", current_out, " = ", gd, "' is superior")
      D$time <- as.integer(as.character(D$time) == gd)
      gd <- "high"
    } else if (out_type == "ordered") {
      lam <- 1
      cmt <- "Lambda = 1 by default"
    } else {
      # For survival, numeric, or integer types
      if (is.na(gd)) {
        gd <- "high"
        cmt <- "Direction = high chosen by default"
      }
    }
    #vardef <- data.frame(lambda = lam,good = gd,comment = cmt,stringsAsFactors = FALSE)
    #vardef <- data.table(lambda = lam, good = gd, comment = cmt)

    D[, .(p_i_d, arm, name, time, status)]
    D_meta <- data.table(name = current_out,type = out_type,lambda = lam,direction = gd,comment = cmt)

    #D_meta <- data.frame(name=current_out,type=out_type)
    #D_meta <- cbind(D_meta,vardef)

    D[, c("prob", "weight", "ord") := {
      res <- c_ecdf_plus(time, status)
      list(res$p, res$w, res$o)
    }, by = arm]


    # 1. Pre-allocate the final params table with defaults (0)
    # Doing this as a data.table allows us to update it instantly inside the loop
    tail_params <- data.table(arm = c("X", "Y"),mlesi = 0,mlesh = 0,extra_n = 0)

    D[,p_i_dtail:=0]

    if (is_surv) {

      # 1. Calculate tail indicator in place (grouped by arm, zero memory copies)
      D[, p_i_dtail := as.integer(time > max(time * status, na.rm = TRUE)), by = arm]

      # Replicate `ddd[ddd$ord, ]` - we sort the whole table by arm, then by ord
      setorder(D, arm, ord)

      # 2. Loop over arms to prevent duplicated code
      for (a in c("X", "Y")) {

        # Count tail events for this specific arm
        ext <- D[arm == a, sum(p_i_dtail == 1, na.rm = TRUE)]

        if (ext > 0) {
          # Extract current arm data once for speed
          arm_data <- D[arm == a]

          # Calculate GPD parameters
          ppp <- param_gpd(arm_data$time, arm_data$status)

          # 2. Update tail_params IN PLACE for the current arm
          # No if-statements needed. It just finds the row where arm == a and updates it.
          tail_params[arm == a, `:=`(
            mlesi   = ppp$mlesi,
            mlesh   = ppp$mlesh,
            extra_n = ext
          )]


          # Run tail calculations using the non-tail subset
          non_tail <- arm_data[p_i_dtail == 0]
          ta <- c_extendtail(non_tail$time, non_tail$prob, ppp$mlesi, ppp$mlesh, ext)

          # Calculate total weight for this arm
          sum_w <- sum(arm_data$weight, na.rm = TRUE)

          # 3. Update the tail rows IN PLACE in the main dataset
          D[arm == a & p_i_dtail == 1, `:=`(
            time   = tail(ta$uni, n = ext),
            prob   = tail(ta$w, n = ext),
            weight = (1 - sum_w) / ext
          )]
        }
      }

      # 4. Final sort by p_i_d (replaces arrange(p_i_d))
      setorder(D, p_i_d)

    }
    tail_params[, name := current_out]



    ##############################
    ### Favorable and Unfavorable
    ##############################
    # 1. Fast subsetting and extraction (assuming todoslosdatos is a data.table)
    dx <- D[arm == "X"]
    dy <- D[arm == "Y"]

    # Extract X vectors
    x <- as.numeric(dx$time)
    ssx <- sum(dx$status,na.rm=TRUE)
    xw <- dx$weight
    or_x <- dx$ord
    x_sorted <- x[or_x]
    px_sorted <- dx$prob[or_x]

    # Extract Y vectors
    y <- as.numeric(dy$time)
    ssy <- sum(dy$status,na.rm=TRUE)
    yw <- dy$weight
    or_y <- dy$ord
    y_sorted <- y[or_y]
    py_sorted <- dy$prob[or_y]

    lambda <- D_meta$lambda

    # 2. FAVORABLE
    # Evaluate the ECDF once and reuse it for both A and ify_A
    pred_y_plus <- 1 - c_ecdf_predict_less(y + lambda, x_sorted, px_sorted)

    A <- sum(yw * pred_y_plus, na.rm = TRUE)
    ifx_A <- c_ecdf_predict(x - lambda, y_sorted, py_sorted) - A
    ify_A <- pred_y_plus - A # Reusing the stored calculation


    # 3. UNFAVORABLE
    # Evaluate the ECDF once and reuse it for both B and ifx_B
    pred_x_plus <- 1 - c_ecdf_predict_less(x + lambda, y_sorted, py_sorted)

    B <- sum(xw * pred_x_plus, na.rm = TRUE)
    ifx_B <- pred_x_plus - B # Reusing the stored calculation
    ify_B <- c_ecdf_predict(y - lambda, x_sorted, px_sorted) - B


    # 4. Swap logic (cleaned up)
    if (D_meta$direction != "high") {
      # Simple temporary variable swap
      tmp <- A; A <- B; B <- tmp
      tmp <- ifx_A; ifx_A <- ifx_B; ifx_B <- tmp
      tmp <- ify_A; ify_A <- ify_B; ify_B <- tmp
    }


    # 5. Assign 'f' and 'u' back to the main dataset IN PLACE
    D[arm == "X", `:=`(f = ifx_A, u = ifx_B)]
    D[arm == "Y", `:=`(f = ify_A, u = ify_B)]

    D_meta[, f := A]
    D_meta[, u := B]

    D_meta[, ssx := ssx]
    D_meta[, ssy := ssy]

    #D_meta$f <- A
    #D_meta$u <- B
    if (out_type == "ordered") {
      D$timelab <- as.character(D$time)
      D$time <- as.numeric(D$time)
    }


    #D_meta$ssx <-ssx
    #D_meta$ssy <-ssy

    list(D = D,D_meta = D_meta,D_tail=tail_params)

  })

  sow_data     <- rbindlist(lapply(processed_vars, "[[", "D"), fill = TRUE)
  sow_metadata <- rbindlist(lapply(processed_vars, "[[", "D_meta"), fill = TRUE)
  sow_tail     <- rbindlist(lapply(processed_vars, "[[", "D_tail"), fill = TRUE)

  # --- NEW: Enforce ordering based on metadata ---
  ordered_names <- unique(sow_metadata$name)
  sow_data[, name := factor(name, levels = ordered_names)]
  sow_data[, arm := factor(arm, levels = c("X","Y"))]

  structure(
    list(
      data=sow_data,
      metadata = sow_metadata,
      tail = sow_tail
      #data=setDF(sow_data),
      #metadata = setDF(sow_metadata),
      #tail = setDF(sow_tail)
    ),
    class = "infuse",
    groups = groups_cmt
  )
}
