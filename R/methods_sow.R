#' Plot Method for sow Objects
#' @param x An object of class \code{sow}.
#' @param ... Additional arguments passed to ggplot2.
#' @import ggplot2
#' @export
plot.sow <- function(x, ...) {

  # 1. Prepare data for ggplot
  plot_list <- lapply(names(x$Variables), function(v_name) {
    v <- x$Variables[[v_name]]

    # Combine X and Y data for this variable
    d <- rbind(
      data.frame(uni = v$x_uni,w = v$x_w, group = "Treated (X)", outcome = v_name),
      data.frame(uni = v$y_uni,w = v$y_w, group = "Control (Y)", outcome = v_name)
    )

    # Add anchor points at (min_x, 0) to ensure the line drops to the bottom
    anchors <- data.frame(
      uni = rep(min(d$uni, na.rm = TRUE), 2),
      w = c(0, 0),
      group = c("Treated (X)", "Control (Y)"),
      outcome = v_name
    )

    return(rbind(anchors, d))
  })

  df_plot <- do.call(rbind, plot_list)

  # 2. Generate the plot
  p <- ggplot(df_plot, aes(x = uni,y=w, color = group)) +
    geom_step(linewidth = 1) +
    facet_wrap(~outcome, scales = "free_x") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
      title = "Empirical Cumulative Distribution Functions",
      subtitle = "Visual check of outcome distributions by group",
      x = "Observed Value",
      y = "Cumulative Probability",
      color = "Study Arm"
    )

  return(p)
}


#' Print Method for sow Objects
#' @param x An object of class \code{sow}.
#' @param ... Additional arguments (ignored).
#' @export
print.sow <- function(x, ...) {
  n_vars <- length(x$Variables)
  cat("\n--- infuse: Sown Data Object ---\n")
  cat(sprintf("Outcomes processed: %d (%s)\n",
              n_vars, paste(names(x$Variables), collapse = ", ")))
  cat(sprintf("Total Sample Size: %d Treated (X), %d Control (Y)\n Overall N = %d\n",
              length(x$xid), length(x$yid) ,length(x$xid) + length(x$yid) ))
  cat("------------------------------------------------------------------------------------------------\n")
  #cat("Next step: harvest(fit)\n")

  temp <- lapply(x$Variables,function(xx){
    #data.frame(Variable=xx$name,Type = xx$type,Direction = xx$good,Lambda = xx$lambda)
    cat(sprintf("\033[1m%s\033[22m (%s); \033[4m%s\033[24m values superior; Lambda = %g",
                xx$name, xx$type,xx$good,xx$lambda))
    cat(xx$comment,"\n")
  })
  cat("------------------------------------------------------------------------------------------------\n")

  invisible(x)
}



#' Summary Method for sow Objects
#' @param object An object of class \code{sow}.
#' @param ... Additional arguments (ignored).
#' @export
summary.sow <- function(object, ...) {
  cat("\nStatistical Summary of Sown Data\n")
  cat(rep("=", 70), "\n", sep = "")

  # Create a summary table for each variable
  summary_list <- lapply(names(object$Variables), function(v_name) {
    v <- object$Variables[[v_name]]
    data.frame(
      Outcome   = v_name,
      Type      = v$type,
      Direction = v$good,
      Mean_X    = round(mean(v$x, na.rm = TRUE), 3),
      Mean_Y    = round(mean(v$y, na.rm = TRUE), 3),
      N_Event_X = sum(v$xs, na.rm = TRUE),
      N_Event_Y = sum(v$ys, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  summary_table <- do.call(rbind, summary_list)
  print(summary_table, row.names = FALSE)

  # Check if GPD was used
  gpd_check <- any(sapply(object$Variables, function(v) v$tail_params$x$extra_n > 0))
  if (gpd_check) {
    cat("\nNote: Survival tail extension (GPD) was applied to one or more variables.\n")
  }

  cat(rep("=", 70), "\n", sep = "")
}
