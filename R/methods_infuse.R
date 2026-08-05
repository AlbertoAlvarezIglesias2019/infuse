#' Plot Method for infuse Objects
#'
#' @description
#' Visualizes the Empirical Cumulative Distribution Functions (ECDF) for all outcomes
#' processed by \code{infuse}. The plot automatically includes clinical threshold
#' (\code{lambda}) annotations and highlights simulated survival tail extensions.
#'
#' @param x An object of class \code{infuse}.
#' @param ... Additional arguments passed to \code{ggplot2}.
#'
#' @details
#' The plot displays the ECDF for the Treated (X) and Control (Y) arms across all
#' faceted outcomes. If a survival outcome underwent GPD tail extension, those
#' simulated points are connected via a black step line to distinguish them from
#' observed data.
#'
#' @import ggplot2
#' @import data.table
#' @export
plot.infuse <- function(x, ...) {

  # Ensure ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is required to plot infuse objects. Please install it.")
  }
  library(ggplot2)
  library(data.table)

  dt <- copy(x$data)
  dt$ttt <- dt$time
  dt$ttt[dt$status_ext==1] <- dt$time_ext[dt$status_ext==1]
  dt$sss <- dt$status
  dt$sss[dt$status_ext==1] <- 1

  dt[, c("prob", "ord") := {
    res <- c_ecdf_plus(ttt, sss)
    list(res$p, res$o)
  }, by = .(name, arm)]


  # --- NEW: Anchor points for 0 and 1 tails ---
  # We expand the time range by 5% on each side to show the horizontal "tails"
  #rng <- diff(range(dt$ttt, na.rm = TRUE))
  #rng <- if (length(rng) == 0 || is.na(rng) || rng == 0) 1 else rng

  #dt_bounds <- dt[, .(
  #  t_min = min(ttt, na.rm = TRUE) - rng * 0.05,
  #  t_max = max(ttt, na.rm = TRUE) + rng * 0.05
  #), by = .(name, arm)]

  dt_bounds <- dt[, {
    r <- diff(range(ttt, na.rm = TRUE))
    rng <- if (length(r) == 0 || is.na(r) || r <= 0) 1 else r

    .SD[, .(
      t_min = min(ttt, na.rm = TRUE) - rng * 0.05,
      t_max = max(ttt, na.rm = TRUE) + rng * 0.05
    ), by = arm]
  }, by = name]


  anchors <- rbind(
    dt_bounds[, .(name, arm, ttt = t_min, prob = 0, status_ext = 0)],
    dt_bounds[, .(name, arm, ttt = t_max, prob = 1, status_ext = 0)]
  )

  dt <- rbindlist(list(dt, anchors), fill = TRUE)

  # Ensure data is sorted by arm and time for proper ECDF step plotting
  setorder(dt, name, arm, ttt,prob)

  # --- NEW: Build Annotation Data for Lambda ---
  # 1. Find the horizontal center of each facet
  lambda_dt <- dt[, .(
    x_min = min(ttt, na.rm = TRUE),
    x_max = max(ttt, na.rm = TRUE),
    lambda = mean(lambda,na.rm=TRUE)
  ), by = name]

  # 2. Merge the lambda values from metadata
  meta_dt <- x$metadata
  lambda_dt <- merge(lambda_dt, meta_dt[, .(name)], by = "name")

  # 3. Calculate segment coordinates (centered) and put them at the top (y = 1.04)
  lambda_dt[, x_mid := (x_min + x_max) / 2]
  lambda_dt[, x_start := x_mid - (lambda / 2)]
  lambda_dt[, x_end := x_mid + (lambda / 2)]
  lambda_dt[, y_pos := 1.04]
  lambda_dt[, label_text := paste0("\u03bb = ", lambda)] # \u03bb is the Greek letter λ


  # 1. Base faceted plot
  p <- ggplot(dt, aes(x = ttt, y = prob, color = arm, group = arm))+

    # --- NEW: Add the horizontal reference line at y = 1 ---
    geom_hline(yintercept = 1, linetype = "dashed", color = "#B0B0B0", linewidth = 0.6) +
    geom_step(linewidth = 1,na.rm = TRUE) +

    # --- NEW: Draw Lambda Segment & Text (Only if lambda > 0) ---
    geom_segment(
      data = lambda_dt[lambda > 0],
      aes(x = x_start, xend = x_end, y = y_pos, yend = y_pos),
      inherit.aes = FALSE, # Prevent it from looking for 'arm' or 'prob'
      color = "#555555",
      linewidth = 0.6,
      arrow = arrow(ends = "both", length = unit(0.15, "cm"))
    ) +
    geom_text(
      data = lambda_dt[lambda > 0],
      aes(x = x_mid, y = y_pos + 0.035, label = label_text),
      inherit.aes = FALSE,
      size = 3,
      color = "#333333",
      fontface = "italic"
    ) +

    # --- NEW: Expand Y-axis to make room for annotations ---
    scale_y_continuous(
      limits = c(0, 1.1),          # Forces extra space at the top
      breaks = seq(0, 1, by = 0.25) # Keeps the y-axis ticks clean
    )+
    scale_color_manual(values = c("X" = "#00BFC4", "Y" = "#F8766D"), name = "Arm",
                       labels = c("X" = "X (treated)", "Y" = "Y (controls)")) +
    facet_wrap(~ name, scales = "free_x") +
    labs(
      title = "Win/Loss Probabilities by Outcome",
      x = "Value / Time",
      y = "Probability"
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
          # Makes main title 50% larger
          plot.title = element_text(size = rel(1.2))) +
    labs(
      title = "Empirical Cumulative Distribution Functions",
      subtitle = "Visual check of outcome distributions by group",
      x = "Observed Value",
      y = "Cumulative Probability",
      color = "Study Arm"
    )



  # 2. Highlight the extended tail in black (across all facets)
  if (any(dt$status_ext == 1)) {

    # Safely extract tail points + the final non-tail point, grouped by BOTH name and arm
    tail_data <- dt[, {
      idx_non_tail <- which(status_ext == 0)
      idx_tail <- which(status_ext == 1)

      # Get the index of the last observed point before the tail begins
      last_obs_idx <- if(length(idx_non_tail) > 0) idx_non_tail[which.max(ttt[idx_non_tail])] else integer(0)

      .SD[c(last_obs_idx, idx_tail)]
    }, by = .(name, arm)]

    # Overlay the black tail connecting line
    p <- p + geom_step(data = tail_data, color = "black", linewidth = 1)

    # Add points to make the simulated extended tail obvious
    #p <- p + geom_point(data = dt[p_i_dtail == 1], color = "black", size = 1, alpha = 0.6)
  }

  return(p)
}


#' Print Method for infuse Objects
#'
#' @param x An object of class \code{infuse}.
#' @param ... Additional arguments (ignored).
#' @export
print.infuse <- function(x, ...) {

  dt <- unique(x$data[, .(name, type, lambda, direction)])[x$metadata, on = "name"]

  #dt <- copy(x$metadata)

  cat("Object of class 'infuse' (Prepared Win/Loss Data)\n")
  cat(rep("=", 46), "\n", sep = "")
  cat("\n--- Overall summary ---\n")
  # 1. High-level summary
  n_outcomes <- nrow(dt)
  arms <- sort(unique(x$data$arm))

  cat(sprintf("Outcomes processed: %d (%s)\n",
              n_outcomes, paste(dt$name, collapse = ", ")))
  cat(sprintf(attr(x, "groups")))
  tt <- x$data[, .N, by = .(name, arm)]
  nx <- tt$N[tt$arm=="X"][1]
  ny <- tt$N[tt$arm=="Y"][1]
  cat(sprintf("\nSample Size: %d Treated (X), %d Control (Y)\nOverall N = %d\n",
              nx, ny ,nx+ny ))

  # 2. Outcome Specifications & Probabilities
  cat("\n--- Outcome Summary ---\n")

  setnames(dt, old = c("name","type","direction","lambda","comment"), new = c("Variable", "Type","Direction","Lambda","Comment"))
  print(dt[, .(Variable, Type, Direction, Lambda, Comment)], row.names = FALSE,class = FALSE)


  # 3. Tail Extension Summary (Only shows if tail extensions happened)
  # Subset safely whether it's a data.table or data.frame
  tail_ext <- as.data.frame(x$tail)
  tail_ext <- tail_ext[tail_ext$extra_n > 0, ]

  if (nrow(tail_ext) > 0) {
    cat("\n--- Survival Tail Extensions (GPD) ---\n")
    print_tail <- data.frame(
      Outcome     = tail_ext$name,
      Arm         = tail_ext$arm,
      Simulated_N = tail_ext$extra_n,
      stringsAsFactors = FALSE
    )
    print(print_tail, row.names = FALSE, right = FALSE)
  }

  # 4. Helpful footer
  cat("\nNote: Access raw patient data via `$data` or visualize using `plot()`.\n")

  # Invisibly return the object so it doesn't print twice
  invisible(x)


}



#' Summary Method for infuse Objects
#'
#' @param x An object of class \code{infuse}.
#' @param ... Additional arguments (ignored).
#' @export
summary.infuse <- function(x, ...) {
  cat("\nStatistical Summary of Infused Data\n")
  cat(rep("=", 54), "\n", sep = "")
  dt <- copy(x$data)
  #setorder(dt, name, arm, time,prob)

  tt <- dt[, .(mean_value = mean(time, na.rm = TRUE)), by = .(name, arm)]
  tt_wide <- dcast(tt[name %in% x$metadata$name], name ~ arm, value.var = "mean_value")

  setnames(tt_wide,c("name","X","Y"),c("Variable","Mean (X)","Mean (Y)"))
  print(tt_wide[, .(Variable,`Mean (X)`,`Mean (Y)`)], row.names = FALSE,class=FALSE)


  # Check if GPD was used
  if (any(x$tail$extra_n>0)) {
    tt <- x$tail[extra_n>0]
    cat("\nNote: Survival tail extension (GPD) was applied to:",
        paste(unique(tt$name),collapse=", "),"\n")

  }
  cat(rep("=", 54), "\n", sep = "")

  cat("\n--- Point estimates for Favourable and Unfavorable pairs ---\n")
  cat(rep("-", 60), "\n", sep = "")

  #res <- copy(x$metadata)
  #res <- res[,.(name,f,u)]
  res <- unique(x$data[, .(name, f_pe, u_pe)])
  setnames(res,c("name","f_pe","u_pe"),c("Variable","Favorable","Unfavorable"))

  print(as.data.frame(res), row.names = FALSE)

  #cat(rep("-", 45), "\n", sep = "")
  #cat(sprintf("Sample Size: %d Treated (X), %d Control (Y)\n",
  #            length(unique(object$IF_x$ID)), length(unique(object$IF_y$ID))))

}


#' Plot Distribution of Influence Values
#'
#' @description
#' Generates diagnostic plots to visualize the distribution of influence function values
#' across groups. Useful for identifying high-impact observations or outliers.
#'
#' @param x An object of class \code{infuse}.
#' @param vari The influence metric to plot: \code{"f"} (Favorable) or \code{"u"} (Unfavorable).
#' @param ... Additional arguments passed to methods.
#' @export
plotinf <- function(x, vari = "f", ...) {
  UseMethod("plotinf")
}

#' @rdname plotinf
#' @export
plotinf.infuse <- function(x,vari = "f", ...) {

  dt <- x$data[,.(name,arm,f_ifval,u_ifval)]
  setnames(dt,c("f_ifval","u_ifval"),c("f","u"))

  dt <- dt[, IF := get(vari)]
  # Note: Influence for GNNT is usually handled via Delta Method on Net Benefit,
  # but for visualization, showing the 'n' influence is most informative.
  setnames(dt,c("name","arm"),c("Variable","Group"))
  col_name <- ifelse(vari == "f","Favorable","Unfavorable")

  p <- ggplot(dt, aes(x = Group, y = IF, color = Group)) +
    geom_jitter(width = 0.2, alpha = 0.3) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
    stat_summary(fun = "mean", geom = "point", shape = 23, size = 3, fill = "white") +
    facet_wrap(~Variable, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = paste("Influence Distribution"),
         subtitle = "Identifying high-impact observations",
         y = paste("Influence Value (", col_name, ")", sep=""), x = "")

  #return(plotly::ggplotly(p) |> plotly::layout(showlegend = TRUE))
  return(p)
}
