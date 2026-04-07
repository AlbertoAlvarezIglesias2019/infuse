library(plotly)

#' Print Method for harvest Objects
#' @param x An object of class \code{harvest}.
#' @param ... Additional arguments (ignored).
#' @export
print.harvest <- function(x, ...) {
  cat("\n--- infuse: Harvested Win/Loss Results ---\n")
  cat(sprintf("Outcomes processed: %d\n\n", nrow(x$theta)))

  # Prepare a clean table for printing
  out_tab <- x$theta[, c("Variable", "f", "u", "n","w")]
  out_tab$w <- exp(out_tab$w)
  colnames(out_tab) <- c("Outcome", "Wins(f)", "Losses(u)", "NetBenefit(n)", "WinRatio(w)")

  print(out_tab, row.names = FALSE)

  #cat("\nNote: Win Ratio (w) is stored as log(f/u) for internal inference.\n")
  #cat("Run summary() for IF variance or reap() for p-values/CIs.\n")
  invisible(x)
}


#' Summary Method for harvest Objects
#' @param object An object of class \code{harvest}.
#' @param ... Additional arguments (ignored).
#' @export
summary.harvest <- function(object, ...) {
  cat("\nDetailed Harvest Summary (Point Estimates)\n")
  cat(rep("-", 45), "\n", sep = "")

  # Calculate IF variance per variable for Net Benefit (n)
  #var_n_x <- aggregate(n ~ Variable, data = object$IF_x, FUN = var)
  #var_n_y <- aggregate(n ~ Variable, data = object$IF_y, FUN = var)

  # Merge with theta
  res <- object$theta
  names(res)[4] <- "log(w)"
  res[,1:4] <- round(object$theta[,1:4],5)
  #res$Var_IF_x <- round(var_n_x$n, 5)
  #res$Var_IF_y <- round(var_n_y$n, 5)

  print(res, row.names = FALSE)

  cat(rep("-", 45), "\n", sep = "")
  cat(sprintf("Sample Size: %d Treated (X), %d Control (Y)\n",
              length(unique(object$IF_x$ID)), length(unique(object$IF_y$ID))))
}




#' Plot Method for harvest Objects
#' @param x An object of class \code{harvest}.
#' @param metric One of \code{"netb"} (Net Benefit), \code{"wr"} (Win Ratio), or \code{"gnnt"} (Generalized NNT).
#' @param type One of \code{"estimate"} (Bar chart of results) or \code{"influence"} (Diagnostic dot plot).
#' @param ... Additional arguments (ignored).
#' @import ggplot2
#' @importFrom plotly ggplotly layout
#' @export
plot.harvest <- function(x, metric = c("ntb", "wr", "gnnt"), type = c("estimate", "influence"), ...) {

  type <- match.arg(type)
  metric <- match.arg(metric)

  # --- 1. Map Metric to Column and Label ---
  # 'col' refers to the column name in theta/IF dataframes
  target <- switch(metric,
                   "ntb" = list(col = "n", lab = "Net Treatment Benefit", line = 0),
                   "wr"   = list(col = "w", lab = "Log ( Win Ratio )", line = 0),
                   "gnnt" = list(col = "n", lab = "GNNT", line = 0)) # GNNT derived from n

  # --- 2. Logic for Estimates (Bar Chart) ---
  if (type == "estimate") {
    plot_df <- x$theta

    # Special handling for GNNT (1/Net Benefit)
    if (metric == "gnnt") {
      plot_df$val <- 1 / plot_df$n
    } else {
      plot_df$val <- plot_df[[target$col]]
    }

    p <- ggplot(plot_df, aes(x = Variable, y = val, fill = val > target$line)) +
      geom_bar(stat = "identity", alpha = 0.7, color = "black") +
      geom_hline(yintercept = target$line, linetype = "dashed")

    #if (metric=="wr") p <- p + scale_y_continuous(trans = "exp")

    p <- p +
      coord_flip() +
      scale_fill_manual(values = c("TRUE" = "#2c7fb8", "FALSE" = "#f03b20"),
                        labels = c("TRUE" = "Favorable", "FALSE" = "Unfavorable")) +
      theme_minimal() +
      labs(title = paste(target$lab, "by Outcome"),
           subtitle = "Point estimates per variable",
           x = "", y = target$lab, fill = "Direction")


    return(plotly::ggplotly(p) %>% plotly::layout(showlegend = TRUE))

    # --- 3. Logic for Influence (Dot Plot) ---
  } else {
    # Extract the correct influence column
    col_name <- target$col

    df_if <- rbind(
      data.frame(IF = x$IF_x[[col_name]], Group = "Treated (X)", Variable = x$IF_x$Variable),
      data.frame(IF = x$IF_y[[col_name]], Group = "Control (Y)", Variable = x$IF_y$Variable)
    )

    # Note: Influence for GNNT is usually handled via Delta Method on Net Benefit,
    # but for visualization, showing the 'n' influence is most informative.

    p <- ggplot(df_if, aes(x = Group, y = IF, color = Group)) +
      geom_jitter(width = 0.2, alpha = 0.3) +
      geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
      stat_summary(fun = "mean", geom = "point", shape = 23, size = 3, fill = "white") +
      facet_wrap(~Variable, scales = "free_y") +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(title = paste(target$lab, "Influence Distribution"),
           subtitle = "Identifying high-impact observations",
           y = paste("Influence Value (", col_name, ")", sep=""), x = "")

    return(plotly::ggplotly(p) %>% plotly::layout(showlegend = TRUE))
  }
}



