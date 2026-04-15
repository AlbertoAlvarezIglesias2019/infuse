#' Plot Method for reap Class
#'
#' @param x An object of class \code{reap}.
#' @param ... Additional arguments (not used).
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbarh facet_wrap theme_bw theme element_text geom_vline scale_color_manual labs
#' @importFrom plotly ggplotly
#' @export
plot.reap <- function(x, ...) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package \"plotly\" needed. Please install it with install.packages('plotly')")
  }

  # 1. Pre-process data for plotting
  # We filter for only the two metrics requested
  plot_data <- x %>%
    dplyr::filter(Metric %in% c("Net Treatment Benefit", "Win Ratio")) %>%
    dplyr::mutate(
      # Distinct color group for NPO
      Is_NPO = ifelse(Variable == "NPO", "NPO", "Outcome"),
      # Order Variable so NPO is at the bottom (first level of factor for Y-axis)
      Variable = factor(Variable, levels = rev(levels(Variable)))
    )

  # 2. Build the ggplot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Estimate, y = Variable, color = Is_NPO)) +
    # Draw the point and the horizontal interval
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = Lower, xmax = Upper), height = 0.3, linewidth = 1) +
    # Create side-by-side facets for the two metrics
    ggplot2::facet_wrap(~Metric, scales = "free_x") +
    # Add vertical reference lines (0 for NTB, 1 for Win Ratio)
    ggplot2::geom_vline(data = data.frame(Metric = "Net Treatment Benefit", val = 0),
                        ggplot2::aes(xintercept = val), linetype = "dashed", alpha = 0.6) +
    ggplot2::geom_vline(data = data.frame(Metric = "Win Ratio", val = 1),
                        ggplot2::aes(xintercept = val), linetype = "dashed", alpha = 0.6) +
    # Colors and Theme
    ggplot2::scale_color_manual(values = c("NPO" = "#D55E00", "Outcome" = "#0072B2")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 12),           # Overall large font
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 12, face = "bold"), # Facet headers
      legend.position = "none"                          # Hide legend (colors are intuitive)
    ) +
    ggplot2::labs(x = "Estimate (95% Confidence Interval)",
                  y = NULL)

  # 3. Convert to interactive plotly
  plotly::ggplotly(p)
}


#' Generate a Formatted Summary Table
#'
#' @param x An object to be tabulated.
#' @param ... Additional arguments passed to methods.
#' @export
tab <- function(x, ...) {
  UseMethod("tab")
}


#' Generate a Formatted Summary Table for reap Class
#'
#' @param x An object of class \code{reap}.
#' @param ... Additional arguments (not used).
#' @importFrom dplyr filter mutate select case_when
#' @importFrom tidyr pivot_wider
#' @importFrom kableExtra kbl kable_styling add_header_above row_spec
#' @export
tab.reap <- function(x, ...) {

  # 1. Helper to format Confidence Intervals
  fmt_ci <- function(l, u, type = "default") {
    if(type == "pct") {
      paste0("[", sprintf("%.1f%%", l * 100), ", ", sprintf("%.1f%%", u * 100), "]")
    } else if (type == "wr") {
      paste0("[", sprintf("%.2f", l), ", ", sprintf("%.2f", u), "]")
    } else {
      # For GNNT, handling potentially flipped intervals due to NTB ~ 0
      paste0("[", round(l, 0), ", ", round(u, 0), "]")
    }
  }

  # 2. Reshape and Format Data
  # We extract NTB, WR, and GNNT and pivot them into blocks
  df_wide <- x %>%
    dplyr::filter(Metric %in% c("Net Treatment Benefit", "Win Ratio", "GNNT")) %>%
    dplyr::mutate(
      Est_Fmt = dplyr::case_when(
        Metric == "Net Treatment Benefit" ~ sprintf("%.1f%%", Estimate * 100),
        Metric == "Win Ratio" ~ sprintf("%.2f", Estimate),
        Metric == "GNNT" ~ as.character(round(Estimate, 0))
      ),
      CI_Fmt = dplyr::case_when(
        Metric == "Net Treatment Benefit" ~ fmt_ci(Lower, Upper, "pct"),
        Metric == "Win Ratio" ~ fmt_ci(Lower, Upper, "wr"),
        Metric == "GNNT" ~ fmt_ci(Lower, Upper, "gnnt")
      )
    ) %>%
    dplyr::select(Variable, Metric, Est_Fmt, CI_Fmt, PValue) %>%
    tidyr::pivot_wider(
      names_from = Metric,
      values_from = c(Est_Fmt, CI_Fmt, PValue),
      names_glue = "{Metric}_{.value}"
    )

  # 3. Organize column order (NTB Block -> WR Block -> GNNT Block)
  df_final <- df_wide %>%
    dplyr::select(
      Variable,
      `Net Treatment Benefit_Est_Fmt`, `Net Treatment Benefit_CI_Fmt`, `Net Treatment Benefit_PValue`,
      `Win Ratio_Est_Fmt`, `Win Ratio_CI_Fmt`, `Win Ratio_PValue`,
      `GNNT_Est_Fmt`, `GNNT_CI_Fmt`, `GNNT_PValue`
    )

  # Ensure NPO is at the very bottom
  vars <- as.character(df_final$Variable)
  if (length(vars)>1) {
    non_npo <- vars[vars != "NPO"]
    df_final <- df_final[match(c(non_npo, "NPO"), df_final$Variable), ]
  }


  # Rename columns for internal kable display (these are the sub-headers)
  colnames(df_final) <- c("Outcome", rep(c("Estimate", "95% CI", "P-value"), 3))

  # 4. Create the kableExtra table
  kef <- df_final %>%
    kableExtra::kbl(booktabs = TRUE, align = "rccccccccc", escape = FALSE) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = FALSE,
      font_size = 12
    ) %>%
    # Add the Top-Level blocks
    kableExtra::add_header_above(c(
      " " = 1,
      "Net Treatment Benefit" = 3,
      "Win Ratio" = 3,
      "GNNT" = 3
    ), bold = TRUE, font_size = 14)  %>%
    # Vertical Lines (after Outcome, after NTB block, after WR block)
    kableExtra::column_spec(1, border_right = TRUE, bold = TRUE) %>%
    kableExtra::column_spec(4, border_right = TRUE) %>%
    kableExtra::column_spec(7, border_right = TRUE)

  if (length(vars)>1) {
    kef <- kef   %>%
      # Special styling for the NPO row (last row)
      kableExtra::row_spec(
        nrow(df_final),
        bold = TRUE,
        background = "#F9F9F9",
        color = "#D55E00"
      )%>%
      # ADDING THE FOOTNOTE HERE
      kableExtra::footnote(
        general = "NPO = Non-prioritised Outcome summary.",
        general_title = "Note: ",
        footnote_as_chunk = TRUE    )
  }

  kef

}


