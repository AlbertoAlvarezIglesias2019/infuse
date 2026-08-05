#' Plot Method for realise Class
#'
#' @description
#' Generates a forest-style plot of the Win Statistics (Net Treatment Benefit
#' and Win Ratio) for all outcomes and the NPO.
#'
#' @param x An object of class \code{realise}.
#' @param ... Additional arguments (not used).
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar facet_wrap theme_bw theme element_text geom_vline scale_color_manual labs
#' @importFrom data.table data.table
#' @export
plot.realise <- function(x, ...) {

  # 1. Pre-process data for plotting
  # We filter for only the two metrics requested
  plot_data <- x[Metric %in% c("Net Treatment Benefit", "Win Ratio"), .(
    Variable = factor(Variable, levels = rev(levels(Variable))),
    Metric,
    Estimate,
    Lower,
    Upper,
    Is_NPO = fifelse(Variable == "NPO", "NPO", "Outcome")
  )]

  #plot_data <- x |>
  #  dplyr::filter(Metric %in% c("Net Treatment Benefit", "Win Ratio")) |>
  #  dplyr::mutate(
  #    # Distinct color group for NPO
  #    Is_NPO = ifelse(Variable == "NPO", "NPO", "Outcome"),
  #    # Order Variable so NPO is at the bottom (first level of factor for Y-axis)
  #    Variable = factor(Variable, levels = rev(levels(Variable)))
  #  )

  # 2. Build the ggplot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Estimate, y = Variable, color = Is_NPO)) +
    # Draw the point and the horizontal interval
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = Lower, xmax = Upper), width = 0.3, linewidth = 1) +
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
  p
  #plotly::ggplotly(p)
}


#' Generate a Formatted Summary Table
#'
#' @param x An object to be tabulated.
#' @param ... Additional arguments passed to methods.
#' @export
tab <- function(x, ...) {
  UseMethod("tab")
}


#' Generate a Formatted Summary Table for realise Class
#'
#' @description
#' Creates a publication-quality HTML table using \code{kableExtra}, grouping
#' Net Treatment Benefit, Win Ratio, and GNNT into logical blocks.
#'
#' @param x An object of class \code{realise}.
#' @param ... Additional arguments (not used).
#' @importFrom data.table dcast fcase
#' @importFrom kableExtra kbl kable_styling add_header_above row_spec column_spec footnote
#' @export
tab.realise <- function(x, ...) {

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


  # 1. Filter and Format in one step
  df_formatted <- x[Metric %in% c("Net Treatment Benefit", "Win Ratio", "GNNT"), .(
    Variable,
    Metric,
    Est_Fmt = fcase(
      Metric == "Net Treatment Benefit", sprintf("%.1f%%", Estimate * 100),
      Metric == "Win Ratio",             sprintf("%.2f", Estimate),
      Metric == "GNNT",                  as.character(round(Estimate, 0))
    ),
    CI_Fmt = fcase(
      Metric == "Net Treatment Benefit", fmt_ci(Lower, Upper, "pct"),
      Metric == "Win Ratio",             fmt_ci(Lower, Upper, "wr"),
      Metric == "GNNT",                  fmt_ci(Lower, Upper, "gnnt")
    ),
    PValue
  )]

  # 2. Pivot wider
  # value.var allows multiple columns, and dcast handles the naming pattern automatically
  df_wide <- dcast(
    df_formatted,
    Variable ~ Metric,
    value.var = c("Est_Fmt", "CI_Fmt", "PValue"),
    sep = "_"
  )

  # 3. Organize column order
  # We use backticks because of the spaces in the Metric names
  df_final <- df_wide[, .(
    Variable,
    `Est_Fmt_Net Treatment Benefit`, `CI_Fmt_Net Treatment Benefit`, `PValue_Net Treatment Benefit`,
    `Est_Fmt_Win Ratio`, `CI_Fmt_Win Ratio`, `PValue_Win Ratio`,
    `Est_Fmt_GNNT`, `CI_Fmt_GNNT`, `PValue_GNNT`
  )]


  # Ensure NPO is at the very bottom
  vars <- as.character(df_final$Variable)
  if (length(vars)>1) {
    non_npo <- vars[vars != "NPO"]
    df_final <- df_final[match(c(non_npo, "NPO"), df_final$Variable), ]
  }


  # Rename columns for internal kable display (these are the sub-headers)
  colnames(df_final) <- c("Outcome", rep(c("Estimate", "95% CI", "P-value"), 3))

  # 4. Create the kableExtra table
  kef <- df_final |>
    kableExtra::kbl(booktabs = TRUE, align = "rccccccccc", escape = FALSE) |>
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = FALSE,
      font_size = 12
    ) |>
    # Add the Top-Level blocks
    kableExtra::add_header_above(c(
      " " = 1,
      "Net Treatment Benefit" = 3,
      "Win Ratio" = 3,
      "GNNT" = 3
    ), bold = TRUE, font_size = 14)  |>
    # Vertical Lines (after Outcome, after NTB block, after WR block)
    kableExtra::column_spec(1, border_right = TRUE, bold = TRUE) |>
    kableExtra::column_spec(4, border_right = TRUE) |>
    kableExtra::column_spec(7, border_right = TRUE)

  if (length(vars)>1) {
    kef <- kef   |>
      # Special styling for the NPO row (last row)
      kableExtra::row_spec(
        nrow(df_final),
        bold = TRUE,
        background = "#F9F9F9",
        color = "#D55E00"
      )|>
      # ADDING THE FOOTNOTE HERE
      kableExtra::footnote(
        general = "NPO = Non-prioritised Outcome summary.",
        general_title = "Note: ",
        footnote_as_chunk = TRUE    )
  }

  kef

}



#' @export
knit_print.realise <- function(x, inline = FALSE, ...) {

  dt_clean <- copy(x)

  # Custom formatting function based on Metric type
  fmt <- function(val, metric) {
    fifelse(metric == "GNNT", sprintf("%.0f", val), sprintf("%.4f", val))
  }

  dt_clean[, `:=`(
    Estimate_fmt = fmt(Estimate, Metric),
    SE_fmt       = sprintf("%.4f", Se),
    CI_fmt       = paste0("[", fmt(Lower, Metric), ", ", fmt(Upper, Metric), "]"),
    PValue_fmt   = ifelse(is.na(PValue) | PValue == "<NA>", "—", PValue)
  )]

  # 2. Select final display columns (dropping 'Variable' for section headers)
  display_dt <- dt_clean[, .(
    Metric,
    `Estimate` = Estimate_fmt,
    `SE`       = SE_fmt,
    `95% CI`   = CI_fmt,
    `p-value`  = PValue_fmt
  )]

  # 3. Get row counts for group header indexing
  var_groups <- setNames(
    as.numeric(table(x$Variable)[unique(x$Variable)]),
    unique(x$Variable)
  )

  # 4. Generate kableExtra table
  out <- display_dt |>
    kbl(
      align = c("l", "r", "r", "c", "r"),
      caption = "<b>Table 1:</b> Treatment Effect Metrics Across Clinical Outcomes"
    ) |>
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width = FALSE,
      position = "left",
      font_size = 14
    ) |>
    pack_rows(index = var_groups, label_row_css = "background-color: #e9ecef; color: #212529; font-weight: bold;") |>
    row_spec(0, bold = TRUE, background = "#343a40", color = "white") |>
    column_spec(1, width = "14em") |>
    column_spec(4, width = "12em")

  out
}

#' @export
print.realise <- function(x, ...) {
  print(knit_print.realise(x, ...))
}

