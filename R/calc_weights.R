calc_weights <- function(uni, w_vec, raw_vals) {
  ord <- order(raw_vals)
  cum_w <- c_ecdf_values(uni, w_vec, raw_vals[ord])
  # Calculate differences to get individual point weights
  w <- cum_w - dplyr::lag(cum_w, default = 0)
  return(w[order(ord)])
}
