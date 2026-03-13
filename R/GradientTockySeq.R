#' Gradient Tocky Sequence (Slerp Model)
#'
#' Projects cells onto a piecewise Slerp manifold defined by the landmarks B -> BR -> R.
#' Cells with negative correlations to all three landmarks are classified as "Timer Negative" (NA),
#' as they lie in the ordination space opposite to the Tocky trajectory.
#'
#' @param res List object returned by `CanonicalTockySeq` (must contain `$cell_scores` coordinates).
#' @param B Vector for the 'New' (Blue) landmark.
#' @param BR Vector for the 'Persistent' (Blue-Red) landmark.
#' @param R Vector for the 'Arrested' (Red) landmark.
#' @param filter_negative Logical. If TRUE, assigns NA to cells with negative projection to all three landmarks.
#'
#' @return A data frame containing:
#'   \item{angle}{Tocky Time (0-90)}
#'   \item{intensity}{Raw Euclidean norm (distance from origin)}
#'   \item{norm_intensity}{Intensity normalized against the Slerp trajectory}
#'   \item{similarity}{Cosine similarity to the manifold}
#' @export
#' @importFrom stats optimize

GradientTockySeq <- function(res, B, BR, R, filter_negative = TRUE) {
  n <- length(B)
  cell_mat <- as.matrix(res$cell_scores[, 1:n])
  n_cells <- nrow(cell_mat)
  
  is_timer_neg <- rep(FALSE, n_cells)
  if (filter_negative) {
    is_timer_neg <- (cell_mat %*% B < 0) & (cell_mat %*% BR < 0) & (cell_mat %*% R < 0)
  }

  norm_B <- get_norm(B); unit_B <- B / norm_B
  norm_BR <- get_norm(BR); unit_BR <- BR / norm_BR
  norm_R <- get_norm(R); unit_R <- R / norm_R
  
  omega1 <- acos(pmax(pmin(sum(unit_B * unit_BR), 1), -1))
  omega2 <- acos(pmax(pmin(sum(unit_BR * unit_R), 1), -1))

  t_values <- rep(NA, n_cells)
  sim_values <- rep(NA, n_cells)
  valid_idx <- which(!is_timer_neg)
  
  if (length(valid_idx) > 0) {
    subset_res <- apply(cell_mat[valid_idx, , drop=FALSE], 1, get_best_t_for_cell,
                        B = B, BR = BR, R = R, omega1 = omega1, omega2 = omega2)
    
    t_values[valid_idx] <- subset_res["t", ]
    sim_values[valid_idx] <- subset_res["sim", ]
  }

  cell_norms <- sqrt(rowSums(cell_mat^2))
  expected_mags <- sapply(t_values, function(t) {
    if (is.na(t)) return(NA)
    ref <- if(t <= 0.5) fast_slerp(B, BR, t*2, omega1) else fast_slerp(BR, R, (t-0.5)*2, omega2)
    get_norm(ref)
  })
  
  return(data.frame(
    angle = t_values * 90,
    intensity = cell_norms,
    norm_intensity = cell_norms / expected_mags,
    similarity = sim_values
  ))
}
