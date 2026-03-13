#' Plot Canonical Tocky Ordination with Continuous Gradient Legend
#'
#' Visualizes the constrained ordination space in a 3-panel layout.
#' Panel 1: Axis 1 vs 2. Panel 2: Axis 2 vs 3. Panel 3: Continuous Gradient Legend.
#'
#' @param tocky_res List object returned by `CanonicalTockySeq` containing `$cell_scores` and `$biplot`.
#' @param gradient_res List object returned by `GradientTockySeq` containing `$angle` Gradient Tocky Time scores (0-90).
#' @param alpha_level Numeric transparency level for cells (0-1). Default is 0.2.
#'
#' @importFrom grDevices colorRampPalette rgb as.raster col2rgb
#' @importFrom graphics plot abline arrows text par rasterImage rect title layout segments
#' @export
plotCanonicalTocky <- function(tocky_res, gradient_res, alpha_level = 0.2) {
  
  tocky_time <- gradient_res$angle
   layout(matrix(c(1, 2, 3), nrow = 1, ncol = 3), widths = c(1, 1, 0.8))
  
  n_colors <- 100
  rbPal <- colorRampPalette(c("blue", "purple", "red"))
  color_map <- rbPal(n_colors)
  
  rgb_vals <- col2rgb(color_map)
  transparent_map <- rgb(rgb_vals[1,], rgb_vals[2,], rgb_vals[3,], maxColorValue = 255, alpha = alpha_level * 255)
  
  grey_transparent <- rgb(0.8, 0.8, 0.8, alpha = alpha_level)
  cell_colors <- rep(grey_transparent, nrow(tocky_res$cell_scores))
  
  valid_idx <- !is.na(tocky_time)
  
  if (sum(valid_idx) > 0) {
    vals <- pmax(0, pmin(90, tocky_time[valid_idx]))
    idx <- floor((vals / 90) * (n_colors - 1)) + 1
    cell_colors[valid_idx] <- transparent_map[idx]
  }

  #PANEL 1: Axis 1 vs 2
  par(mar = c(5, 4, 4, 1))
  
  coords_12 <- tocky_res$cell_scores[, c(1, 2)]
  vecs_12   <- tocky_res$biplot[, c(1, 2)]
  
  scale_12 <- (max(abs(coords_12)) / max(abs(vecs_12))) * 0.5
  vecs_12_sc <- vecs_12 * scale_12
  
  plot(coords_12, col = cell_colors, pch = 16, cex = 1.0, # Bigger points
       xlab = "Axis 1", ylab = "Axis 2", main = "Projection: 1 vs 2",
       las = 1, bty = "l")
  abline(h = 0, v = 0, lty = 2, col = "grey80")
  arrows(0, 0, vecs_12_sc[,1], vecs_12_sc[,2], length = 0.1, lwd = 2, col = "black")
  text(vecs_12_sc[,1]*1.2, vecs_12_sc[,2]*1.2, rownames(vecs_12),
       col="black", font=2, cex=1.0)

  #PANEL 2: Axis 2 vs 3
  coords_23 <- tocky_res$cell_scores[, c(2, 3)]
  vecs_23   <- tocky_res$biplot[, c(2, 3)]
  
  scale_23 <- (max(abs(coords_23)) / max(abs(vecs_23))) * 0.5
  vecs_23_sc <- vecs_23 * scale_23
  
  plot(coords_23, col = cell_colors, pch = 16, cex = 1.0, # Bigger points
       xlab = "Axis 2", ylab = "Axis 3", main = "Projection: 2 vs 3",
       las = 1, bty = "l")
  abline(h = 0, v = 0, lty = 2, col = "grey80")
  arrows(0, 0, vecs_23_sc[,1], vecs_23_sc[,2], length = 0.1, lwd = 2, col = "black")
  text(vecs_23_sc[,1]*1.2, vecs_23_sc[,2]*1.2, rownames(vecs_23),
       col="black", font=2, cex=1.0)

  #Legend
  par(mar = c(5, 0, 4, 2))
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
  title("Tocky Time", line = 1)
  
  legend_raster <- as.raster(matrix(rev(transparent_map), ncol=1))
  
  rasterImage(legend_raster, 0.2, 0.2, 0.3, 0.8)
  
  tick_y <- c(0.2, 0.5, 0.8)
  tick_lab <- c("0 (New)", "45 (Persistent)", "90 (Arrested)")
  tick_col <- c("blue", "purple", "red") # Text stays solid for readability
  
  segments(0.3, tick_y, 0.35, tick_y, lwd = 1)
  text(0.38, tick_y, labels = tick_lab, pos = 4, col = tick_col, cex = 1.1, font=2)
  
  rect(0.2, 0.05, 0.3, 0.1, col = grey_transparent, border = NA)
  text(0.38, 0.075, labels = "Timer Neg", pos = 4, col = "grey50", cex = 1.1)
}



#' Plot Gene Expression Dynamics by Group (Fixed Color Mapping & Robust Limits)
#'
#' Visualizes gene expression trends along the Gradient Tocky Time axis (0-90 degrees).
#' Supports both data frame output from GradientTocky_Slerp and simple numeric vectors.
#'
#' @param expression_data Sparse matrix or data frame (Genes x Cells).
#' @param gradient_res Either the data frame returned by `GradientTocky_Slerp` OR a named numeric vector of angles.
#' @param gene Character string. The gene symbol to plot.
#' @param groups Optional Character or Factor vector defining cell groups (e.g., "WT", "KO").
#' @param group_cols Optional vector of colors. Can be named (c("WT"="black")) or unnamed.
#' @param span Numeric. Smoothing span for LOESS (default 0.8).
#' @param jitter_amount Numeric. Amount of vertical jitter for points (default 0.2).
#' @param pt_alpha Numeric. Transparency of single cell points (0-1).
#' @param m Numeric. Magnification factor for y-max.
#'
#' @importFrom graphics plot rect lines legend grid par
#' @importFrom stats loess predict
#' @importFrom grDevices rgb adjustcolor rainbow
#' @export
plotGeneDynamics <- function(expression_data, gradient_res, gene,
                             groups = NULL,
                             group_cols = NULL,
                             span = 0.8, jitter_amount = 0.2,
                             pt_alpha = 0.2, m = 1.15) {
  
  if (!gene %in% rownames(expression_data)) stop(paste("Gene", gene, "not found in expression matrix."))
  
  if (is.data.frame(gradient_res)) {
    if (!"angle" %in% colnames(gradient_res)) stop("gradient_res data frame must contain an 'angle' column.")
    tocky_time <- gradient_res$angle
    names(tocky_time) <- rownames(gradient_res)
  } else {
    tocky_time <- gradient_res
  }
  
  common_cells <- intersect(colnames(expression_data), names(tocky_time))
  
  if (length(common_cells) == 0) {
    warning("No common cells found between expression matrix and gradient data.")
    return(NULL)
  }
  
  if (!is.null(groups)) {
    if (!is.null(names(groups))) {
      common_cells <- intersect(common_cells, names(groups))
      groups <- groups[common_cells]
    } else if (length(groups) == length(tocky_time)) {
      names(groups) <- names(tocky_time)
      groups <- groups[common_cells]
    } else {
      stop("Length of 'groups' does not match data and no names were provided.")
    }
  }
  
  expr_values <- as.numeric(expression_data[gene, common_cells])
  times       <- as.numeric(tocky_time[common_cells])
  
  valid_idx <- !is.na(times) & !is.na(expr_values)
  if (!is.null(groups)) valid_idx <- valid_idx & !is.na(groups)
  
  plot_x <- times[valid_idx]
  plot_y <- expr_values[valid_idx]
  
  if (length(plot_x) == 0) {
    plot(1, 1, type="n", xlim=c(0, 90), ylim=c(0, 1),
         xlab="Gradient Tocky Time", ylab="Expression",
         main=paste(gene, "(No Valid Cells)"))
    text(45, 0.5, "No valid cells for this gene/group", col="red")
    return(NULL)
  }
  
  if (!is.null(groups)) {
    if (is.factor(groups)) {
      plot_group <- droplevels(groups[valid_idx])
    } else {
      plot_group <- as.factor(groups[valid_idx])
    }
    unique_groups <- levels(plot_group)
  } else {
    plot_group <- NULL
    unique_groups <- NULL
  }
  
   if (is.null(plot_group)) {
    pt_colors <- rgb(0.2, 0.2, 0.2, pt_alpha)
    line_colors <- "black"
    unique_groups <- "Cells"
  } else {
    if (is.null(group_cols)) {
      default_pal <- c("red", "blue", "forestgreen", "orange", "purple")
      if(length(unique_groups) > length(default_pal)) {
        group_cols <- rainbow(length(unique_groups))
      } else {
        group_cols <- default_pal[1:length(unique_groups)]
      }
      names(group_cols) <- unique_groups
    } else {
      if (is.null(names(group_cols))) {
        if(length(group_cols) < length(unique_groups)){
          stop("Not enough colors provided for the number of groups.")
        }
        names(group_cols)[1:length(unique_groups)] <- unique_groups
      } else {
        missing_grps <- setdiff(unique_groups, names(group_cols))
        if(length(missing_grps) > 0) {
          stop(paste("Missing colors for groups:", paste(missing_grps, collapse=", ")))
        }
      }
    }
    pt_colors_base <- group_cols[as.character(plot_group)]
    pt_colors <- adjustcolor(pt_colors_base, alpha.f = pt_alpha)
  }
  
  plot_y_jittered <- jitter(plot_y, amount = jitter_amount)
  plot_y_jittered[plot_y_jittered < 0] <- 0
  
  ymax_raw <- max(plot_y_jittered, na.rm=TRUE)
  if (!is.finite(ymax_raw) || ymax_raw == 0) {
    ymax <- 1.0
  } else {
    ymax <- ymax_raw * m
  }
  
  plot(plot_x, plot_y_jittered,
       pch = 16, cex = 0.6, col = pt_colors,
       xlab = "Gradient Tocky Time (Angle)",
       ylab = paste(gene, "Expression (Log)"),
       main = gene,
       las = 1, xlim = c(0, 90),
       ylim = c(0, ymax))
  
  rect(0, -10, 30, ymax*2,  col = rgb(0, 0, 1, 0.05), border = NA)
  rect(30, -10, 60, ymax*2, col = rgb(0.5, 0, 0.5, 0.05), border = NA)
  rect(60, -10, 90, ymax*2, col = rgb(1, 0, 0, 0.05), border = NA)
  
  legend_text <- c()
  legend_cols_vec <- c()
  
  draw_curve <- function(x_vec, y_vec, col_hex, label_txt) {
    if (length(x_vec) < 10) return(NULL)
    
    loess_fit <- tryCatch(loess(y_vec ~ x_vec, span = span), error = function(e) NULL)
    if (is.null(loess_fit)) return(NULL)
    
    x_seq <- seq(min(x_vec), max(x_vec), length.out = 100)
    y_pred <- predict(loess_fit, newdata = data.frame(x_vec = x_seq))
    y_pred[y_pred < 0] <- 0
    
    lines(x_seq, y_pred, col = col_hex, lwd = 3)
    
    pct <- round(sum(y_vec > 0) / length(y_vec) * 100, 0)
    return(paste0(label_txt, " (", pct, "%)"))
  }
  
  if (is.null(plot_group)) {
    lab <- draw_curve(plot_x, plot_y, line_colors, "Cells")
    if(!is.null(lab)) legend("topleft", legend=lab, col="black", lwd=3, bty="n")
  } else {
    for (grp in unique_groups) {
      idx <- plot_group == grp
      grp_col <- group_cols[grp]
      
      lab <- draw_curve(plot_x[idx], plot_y[idx], grp_col, grp)
      
      if (!is.null(lab)) {
        legend_text <- c(legend_text, lab)
        legend_cols_vec <- c(legend_cols_vec, grp_col)
      }
    }
    
    if(length(legend_text) > 0) {
      legend("topleft", legend = legend_text, col = legend_cols_vec,
             lwd = 3, bty = "n", cex = 0.8)
    }
  }
  grid()
}


#' Plot Pseudotime Heatmap (Ordered by Peak Timing)
#'
#' Visualizes gene expression along the Tocky trajectory.
#' Supports ordering genes by their "Peak Time" (Cascade) or Hierarchical Clustering.
#'
#' @param object Seurat object or expression matrix.
#' @param gradient_res Data frame from GradientTocky or named vector of angles.
#' @param genes Character vector of genes to plot.
#' @param n_bins Integer. Number of bins (default 100).
#' @param ordering_method Character. "peak" (sort by max expression time) or "cluster" (hclust).
#' @param span Numeric. Smoothing span (default 0.5).
#' @param scale_rows Logical. Z-score rows (default TRUE).
#'
#' @importFrom stats loess predict hclust dist sd complete.cases
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom Seurat GetAssayData
#' @export
PlotTockyHeatmap <- function(object, gradient_res, genes,
                             n_bins = 100,
                             ordering_method = c("peak", "cluster"),
                             span = 0.5, scale_rows = TRUE) {
  
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Please install 'pheatmap'.")
  ordering_method <- match.arg(ordering_method)

  if (is.data.frame(gradient_res)) {
    tocky_time <- gradient_res$angle
    names(tocky_time) <- rownames(gradient_res)
  } else {
    tocky_time <- gradient_res
  }
  
  if (inherits(object, "Seurat")) {
    counts <- Seurat::GetAssayData(object, layer = "data")
  } else {
    counts <- object
  }
  
  common_cells <- intersect(colnames(counts), names(tocky_time))
  valid_idx <- !is.na(tocky_time[common_cells])
  cells_use <- common_cells[valid_idx]
  
  if (length(cells_use) < 10) stop("Too few cells.")
  
  time_vec <- tocky_time[cells_use]
  valid_genes <- intersect(genes, rownames(counts))
  if (length(valid_genes) == 0) stop("No valid genes found.")
  
  expr_mat <- as.matrix(counts[valid_genes, cells_use])
  
  breaks <- seq(0, 90, length.out = n_bins + 1)
  bin_centers <- (breaks[1:n_bins] + breaks[2:(n_bins+1)]) / 2
  cell_bins <- cut(time_vec, breaks = breaks, labels = FALSE, include.lowest = TRUE)
  
  binned_mat <- matrix(NA, nrow = nrow(expr_mat), ncol = n_bins)
  rownames(binned_mat) <- rownames(expr_mat)
  colnames(binned_mat) <- round(bin_centers, 1)
  
  cat(sprintf("Binning %d cells into %d angular intervals...\n", length(cells_use), n_bins))
  
  for (i in 1:n_bins) {
    cells_in_bin <- which(cell_bins == i)
    if (length(cells_in_bin) > 0) {
      if (length(cells_in_bin) == 1) {
        binned_mat[, i] <- expr_mat[, cells_in_bin]
      } else {
        binned_mat[, i] <- rowMeans(expr_mat[, cells_in_bin])
      }
    }
  }
  
  smoothed_mat <- binned_mat
  if (span > 0) {
    cat("Smoothing binned trajectories...\n")
    x_grid <- 1:n_bins
    for (i in 1:nrow(binned_mat)) {
      y_vals <- binned_mat[i, ]
      valid_bins <- !is.na(y_vals)
      if (sum(valid_bins) < 5) next
      fit <- tryCatch(loess(y_vals[valid_bins] ~ x_grid[valid_bins], span = span),
                      error = function(e) NULL)
      if (!is.null(fit)) {
        smoothed_mat[i, ] <- predict(fit, newdata = x_grid)
      }
    }
  }
  
  smoothed_mat <- smoothed_mat[complete.cases(smoothed_mat), , drop = FALSE]
  if (nrow(smoothed_mat) < 2) {
    warning("Fewer than 2 genes passed smoothing. Skipping plot.")
    return(NULL)
  }

  plot_mat <- smoothed_mat
  if (scale_rows) {
    row_means <- rowMeans(plot_mat)
    row_sds   <- apply(plot_mat, 1, sd)
    plot_mat  <- (plot_mat - row_means) / replace(row_sds, row_sds==0, 1)
    plot_mat[plot_mat > 3]  <- 3
    plot_mat[plot_mat < -3] <- -3
  }
  
  if (ordering_method == "peak") {
    peak_col <- apply(plot_mat, 1, which.max)
    plot_mat <- plot_mat[order(peak_col, decreasing = FALSE), ]
    cluster_rows_flag <- FALSE
    main_title <- "Tocky Cascade (Ordered by Peak Time)"
  } else {
    cluster_rows_flag <- TRUE
    main_title <- "Tocky Heatmap (Hierarchical Clustering)"
  }

  heat_colors <- colorRampPalette(c("#7b3294", "#f7f7f7", "#d95f0e"))(100)

  pheatmap::pheatmap(plot_mat,
                     cluster_rows = cluster_rows_flag,
                     cluster_cols = FALSE,
                     show_colnames = FALSE,
                     main = main_title,
                     color = heat_colors,
                     border_color = NA,
                     fontsize_row = 8)
}
