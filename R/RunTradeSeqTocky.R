#' Run TradeSeq Differential Dynamics Analysis on Tocky Trajectory
#'
#' A wrapper function to perform differential dynamic expression testing using tradeSeq.
#' Supports Seurat v5+ (uses 'layer') and matrix inputs.
#'
#' @param object A Seurat object OR a sparse expression matrix (Genes x Cells).
#' @param gradient_res The data frame returned by `GradientTocky_Slerp` (must contain `$angle`).
#' @param group_by If `object` is Seurat: Character string of the metadata column name (e.g., "treatment").
#'                 If `object` is Matrix: A named factor/character vector of groups matching column names.
#' @param genes Character vector of genes to test (e.g., output of `SelectTockyGenes`).
#' @param n_knots Integer. Number of knots for the GAM (default 5).
#' @param n_cores Integer. Number of cores for parallel processing (default 1).
#'
#' @return A list containing:
#'   \item{results}{Data frame of differential dynamics results (p-values, FDR).}
#'   \item{gam_fit}{The raw tradeSeq GAM object (for downstream plotting).}
#'
#' @importFrom tradeSeq fitGAM conditionTest
#' @importFrom BiocParallel MulticoreParam SerialParam SnowParam
#' @importFrom stats p.adjust
#' @importFrom methods is
#' @export
RunTradeSeqTocky <- function(object, gradient_res, group_by, genes,
                             n_knots = 5, n_cores = 1) {
  
  if (!requireNamespace("tradeSeq", quietly = TRUE)) {
    stop("Package 'tradeSeq' is required. Please install it.")
  }
  
  if (methods::is(object, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required for Seurat object input.")
    
    tryCatch({
      counts_mat <- Seurat::GetAssayData(object, layer = "counts")
    }, error = function(e) {
      counts_mat <- Seurat::GetAssayData(object, slot = "counts")
    })
    
     if (!group_by %in% colnames(object@meta.data)) {
      stop(paste("Column", group_by, "not found in Seurat metadata."))
    }
    group_vec <- object@meta.data[[group_by]]
    names(group_vec) <- colnames(object)
    
  } else {
    counts_mat <- object
    group_vec <- group_by
    
    if (is.null(names(group_vec))) {
      if (length(group_vec) == ncol(counts_mat)) {
        names(group_vec) <- colnames(counts_mat)
        warning("group_by vector was unnamed. Assumed order matches expression matrix columns.")
      } else {
        stop("group_by must be a named vector or match the number of matrix columns.")
      }
    }
  }
  
  if (is.data.frame(gradient_res)) {
    if (!"angle" %in% colnames(gradient_res)) stop("gradient_res must contain '$angle'.")
    tocky_vals <- gradient_res$angle
    names(tocky_vals) <- rownames(gradient_res)
  } else {
    tocky_vals <- gradient_res
  }

  common_cells <- intersect(colnames(counts_mat), names(tocky_vals))
  common_cells <- intersect(common_cells, names(group_vec))
  
  if (length(common_cells) == 0) stop("No common cells found between expression, gradient, and groups.")
  
    tocky_subset <- tocky_vals[common_cells]
    group_subset <- group_vec[common_cells]
    valid_idx <- !is.na(tocky_subset) & !is.na(group_subset)
  
  if (sum(valid_idx) < 10) stop("Fewer than 10 valid (Timer Positive) cells remain. Cannot run GAM.")
  
  final_cells <- common_cells[valid_idx]
  
  valid_genes <- intersect(genes, rownames(counts_mat))
  if(length(valid_genes) < length(genes)) {
    warning(paste("Dropped", length(genes) - length(valid_genes), "genes not found in expression matrix."))
  }
  
  y_counts <- as.matrix(counts_mat[valid_genes, final_cells])
  t_time   <- tocky_subset[final_cells]
  f_groups <- as.factor(group_vec[final_cells])
  
  cat(sprintf("Running tradeSeq on %d genes across %d cells (%d groups)...\n",
              nrow(y_counts), ncol(y_counts), length(levels(f_groups))))

  pseudotime_mat <- matrix(t_time, ncol = 1)
  rownames(pseudotime_mat) <- final_cells
  colnames(pseudotime_mat) <- "Lineage1"
  
  weights_mat <- matrix(1, nrow = length(final_cells), ncol = 1)
  rownames(weights_mat) <- final_cells
  colnames(weights_mat) <- "Lineage1"
  
  if (n_cores > 1) {
    if (.Platform$OS.type == "windows") {
      bpparam <- BiocParallel::SnowParam(workers = n_cores)
    } else {
      bpparam <- BiocParallel::MulticoreParam(workers = n_cores)
    }
    parallel_flag <- TRUE
  } else {
    bpparam <- BiocParallel::SerialParam()
    parallel_flag <- FALSE
  }
  
  gam_fit <- tradeSeq::fitGAM(counts = y_counts,
                              pseudotime = pseudotime_mat,
                              cellWeights = weights_mat,
                              conditions = f_groups,
                              nknots = n_knots,
                              parallel = parallel_flag,
                              BPPARAM = bpparam,
                              verbose = TRUE)
  
  # 5. Run Condition Test (Differential Dynamics)
  cat("Testing for differential dynamics (Condition Test)...\n")
  cond_res <- tradeSeq::conditionTest(gam_fit, global = TRUE)
  
  cond_res$padj <- p.adjust(cond_res$pvalue, method = "fdr")
  cond_res <- cond_res[order(cond_res$padj), ]
  
  return(list(
    results = cond_res,
    gam_fit = gam_fit
  ))
}


#' Select Genes Dynamic along Tocky Time
#'
#' Rapidly filters for genes that show significant variation along the
#' Gradient Tocky Time axis using a linear spline regression.
#' Automatically excludes cells with NA angles (Timer Negative).
#'
#' @param expression_data Sparse matrix or data frame (Genes x Cells).
#' @param gradient_res Either the data frame returned by `GradientTocky_Slerp` OR a named numeric vector of angles.
#' @param top_n Integer. Number of genes to return (default 2000).
#' @param min_expr Minimum proportion of cells expressing the gene (0-1). Default 0.05.
#'
#' @return A character vector of gene names sorted by dynamic score (R-squared).
#' @importFrom splines ns
#' @importFrom stats lm
#' @importFrom utils head
#' @importFrom Matrix rowMeans
#' @export
SelectTockyGenes <- function(expression_data, gradient_res, top_n = 2000, min_expr = 0.05) {
  
  if (is.data.frame(gradient_res)) {
    if (!"angle" %in% colnames(gradient_res)) stop("gradient_res data frame must contain an 'angle' column.")
    tocky_time <- gradient_res$angle
    names(tocky_time) <- rownames(gradient_res)
  } else {
    tocky_time <- gradient_res
  }
  
  if (!requireNamespace("splines", quietly = TRUE)) {
    stop("Package 'splines' is required for this function.")
  }
  
  common_cells <- intersect(colnames(expression_data), names(tocky_time))
  
  if(length(common_cells) < 10) {
    stop("Fewer than 10 common cells found between expression data and gradient results.")
  }
  
  sub_time_full <- tocky_time[common_cells]
  
  valid_t_idx <- !is.na(sub_time_full)
  valid_cells <- common_cells[valid_t_idx]
  sub_time    <- sub_time_full[valid_t_idx]
  
  if(length(sub_time) < 10) {
    stop("After removing Timer Negative (NA) cells, fewer than 10 cells remain. Cannot fit splines.")
  }
  
   sub_expr <- expression_data[, valid_cells]
  
  n_cells <- length(valid_cells)
  
  if (inherits(sub_expr, "sparseMatrix")) {
    gene_freq <- Matrix::rowMeans(sub_expr > 0)
  } else {
    gene_freq <- rowMeans(sub_expr > 0)
  }
  
  valid_genes <- names(gene_freq)[gene_freq >= min_expr]
  
  if(length(valid_genes) == 0) {
    stop("No genes passed the expression threshold.")
  }
  
  cat(sprintf("  Initial screen: %d genes passed expression threshold (%.1f%%).\n",
              length(valid_genes), min_expr*100))
  
  cat(sprintf("  Scoring dynamics for %d genes...\n", length(valid_genes)))
  
  time_basis <- splines::ns(sub_time, df = 3)
  
  get_score <- function(g) {
    val <- as.numeric(sub_expr[g, ])
    tryCatch({
      fit <- lm(val ~ time_basis)
      summary(fit)$r.squared
    }, error = function(e) 0)
  }
  
  scores <- vapply(valid_genes, get_score, numeric(1))
  
  scores[is.na(scores)] <- 0
  sorted_genes <- names(sort(scores, decreasing = TRUE))
  
  n_select <- min(top_n, length(sorted_genes))
  selected <- head(sorted_genes, n_select)
  
  cat(sprintf("  Selected top %d genes.\n", length(selected)))
  
  return(selected)
}
