#' Canonical Redundancy Analysis for Tocky Differentiation using single-cell RNA-sequencing data (CanonicalTockySeq)
#'
#' Performs a specialized Canonical Redundancy Analysis (RDA) designed for single-cell
#' Tocky data. In this architecture, genes are treated as "sites" and cells as "species",
#' constrained by Tocky-defined developmental stages (Blue, Blue-Red, Red).
#'
#' @details
#' Implements the two-step RDA procedure described by Legendre & Legendre (2012):
#' 1. Multivariate regression of Expression (X) on Constraints (Z).
#' 2. PCA of the fitted values to derive canonical axes.
#'
#' @param X Numeric matrix of gene expression (Genes x Cells).
#' @param Z Numeric matrix of explanatory variables (Genes x Signatures). Typically
#'   contains Tocky-specific gene signatures (e.g., Blue, Blue-Red, Red).
#'
#' @return A list containing:
#' \describe{
#'   \item{expression_scores}{Data frame of Gene loadings (Genes x k).}
#'   \item{fitted_cell_scores}{Data frame of Fitted Cell scores (Cells x k).}
#'   \item{cell_scores}{Data frame of Cell scores (Cells x k). Use this for downstream GradientTockySeq.}
#'   \item{biplot}{Matrix of biplot vectors for the constraint variables (Z).}
#' }
#'
#' @importFrom irlba irlba
#' @export
CanonicalTockySeq <- function(X, Z) {
  
  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Package 'irlba' is required for fast SVD. Please install it.")
  }
  
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  
  Z_scaled <- scale(Z, center = FALSE, scale = TRUE)
  S <- scale(X)
  
  inv_ZtZ <- solve(t(Z_scaled) %*% Z_scaled)
  projection_coeff <- Z_scaled %*% inv_ZtZ
  ZtS <- t(Z_scaled) %*% S
  S_star <- projection_coeff %*% ZtS
  
  cat("Normalization and projection completed... \n")
  
  k_comps <- min(ncol(Z), nrow(S_star), ncol(S_star))
  cat(paste0("Performing fast partial SVD (irlba) for top ", k_comps, " components...\n"))
  
  svd_decomp <- irlba::irlba(S_star, nv = k_comps, nu = k_comps)
  cat("SVD completed... \n")
  
  U <- svd_decomp$u
  D_alpha <- diag(svd_decomp$d, nrow = length(svd_decomp$d))
  V <- svd_decomp$v
  
  gene_expression_scores <- (U %*% D_alpha)
  fitted_cell_scores <- V %*% D_alpha
  cell_scores <- crossprod(S, U)
  biplot_values <- solve(t(Z) %*% Z) %*% t(Z) %*% gene_expression_scores
  
  fitted_cell_scores <- scale(fitted_cell_scores, center = TRUE, scale = FALSE)
  cell_scores     <- scale(cell_scores, center = TRUE, scale = FALSE)
  
  fitted_cell_scores <- as.data.frame(fitted_cell_scores)
  cell_scores     <- as.data.frame(cell_scores)
  gene_expression_scores <- as.data.frame(gene_expression_scores)
  
  rownames(gene_expression_scores) <- rownames(X)
  rownames(fitted_cell_scores) <- colnames(X)
  rownames(cell_scores)     <- colnames(X)
  
  axis_names <- paste0("Axis", 1:k_comps)
  colnames(gene_expression_scores) <- axis_names
  colnames(fitted_cell_scores) <- axis_names
  colnames(cell_scores)     <- axis_names
  
  biplot_limit <- min(ncol(Z), ncol(biplot_values))
  biplot_values <- biplot_values[, 1:biplot_limit, drop = FALSE]
  colnames(biplot_values) <- paste0("Axis", 1:ncol(biplot_values))
  
  out <- list(
    expression_scores = gene_expression_scores,
    fitted_cell_scores = fitted_cell_scores,
    cell_scores = cell_scores,
    biplot = biplot_values
  )
  
  return(out)
}
