## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(CanonicalTockySeq)


## ----simulate-and-build-------------------------------------------------------
set.seed(2026)
n_genes <- 500
n_cells <- 300

# 1. Assign true developmental time and groups to cells
true_time <- runif(n_cells, 0, pi/2)
cell_groups <- sample(c("WT", "KO"), n_cells, replace = TRUE)
names(cell_groups) <- paste0("Cell_", 1:n_cells)

# 2. Engineer Continuous Temporal Cascades for Top 100 Genes
# We use a Gaussian function so each gene peaks sequentially along the trajectory
gene_peaks <- seq(0, pi/2, length.out = 100)
X_signal <- matrix(0, nrow = 100, ncol = n_cells)
for(i in 1:100) {
  X_signal[i, ] <- exp( - (true_time - gene_peaks[i])^2 / 0.1) * 20
}

# Add 400 random noise genes
X_noise <- matrix(runif(400 * n_cells, 0, 2), nrow = 400, ncol = n_cells)
X_base <- rbind(X_signal, X_noise)

# Simulate a biological knockout effect: KO cells fail to fully upregulate 'Late' genes
ko_idx <- which(cell_groups == "KO")
X_base[70:100, ko_idx] <- X_base[70:100, ko_idx] * 0.5

# Format as non-negative counts for tradeSeq
X <- round(X_base + matrix(rnorm(n_genes * n_cells, sd = 1), nrow = n_genes))
X[X < 0] <- 0
colnames(X) <- names(cell_groups)
rownames(X) <- paste0("Gene_", 1:n_genes)

# 3. Create Gene Constraints (Z) matching the biological landmarks
# Blue (New) peaks at t=0, BlueRed (Persistent) at t=pi/4, Red (Arrested) at t=pi/2
all_gene_peaks <- c(gene_peaks, runif(400, 0, pi/2))
Z <- matrix(0, nrow = n_genes, ncol = 3)
colnames(Z) <- c("Blue", "BlueRed", "Red")
rownames(Z) <- rownames(X)

Z[, "Blue"]    <- exp( - (all_gene_peaks - 0)^2 / 0.2 )
Z[, "BlueRed"] <- exp( - (all_gene_peaks - pi/4)^2 / 0.2 )
Z[, "Red"]     <- exp( - (all_gene_peaks - pi/2)^2 / 0.2 )

# 4. Reconstruct Manifold
tocky_res <- CanonicalTockySeq(X = X, Z = Z)
gradient_res <- GradientTockySeq(
  res = tocky_res, 
  B = tocky_res$biplot["Blue",], 
  BR = tocky_res$biplot["BlueRed",], 
  R = tocky_res$biplot["Red",]
)

## ----pre-filter---------------------------------------------------------------
# Select the top 100 most dynamic genes for a highly impressive cascade
dynamic_genes <- SelectTockyGenes(
  expression_data = X, 
  gradient_res = gradient_res, 
  top_n = 100, 
  min_expr = 0.05
)


## ----run-tradeseq-------------------------------------------------------------
# Fit GAMs and perform statistical testing between WT and KO
tradeSeq_res <- RunTradeSeqTocky(
  object = X,
  gradient_res = gradient_res,
  group_by = cell_groups,
  genes = dynamic_genes,
  n_knots = 5,
  n_cores = 1
)

# View the differential dynamics results (p-values and FDR)
head(tradeSeq_res$results)

## ----plot-cascade, fig.width=7, fig.height=5----------------------------------
# Plot the expression cascade, ordered by the peak Tocky Time
PlotTockyHeatmap(
  object = X,
  gradient_res = gradient_res,
  genes = dynamic_genes,
  ordering_method = "peak",
  n_bins = 50,
  span = 0.5
)


