## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(CanonicalTockySeq)


## ----simulate-data------------------------------------------------------------
set.seed(2026)

n_genes <- 500
n_cells <- 300

# 1. Assign true developmental time to cells
true_time <- runif(n_cells, 0, pi/2)
cell_names <- paste0("Cell_", 1:n_cells)

# 2. Engineer Continuous Temporal Cascades for Top 100 Genes
# We use a Gaussian function so each gene peaks sequentially along the trajectory
gene_peaks <- seq(0, pi/2, length.out = 100)
X_signal <- matrix(0, nrow = 100, ncol = n_cells)
for(i in 1:100) {
  X_signal[i, ] <- exp( - (true_time - gene_peaks[i])^2 / 0.1) * 20
}

# Add 400 random noise genes to mimic background transcription
X_noise <- matrix(runif(400 * n_cells, 0, 2), nrow = 400, ncol = n_cells)
X_base <- rbind(X_signal, X_noise)

# Add statistical noise and format as non-negative counts
X <- round(X_base + matrix(rnorm(n_genes * n_cells, sd = 1), nrow = n_genes))
X[X < 0] <- 0
colnames(X) <- cell_names
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


## ----run-rda------------------------------------------------------------------
# Perform the RDA procedure
tocky_res <- CanonicalTockySeq(X = X, Z = Z)

# The result contains cell_scores and biplot vectors representing the constraints
head(tocky_res$cell_scores)


## ----run-slerp----------------------------------------------------------------
# Extract landmark vectors from the RDA biplot
landmark_B  <- tocky_res$biplot["Blue", ]
landmark_BR <- tocky_res$biplot["BlueRed", ]
landmark_R  <- tocky_res$biplot["Red", ]

# Calculate Tocky Time (Angle) and Tocky Intensity (Radial norm)
gradient_res <- GradientTockySeq(
  res = tocky_res, 
  B = landmark_B, 
  BR = landmark_BR, 
  R = landmark_R, 
  filter_negative = TRUE
)

head(gradient_res)


## ----plot-manifold, fig.width=9, fig.height=4.5-------------------------------
plotCanonicalTocky(tocky_res, gradient_res, alpha_level = 0.5)


## ----gene-dynamics------------------------------------------------------------
# Filter for top dynamic genes
dynamic_genes <- SelectTockyGenes(
  expression_data = X, 
  gradient_res = gradient_res, 
  top_n = 5
)

# Plot the expression dynamics of a top gene using LOESS smoothing
plotGeneDynamics(
  expression_data = X, 
  gradient_res = gradient_res, 
  gene = dynamic_genes[1], 
  span = 0.8
)


