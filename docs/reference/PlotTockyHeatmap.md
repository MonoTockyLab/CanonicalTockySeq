# Plot Pseudotime Heatmap (Ordered by Peak Timing)

Visualizes gene expression along the Tocky trajectory. Supports ordering
genes by their "Peak Time" (Cascade) or Hierarchical Clustering.

## Usage

``` r
PlotTockyHeatmap(
  object,
  gradient_res,
  genes,
  n_bins = 100,
  ordering_method = c("peak", "cluster"),
  span = 0.5,
  scale_rows = TRUE
)
```

## Arguments

- object:

  Seurat object or expression matrix.

- gradient_res:

  Data frame from GradientTocky or named vector of angles.

- genes:

  Character vector of genes to plot.

- n_bins:

  Integer. Number of bins (default 100).

- ordering_method:

  Character. "peak" (sort by max expression time) or "cluster" (hclust).

- span:

  Numeric. Smoothing span (default 0.5).

- scale_rows:

  Logical. Z-score rows (default TRUE).
