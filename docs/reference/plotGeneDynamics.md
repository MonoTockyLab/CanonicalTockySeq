# Plot Gene Expression Dynamics by Group (Fixed Color Mapping & Robust Limits)

Visualizes gene expression trends along the Gradient Tocky Time axis
(0-90 degrees). Supports both data frame output from GradientTocky_Slerp
and simple numeric vectors.

## Usage

``` r
plotGeneDynamics(
  expression_data,
  gradient_res,
  gene,
  groups = NULL,
  group_cols = NULL,
  span = 0.8,
  jitter_amount = 0.2,
  pt_alpha = 0.2,
  m = 1.15
)
```

## Arguments

- expression_data:

  Sparse matrix or data frame (Genes x Cells).

- gradient_res:

  Either the data frame returned by `GradientTocky_Slerp` OR a named
  numeric vector of angles.

- gene:

  Character string. The gene symbol to plot.

- groups:

  Optional Character or Factor vector defining cell groups (e.g., "WT",
  "KO").

- group_cols:

  Optional vector of colors. Can be named (c("WT"="black")) or unnamed.

- span:

  Numeric. Smoothing span for LOESS (default 0.8).

- jitter_amount:

  Numeric. Amount of vertical jitter for points (default 0.2).

- pt_alpha:

  Numeric. Transparency of single cell points (0-1).

- m:

  Numeric. Magnification factor for y-max.
