# Run TradeSeq Differential Dynamics Analysis on Tocky Trajectory

A wrapper function to perform differential dynamic expression testing
using tradeSeq. Supports Seurat v5+ (uses 'layer') and matrix inputs.

## Usage

``` r
RunTradeSeqTocky(
  object,
  gradient_res,
  group_by,
  genes,
  n_knots = 5,
  n_cores = 1
)
```

## Arguments

- object:

  A Seurat object OR a sparse expression matrix (Genes x Cells).

- gradient_res:

  The data frame returned by `GradientTocky_Slerp` (must contain
  `$angle`).

- group_by:

  If `object` is Seurat: Character string of the metadata column name
  (e.g., "treatment"). If `object` is Matrix: A named factor/character
  vector of groups matching column names.

- genes:

  Character vector of genes to test (e.g., output of
  `SelectTockyGenes`).

- n_knots:

  Integer. Number of knots for the GAM (default 5).

- n_cores:

  Integer. Number of cores for parallel processing (default 1).

## Value

A list containing:

- results:

  Data frame of differential dynamics results (p-values, FDR).

- gam_fit:

  The raw tradeSeq GAM object (for downstream plotting).
