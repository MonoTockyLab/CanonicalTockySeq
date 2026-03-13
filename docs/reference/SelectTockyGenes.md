# Select Genes Dynamic along Tocky Time

Rapidly filters for genes that show significant variation along the
Gradient Tocky Time axis using a linear spline regression. Automatically
excludes cells with NA angles (Timer Negative).

## Usage

``` r
SelectTockyGenes(expression_data, gradient_res, top_n = 2000, min_expr = 0.05)
```

## Arguments

- expression_data:

  Sparse matrix or data frame (Genes x Cells).

- gradient_res:

  Either the data frame returned by `GradientTocky_Slerp` OR a named
  numeric vector of angles.

- top_n:

  Integer. Number of genes to return (default 2000).

- min_expr:

  Minimum proportion of cells expressing the gene (0-1). Default 0.05.

## Value

A character vector of gene names sorted by dynamic score (R-squared).
