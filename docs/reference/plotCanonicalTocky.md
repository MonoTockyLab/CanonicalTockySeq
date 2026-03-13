# Plot Canonical Tocky Ordination with Continuous Gradient Legend

Visualizes the constrained ordination space in a 3-panel layout. Panel
1: Axis 1 vs 2. Panel 2: Axis 2 vs 3. Panel 3: Continuous Gradient
Legend.

## Usage

``` r
plotCanonicalTocky(tocky_res, gradient_res, alpha_level = 0.2)
```

## Arguments

- tocky_res:

  List object returned by `CanonicalTockySeq` containing `$cell_scores`
  and `$biplot`.

- gradient_res:

  List object returned by `GradientTockySeq` containing `$angle`
  Gradient Tocky Time scores (0-90).

- alpha_level:

  Numeric transparency level for cells (0-1). Default is 0.2.
