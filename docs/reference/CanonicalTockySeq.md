# Canonical Redundancy Analysis for Tocky Differentiation using single-cell RNA-sequencing data (CanonicalTockySeq)

Performs a specialized Canonical Redundancy Analysis (RDA) designed for
single-cell Tocky data. In this architecture, genes are treated as
"sites" and cells as "species", constrained by Tocky-defined
developmental stages (Blue, Blue-Red, Red).

## Usage

``` r
CanonicalTockySeq(X, Z)
```

## Arguments

- X:

  Numeric matrix of gene expression (Genes x Cells).

- Z:

  Numeric matrix of explanatory variables (Genes x Signatures).
  Typically contains Tocky-specific gene signatures (e.g., Blue,
  Blue-Red, Red).

## Value

A list containing:

- expression_scores:

  Data frame of Gene loadings (Genes x k).

- fitted_cell_scores:

  Data frame of Fitted Cell scores (Cells x k).

- cell_scores:

  Data frame of Cell scores (Cells x k). Use this for downstream
  GradientTockySeq.

- biplot:

  Matrix of biplot vectors for the constraint variables (Z).

## Details

Implements the two-step RDA procedure described by Legendre & Legendre
(2012):

1.  Multivariate regression of Expression (X) on Constraints (Z).

2.  PCA of the fitted values to derive canonical axes.
