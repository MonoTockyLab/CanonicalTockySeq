# Gradient Tocky Sequence (Slerp Model)

Projects cells onto a piecewise Slerp manifold defined by the landmarks
B -\> BR -\> R. Cells with negative correlations to all three landmarks
are classified as "Timer Negative" (NA), as they lie in the ordination
space opposite to the Tocky trajectory.

## Usage

``` r
GradientTockySeq(res, B, BR, R, filter_negative = TRUE)
```

## Arguments

- res:

  List object returned by `CanonicalTockySeq` (must contain
  `$cell_scores` coordinates).

- B:

  Vector for the 'New' (Blue) landmark.

- BR:

  Vector for the 'Persistent' (Blue-Red) landmark.

- R:

  Vector for the 'Arrested' (Red) landmark.

- filter_negative:

  Logical. If TRUE, assigns NA to cells with negative projection to all
  three landmarks.

## Value

A data frame containing:

- angle:

  Tocky Time (0-90)

- intensity:

  Raw Euclidean norm (distance from origin)

- norm_intensity:

  Intensity normalized against the Slerp trajectory

- similarity:

  Cosine similarity to the manifold
