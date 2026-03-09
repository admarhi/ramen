# Compare Alignments

Visual comparison of multiple alignments by plotting the number of
aligned reactions against the fraction of aligned communities

## Usage

``` r
compareAlignments(
  ...,
  names,
  smooth = FALSE,
  se = FALSE,
  min_frac = NULL,
  max_frac = NULL
)
```

## Arguments

- ...:

  MiCoAl objects

- names:

  Character vector supplying names to be used as labels in plot

- smooth:

  Boolean to toggle a smooth line

- se:

  Boolean to toggle std. error bands

- min_frac:

  Numerical value specifying minimum fraction in alignment

- max_frac:

  Numerical value specifying maximum fraction in alignment

## Value

A ggplot for the comparison of multiple MiCoAl objects.
