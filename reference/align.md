# Align Consortium Metabolisms

Computes functional alignment between consortium metabolisms. Dispatches
on the combination of `x` and `y`:

- `align(CM, CM)`: Pairwise alignment of two consortia

- `align(CMS)`: Multiple alignment across all consortia in the set

- `align(CM, CMS)`: Database search – align one consortium against all
  in a set

## Usage

``` r
align(x, y, method = "FOS", ...)

# S4 method for class 'ConsortiumMetabolism,ConsortiumMetabolism'
align(x, y, method = "FOS", computePvalue = FALSE, nPermutations = 999L, ...)

# S4 method for class 'ConsortiumMetabolismSet,missing'
align(
  x,
  y,
  method = "FOS",
  linkage = "complete",
  BPPARAM = BiocParallel::SerialParam(),
  ...
)

# S4 method for class 'ConsortiumMetabolism,ConsortiumMetabolismSet'
align(x, y, method = "FOS", BPPARAM = BiocParallel::SerialParam(), ...)
```

## Arguments

- x:

  A
  [ConsortiumMetabolism](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
  object (query).

- y:

  A
  [ConsortiumMetabolismSet](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismSet.md)
  object (database).

- method:

  Character scalar specifying the metric.

- ...:

  Additional arguments (currently unused).

- computePvalue:

  Logical; compute permutation p-value? Default `FALSE`.

- nPermutations:

  Integer; number of permutations. Default `999L`.

- linkage:

  Character scalar specifying the agglomeration method for hierarchical
  clustering. Passed to [`hclust`](https://rdrr.io/r/stats/hclust.html).
  One of `"complete"` (default), `"average"`, `"single"`, or
  `"ward.D2"`.

- BPPARAM:

  A
  [BiocParallel::BiocParallelParam](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object. Default
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html).

## Value

A
[ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md)
object.

A
[ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md)
object of type `"pairwise"`.

A
[ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md)
object of type `"multiple"`.

A
[ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md)
object of type `"search"`.

## Functions

- `align(x = ConsortiumMetabolism, y = ConsortiumMetabolism)`: Pairwise
  alignment of two
  [ConsortiumMetabolism](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
  objects

- `align(x = ConsortiumMetabolismSet, y = missing)`: Multiple alignment
  across all consortia in a
  [ConsortiumMetabolismSet](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismSet.md)

- `align(x = ConsortiumMetabolism, y = ConsortiumMetabolismSet)`:
  Database search – align a
  [ConsortiumMetabolism](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
  against all consortia in a
  [ConsortiumMetabolismSet](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismSet.md)

## Examples

``` r
cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cma <- align(cm1, cm2)
cma
#> 
#> ── ConsortiumMetabolismAlignment 
#> Name: NA
#> Type: "pairwise"
#> Metric: "FOS"
#> Score: 0
#> Query: "comm_1", Reference: "comm_2"
```
