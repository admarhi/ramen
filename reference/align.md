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
align(
  x,
  y,
  method = "FOS",
  metrics = c("FOS", "jaccard", "brayCurtis", "redundancyOverlap"),
  topK = NULL,
  computePvalue = FALSE,
  nPermutations = 999L,
  BPPARAM = BiocParallel::SerialParam(),
  ...
)
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

  Character scalar specifying the primary metric used to rank hits. One
  of `"FOS"` (default), `"jaccard"`, `"brayCurtis"`,
  `"redundancyOverlap"`, or `"MAAS"`.

- ...:

  Additional arguments (currently unused).

- computePvalue:

  Logical; compute a permutation p-value for the top hit? Default
  `FALSE`. Not supported for `"brayCurtis"` or `"redundancyOverlap"`.

- nPermutations:

  Integer; number of permutations used when `computePvalue = TRUE`.
  Default `999L`.

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

- metrics:

  Character vector of metrics to compute per hit. Defaults to all four
  base metrics (`"FOS"`, `"jaccard"`, `"brayCurtis"`,
  `"redundancyOverlap"`). Restrict to a subset (e.g. `metrics = "FOS"`)
  to skip weighted-assay expansion and speed up large-database searches.
  `"MAAS"` as the primary metric requires all four components.

- topK:

  Integer; if non-`NULL`, truncate the ranked results (and
  `SimilarityMatrix`) to the top `topK` hits. The overall top hit's
  pathway correspondences are always reported regardless of `topK`.
  Default `NULL` (all hits).

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
object of type `"search"`. `PrimaryScore`/`ReferenceName` point to the
top hit; the full ranked hit table is stored in `Scores$ranking` and as
a 1-row `SimilarityMatrix`. Pathway correspondences (`SharedPathways`,
`UniqueQuery`, `UniqueReference`) are for the top hit only.

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

## Note

`method = "brayCurtis"` is defined only for weighted CMs (i.e.
constructed from non-unit fluxes). On unweighted inputs the score
returns `NA`. To suppress this, build CMs from a weighted edge list or
pick a different metric.

`align(CM, CMS)` (database search) currently raises a
not-yet-implemented error; the dispatch is reserved for the upcoming
MinHash-prefiltered search feature.

If `tibble` (or a package that re-exports
[`pillar::align`](https://pillar.r-lib.org/reference/align.html)) is on
the search path, plain `?align` may resolve to
[`pillar::align`](https://pillar.r-lib.org/reference/align.html). Use
`?ramen::align` to land on this page.

For the formal definitions of FOS (Szymkiewicz-Simpson), Jaccard,
Bray-Curtis, RedundancyOverlap, MAAS, coverage ratios, the Hill-1
perplexity construction underlying the EffectiveConsumption /
EffectiveProduction (flux-corrected) and nEffectiveSpeciesConsumption /
nEffectiveSpeciesProduction (effective species counts) assays, and the
"Mathematical formulation" section of
[`vignette("alignment", package = "ramen")`](https://admarhi.github.io/ramen/articles/alignment.md).
The same section discusses cross-product inflation, Hill-1 saturation on
small consortia, flux reversibility, and alternate-optima caveats
relevant to interpreting the scores returned here.

## Examples

``` r
cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cma <- align(cm1, cm2)
cma
#> 
#> ── ConsortiumMetabolismAlignment 
#> Name: "comm_1 vs comm_2"
#> Type: "pairwise"
#> Metric: "FOS"
#> Score: 0.2857
#> Query: "comm_1", Reference: "comm_2"
#> Coverage: query 0.286, reference 0.118
#> Pathways: 2 shared, 5 query-only, 15 reference-only.
```
