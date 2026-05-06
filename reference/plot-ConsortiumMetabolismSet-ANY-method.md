# Plot a ConsortiumMetabolismSet object

Plot a ConsortiumMetabolismSet object

## Usage

``` r
# S4 method for class 'ConsortiumMetabolismSet,ANY'
plot(
  x,
  label_colours = NULL,
  max_nodes = 20,
  label_size = 3,
  showClusterIds = TRUE
)
```

## Arguments

- x:

  A `ConsortiumMetabolismSet` object.

- label_colours:

  Optional tibble with label and colour columns.

- max_nodes:

  Maximum number of dendrogram nodes.

- label_size:

  Numeric label size.

- showClusterIds:

  Logical. If `TRUE` (default), draw the internal cluster identifiers
  used by
  [`extractCluster`](https://admarhi.github.io/ramen/reference/extractCluster.md)
  as small filled circles on top of the dendrogram. Set to `FALSE` for a
  clean dendrogram suitable for figure export.

## Value

A `ggplot` object (returned invisibly).

## Examples

``` r
cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cms <- ConsortiumMetabolismSet(cm1, cm2, name = "test")
#> 
#> ── Creating CMS "test" ─────────────────────────────────────────────────────────
#> ℹ Validating 2 <ConsortiumMetabolism> objects
#> ✔ Validating 2 <ConsortiumMetabolism> objects [12ms]
#> 
#> ℹ Collecting metabolites from 2 consortia
#> ✔ Collecting metabolites from 2 consortia [30ms]
#> 
#> ℹ Re-indexing 4 unique metabolites
#> ✔ Re-indexing 4 unique metabolites [26ms]
#> 
#> ℹ Expanding 2 binary matrices to 4-dimensional space
#> ✔ Expanding 2 binary matrices to 4-dimensional space [22ms]
#> 
#> ℹ Computing 4 x 4 levels matrix
#> ✔ Computing 4 x 4 levels matrix [23ms]
#> 
#> ℹ Computing pairwise overlap (1 pairs via crossprod)
#> ✔ Computing pairwise overlap (1 pairs via crossprod) [22ms]
#> 
#> ℹ Assembling pathway data from 2 consortia
#> ✔ Assembling pathway data from 2 consortia [30ms]
#> 
#> ℹ Building dendrogram from 2 x 2 dissimilarity matrix
#> ✔ Building dendrogram from 2 x 2 dissimilarity matrix [21ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [23ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> CMS "test" created: 2 consortia, 4 metabolites (0.2s)
#> ✔ Collecting 2 consortium graphs [88ms]
#> 
plot(cms)
```
