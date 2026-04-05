# Plot a ConsortiumMetabolismSet object

Plot a ConsortiumMetabolismSet object

## Usage

``` r
# S4 method for class 'ConsortiumMetabolismSet,ANY'
plot(x, label_colours = NULL, max_nodes = 20, label_size = 4)
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

## Value

A `ggplot` object (returned invisibly).

## Examples

``` r
# \donttest{
cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cms <- ConsortiumMetabolismSet(cm1, cm2, name = "test")
#> 
#> ── Creating CMS "test" ─────────────────────────────────────────────────────────
#> ℹ Validating 2 <ConsortiumMetabolism> objects
#> ✔ Validating 2 <ConsortiumMetabolism> objects [11ms]
#> 
#> ℹ Collecting metabolites from 2 consortia
#> ✔ Collecting metabolites from 2 consortia [28ms]
#> 
#> ℹ Re-indexing 6 unique metabolites
#> ✔ Re-indexing 6 unique metabolites [25ms]
#> 
#> ℹ Expanding 2 binary matrices to 6-dimensional space
#> ✔ Expanding 2 binary matrices to 6-dimensional space [21ms]
#> 
#> ℹ Computing 6 x 6 levels matrix
#> ✔ Computing 6 x 6 levels matrix [23ms]
#> 
#> ℹ Computing pairwise overlap (1 pairs via crossprod)
#> ✔ Computing pairwise overlap (1 pairs via crossprod) [21ms]
#> 
#> ℹ Assembling pathway data from 2 consortia
#> ✔ Assembling pathway data from 2 consortia [28ms]
#> 
#> ℹ Building dendrogram from 2 x 2 dissimilarity matrix
#> ✔ Building dendrogram from 2 x 2 dissimilarity matrix [20ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [22ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> ✔ CMS "test" created: 2 consortia, 6 metabolites (0.2s)
#> ℹ Collecting 2 consortium graphs
#> ✔ Collecting 2 consortium graphs [75ms]
#> 
plot(cms)

# }
```
