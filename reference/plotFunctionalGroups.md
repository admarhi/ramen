# Plot Functional Groups Dendrogram

Visualizes the functional groups dendrogram computed by
[`functionalGroups`](https://admarhi.github.io/ramen/reference/functionalGroups.md).
Branches are colored by cluster membership and species labels are placed
below the dendrogram.

## Usage

``` r
plotFunctionalGroups(fg, k = 4, label_size = 6, label_colours = NULL)
```

## Arguments

- fg:

  A list as returned by
  [`functionalGroups`](https://admarhi.github.io/ramen/reference/functionalGroups.md),
  containing at least a `dendrogram` element.

- k:

  Integer scalar specifying the number of clusters to color in the
  dendrogram. Default `4`.

- label_size:

  Numeric scalar specifying the text size for species labels. Default
  `6`.

- label_colours:

  If not `NULL`, a data frame with columns `label` and `colour` mapping
  species names to colors. When `NULL` (default), all labels are drawn
  in black.

## Value

A `ggplot` object.

## See also

[`functionalGroups`](https://admarhi.github.io/ramen/reference/functionalGroups.md)
for computing the functional groups data.

## Examples

``` r
cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cms <- ConsortiumMetabolismSet(
    cm1, cm2, name = "test"
)
#> 
#> ── Creating CMS "test" ─────────────────────────────────────────────────────────
#> ℹ Validating 2 <ConsortiumMetabolism> objects
#> ✔ Validating 2 <ConsortiumMetabolism> objects [11ms]
#> 
#> ℹ Collecting metabolites from 2 consortia
#> ✔ Collecting metabolites from 2 consortia [31ms]
#> 
#> ℹ Re-indexing 7 unique metabolites
#> ✔ Re-indexing 7 unique metabolites [39ms]
#> 
#> ℹ Expanding 2 binary matrices to 7-dimensional space
#> ✔ Expanding 2 binary matrices to 7-dimensional space [24ms]
#> 
#> ℹ Computing 7 x 7 levels matrix
#> ✔ Computing 7 x 7 levels matrix [26ms]
#> 
#> ℹ Computing pairwise overlap (1 pairs via crossprod)
#> ✔ Computing pairwise overlap (1 pairs via crossprod) [23ms]
#> 
#> ℹ Assembling pathway data from 2 consortia
#> ✔ Assembling pathway data from 2 consortia [31ms]
#> 
#> ℹ Building dendrogram from 2 x 2 dissimilarity matrix
#> ✔ Building dendrogram from 2 x 2 dissimilarity matrix [22ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [24ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> CMS "test" created: 2 consortia, 7 metabolites (0.2s)
#> ✔ Collecting 2 consortium graphs [85ms]
#> 
fg <- functionalGroups(cms)
plotFunctionalGroups(fg, k = 2)
#> Loading required namespace: colorspace

```
