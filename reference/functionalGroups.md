# Get Functional Groups

Calculates and returns functional groups based on metabolic reactions.
For `ConsortiumMetabolismSet` objects, this involves analyzing shared
reactions across species to identify clusters of species with similar
metabolic capabilities.

## Usage

``` r
functionalGroups(object, ...)

# S4 method for class 'ConsortiumMetabolismSet'
functionalGroups(object, ...)
```

## Arguments

- object:

  A `ConsortiumMetabolismSet` object.

- ...:

  Additional arguments passed to methods. Supported arguments include:

  `linkage`

  :   Character scalar specifying the agglomeration method for
      hierarchical clustering. Passed to
      [`hclust`](https://rdrr.io/r/stats/hclust.html) as the `method`
      argument. One of `"complete"` (default), `"average"`, `"single"`,
      or `"ward.D2"`.

## Value

A list (returned invisibly) containing:

- `dendrogram`: The dendrogram object

- `similarity_matrix`: Matrix of Jaccard similarities between species

- `incidence_matrix`: Sparse binary species-by-pathway incidence matrix

- `reactions_per_species`: Data frame mapping species to their pathways

## Details

This method computes a Jaccard similarity matrix between species based
on shared pathways, then performs hierarchical clustering. To visualize
the resulting dendrogram, pass the output to
[`plotFunctionalGroups`](https://admarhi.github.io/ramen/reference/plotFunctionalGroups.md).

## See also

[`plotFunctionalGroups`](https://admarhi.github.io/ramen/reference/plotFunctionalGroups.md)
for visualizing the dendrogram.

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
#> ✔ Validating 2 <ConsortiumMetabolism> objects [10ms]
#> 
#> ℹ Collecting metabolites from 2 consortia
#> ✔ Collecting metabolites from 2 consortia [28ms]
#> 
#> ℹ Re-indexing 7 unique metabolites
#> ✔ Re-indexing 7 unique metabolites [24ms]
#> 
#> ℹ Expanding 2 binary matrices to 7-dimensional space
#> ✔ Expanding 2 binary matrices to 7-dimensional space [22ms]
#> 
#> ℹ Computing 7 x 7 levels matrix
#> ✔ Computing 7 x 7 levels matrix [23ms]
#> 
#> ℹ Computing pairwise overlap (1 pairs via crossprod)
#> ✔ Computing pairwise overlap (1 pairs via crossprod) [22ms]
#> 
#> ℹ Assembling pathway data from 2 consortia
#> ✔ Assembling pathway data from 2 consortia [29ms]
#> 
#> ℹ Building dendrogram from 2 x 2 dissimilarity matrix
#> ✔ Building dendrogram from 2 x 2 dissimilarity matrix [20ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [22ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> CMS "test" created: 2 consortia, 7 metabolites (0.2s)
#> ✔ Collecting 2 consortium graphs [77ms]
#> 
fg <- functionalGroups(cms)
plotFunctionalGroups(fg, k = 2)
#> Loading required namespace: colorspace
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the dendextend package.
#>   Please report the issue at <https://github.com/talgalili/dendextend/issues>.

```
