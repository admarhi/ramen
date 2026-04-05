# Get Functional Groups

Calculates and returns functional groups based on metabolic reactions.
For `ConsortiumMetabolismSet` objects, this involves analyzing shared
reactions across species to identify clusters of species with similar
metabolic capabilities.

## Usage

``` r
functionalGroups(object, k = 4, ...)

# S4 method for class 'ConsortiumMetabolismSet'
functionalGroups(object, k = 4, label_size = 6, label_colours = NULL)
```

## Arguments

- object:

  A `ConsortiumMetabolismSet` object.

- k:

  An integer scalar specifying the number of clusters to color in the
  dendrogram.

- ...:

  Additional arguments to be passed to specific methods.

- label_size:

  Numeric scalar specifying label size in the output plot.

- label_colours:

  If not `NULL`, a tibble with columns `label` and `colour`.

## Value

A list (returned invisibly) containing:

- plot: The ggplot2 dendrogram visualization

- dendrogram: The dendrogram object

- similarity_matrix: Matrix of Jaccard similarities between species

- species_combinations: Tibble with pairwise species comparisons

- reactions_per_species: Tibble mapping species to their reactions

The plot is automatically displayed.

## Details

This method is currently implemented for `ConsortiumMetabolismSet`
objects. Future versions will extend functionality to
`ConsortiumMetabolismAlignment` objects to allow for comparative
functional group analysis across different alignments.

## Examples

``` r
# \donttest{
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
#> ℹ Re-indexing 6 unique metabolites
#> ✔ Re-indexing 6 unique metabolites [25ms]
#> 
#> ℹ Expanding 2 binary matrices to 6-dimensional space
#> ✔ Expanding 2 binary matrices to 6-dimensional space [22ms]
#> 
#> ℹ Computing 6 x 6 levels matrix
#> ✔ Computing 6 x 6 levels matrix [23ms]
#> 
#> ℹ Computing pairwise overlap (1 pairs via crossprod)
#> ✔ Computing pairwise overlap (1 pairs via crossprod) [22ms]
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
functionalGroups(cms, k = 2)
#> Loading required namespace: colorspace
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the dendextend package.
#>   Please report the issue at <https://github.com/talgalili/dendextend/issues>.

# }
```
