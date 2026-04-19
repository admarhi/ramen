# Get Overlap Matrix

Returns the pairwise dissimilarity matrix from a
`ConsortiumMetabolismSet` object. Values are `1 - FOS` (Functional
Overlap Score), so 0 indicates identical consortia and 1 indicates no
shared pathways.

## Usage

``` r
overlapMatrix(object)

# S4 method for class 'ConsortiumMetabolismSet'
overlapMatrix(object)
```

## Arguments

- object:

  A `ConsortiumMetabolismSet` object.

## Value

A numeric \\n \times n\\ matrix of pairwise dissimilarities, where \\n\\
is the number of consortia. Row and column names are consortium names.

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
#> ✔ Assembling pathway data from 2 consortia [29ms]
#> 
#> ℹ Building dendrogram from 2 x 2 dissimilarity matrix
#> ✔ Building dendrogram from 2 x 2 dissimilarity matrix [20ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [22ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> CMS "test" created: 2 consortia, 6 metabolites (0.2s)
#> ✔ Collecting 2 consortium graphs [78ms]
#> 
overlapMatrix(cms)
#>           comm_1    comm_2
#> comm_1 0.0000000 0.5454545
#> comm_2 0.5454545 0.0000000
```
