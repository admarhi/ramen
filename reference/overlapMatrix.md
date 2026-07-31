# Get Overlap Matrix

Returns the pairwise Functional Overlap Score (FOS) matrix from a
`ConsortiumMetabolismSet` object. Values lie in `[0, 1]`: the diagonal
is 1 (each consortium fully overlaps itself), 1 indicates identical
consortia, and 0 indicates no shared pathways.

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

A numeric \\n \times n\\ matrix of pairwise FOS overlaps, where \\n\\ is
the number of consortia. The diagonal is 1. Row and column names are
consortium names.

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
#> ✔ Validating 2 <ConsortiumMetabolism> objects [12ms]
#> 
#> ℹ Collecting metabolites from 2 consortia
#> ✔ Collecting metabolites from 2 consortia [30ms]
#> 
#> ℹ Re-indexing 6 unique metabolites
#> ✔ Re-indexing 6 unique metabolites [27ms]
#> 
#> ℹ Expanding 2 binary matrices to 6-dimensional space
#> ✔ Expanding 2 binary matrices to 6-dimensional space [24ms]
#> 
#> ℹ Computing 6 x 6 levels matrix
#> ✔ Computing 6 x 6 levels matrix [25ms]
#> 
#> ℹ Computing pairwise overlap (1 pairs via crossprod)
#> ✔ Computing pairwise overlap (1 pairs via crossprod) [23ms]
#> 
#> ℹ Assembling pathway data from 2 consortia
#> ✔ Assembling pathway data from 2 consortia [39ms]
#> 
#> ℹ Building dendrogram from 2 x 2 overlap matrix
#> ✔ Building dendrogram from 2 x 2 overlap matrix [21ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [23ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> CMS "test" created: 2 consortia, 6 metabolites (0.2s)
#> ✔ Collecting 2 consortium graphs [84ms]
#> 
overlapMatrix(cms)
#>        comm_1 comm_2
#> comm_1    1.0    0.6
#> comm_2    0.6    1.0
```
